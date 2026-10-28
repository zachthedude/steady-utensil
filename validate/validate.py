#!/usr/bin/env python3
"""
validate.py -- numerical validation of the steady-utensil control law.

Line-for-line transliteration of ../include/steady/*.h. The algorithm makes
quantitative claims ("N dB of tremor attenuation"), so they get measured.

Structure:
  PROBE 1  inner joint loop step response        -- is the PID even stable?
  PROBE 2  WFLC lock and prediction accuracy     -- open loop, no plant
  PROBE 3  full closed loop vs. two naive designs
  PROBE 4  predict-ahead sweep                   -- the phase budget

Plant includes the things that actually break tremor cancellers:
  one control period of transport delay, finite authority, encoder
  quantisation, velocity-estimate noise, gyro bias and noise, an
  accelerometer corrupted by the tremor's own linear acceleration, base
  acceleration coupling, Coulomb friction, and motor cogging torque.

Run:  python validate.py
"""

import math
import random

TWO_PI = 2.0 * math.pi
G = 9.80665


def clampf(v, lo, hi):
    return lo if v < lo else (hi if v > hi else v)


def wrap_pi(a):
    while a > math.pi:
        a -= TWO_PI
    while a < -math.pi:
        a += TWO_PI
    return a


# ---------------------------------------------------------------------------
# Controller components (mirrors of the C++)
# ---------------------------------------------------------------------------
class Biquad:
    def __init__(self):
        self.b0, self.b1, self.b2, self.a1, self.a2 = 1.0, 0.0, 0.0, 0.0, 0.0
        self.z1 = self.z2 = 0.0

    def _pw(self, fc, fs, q):
        w0 = TWO_PI * fc / fs
        return math.cos(w0), math.sin(w0) / (2.0 * q)

    def lowpass(self, fc, fs, q=0.70710678):
        cw, al = self._pw(fc, fs, q)
        a0 = 1.0 + al
        self.b0 = ((1.0 - cw) * 0.5) / a0
        self.b1 = (1.0 - cw) / a0
        self.b2 = self.b0
        self.a1 = (-2.0 * cw) / a0
        self.a2 = (1.0 - al) / a0
        return self

    def highpass(self, fc, fs, q=0.70710678):
        cw, al = self._pw(fc, fs, q)
        a0 = 1.0 + al
        self.b0 = ((1.0 + cw) * 0.5) / a0
        self.b1 = -(1.0 + cw) / a0
        self.b2 = self.b0
        self.a1 = (-2.0 * cw) / a0
        self.a2 = (1.0 - al) / a0
        return self

    def process(self, x):
        y = self.b0 * x + self.z1
        self.z1 = self.b1 * x - self.a1 * y + self.z2
        self.z2 = self.b2 * x - self.a2 * y
        return y

    def reset(self, x=0.0):
        den = 1.0 + self.a1 + self.a2
        g = (self.b0 + self.b1 + self.b2) / den if abs(den) > 1e-12 else 0.0
        y = g * x
        self.z1 = y - self.b0 * x
        self.z2 = self.b2 * x - self.a2 * y
        return self


class TiltEstimator:
    ACCEL_GATE = 1.5
    MAX_BIAS = 0.15

    def __init__(self, crossover_hz, zeta=1.0):
        wc = TWO_PI * crossover_hz
        self.k_angle = 2.0 * zeta * wc
        self.k_bias = wc * wc
        self.angle = self.bias = 0.0
        self.gated = False

    def reset(self, angle):
        self.angle, self.bias = angle, 0.0

    def update(self, accel_angle, accel_mag, gyro_rate, dt):
        rate = gyro_rate - self.bias
        self.angle += rate * dt
        if abs(accel_mag - G) < self.ACCEL_GATE:
            e = wrap_pi(accel_angle - self.angle)
            self.angle += self.k_angle * e * dt
            self.bias -= self.k_bias * e * dt
            self.bias = clampf(self.bias, -self.MAX_BIAS, self.MAX_BIAS)
            self.gated = False
        else:
            self.gated = True
        self.angle = wrap_pi(self.angle)
        return rate


class Wflc:
    """Rate-form adaptation so the tuning is sample-rate independent:
       mu_w    [1/s]      weight adaptation
       mu_omega[rad/s^2]  frequency adaptation (amplitude-normalised)
    """
    EPS = 1e-8

    def __init__(self, M=2, mu_w=6.0, mu_omega=250.0,
                 f_init=6.0, f_min=3.0, f_max=14.0, freq_harmonics=1):
        self.M, self.mu_w, self.mu_omega = M, mu_w, mu_omega
        self.Mf = min(freq_harmonics, M)
        self.f_min, self.f_max = f_min, f_max
        self.omega = TWO_PI * f_init
        self.f_init = f_init
        self.w = [0.0] * (2 * M)
        self.phase = 0.0

    def reset(self):
        self.w = [0.0] * (2 * self.M)
        self.phase = 0.0
        self.omega = TWO_PI * self.f_init

    def update(self, s, dt):
        M, ph = self.M, self.phase
        x = [0.0] * (2 * M)
        for r in range(1, M + 1):
            rp = r * ph
            x[r - 1] = math.sin(rp)
            x[M + r - 1] = math.cos(rp)

        y = sum(self.w[i] * x[i] for i in range(2 * M))
        e = s - y

        # Frequency gradient from the FUNDAMENTAL ONLY. Including harmonic r
        # weights the gradient by r, so a 25% second harmonic contributes half
        # the frequency drive of the fundamental it is supposed to be a
        # detail of -- enough to pull omega off lock entirely. Harmonics are
        # still reconstructed by the weight vector; they just do not get a vote
        # on what the frequency is.
        grad = 0.0
        for r in range(1, self.Mf + 1):
            rp = r * ph
            grad += r * (self.w[r - 1] * math.cos(rp) -
                         self.w[M + r - 1] * math.sin(rp))

        amp2 = sum(self.w[r] ** 2 + self.w[M + r] ** 2
                   for r in range(self.Mf))
        self.omega += 2.0 * self.mu_omega * e * grad / (amp2 + self.EPS) * dt
        self.omega = clampf(self.omega, TWO_PI * self.f_min, TWO_PI * self.f_max)

        g = 2.0 * self.mu_w * e * dt
        for i in range(2 * M):
            self.w[i] += g * x[i]

        self.phase = wrap_pi(ph + self.omega * dt)
        return y

    def _eval(self, t_ahead, order):
        M = self.M
        ph = self.phase + self.omega * t_ahead
        out = 0.0
        for r in range(1, M + 1):
            rp = r * ph
            s_, c_ = math.sin(rp), math.cos(rp)
            a, b = self.w[r - 1], self.w[M + r - 1]
            if order == 0:
                out += a * s_ + b * c_
            elif order == 1:
                out += r * self.omega * (a * c_ - b * s_)
            else:
                out -= r * r * self.omega * self.omega * (a * s_ + b * c_)
        return out

    def predict(self, t):        return self._eval(t, 0)
    def predict_rate(self, t):   return self._eval(t, 1)
    def predict_accel(self, t):  return self._eval(t, 2)
    def frequency_hz(self):      return self.omega / TWO_PI
    def amplitude(self):         return math.hypot(self.w[0], self.w[self.M])


class Pid:
    def __init__(self, kp=0.0, ki=0.0, kd=0.0, b=1.0, c=0.0,
                 out_min=-1.0, out_max=1.0, slew=1e9,
                 i_min=-0.20, i_max=0.20, tracking_tc=0.05):
        self.kp, self.ki, self.kd, self.b, self.c = kp, ki, kd, b, c
        self.out_min, self.out_max, self.slew = out_min, out_max, slew
        self.i_min, self.i_max, self.tt = i_min, i_max, tracking_tc
        self.integ = self.out = 0.0
        self.seeded = self.sat = False

    def reset(self, output_now=0.0):
        self.integ = 0.0
        self.out = clampf(output_now, self.out_min, self.out_max)
        self.seeded = self.sat = False

    def update(self, sp, sp_rate, meas, rate, ff, dt):
        dt = dt if dt > 1e-6 else 1e-6
        err = sp - meas
        p = self.kp * (self.b * sp - meas)
        d = self.kd * (self.c * sp_rate - rate)

        if not self.seeded:
            self.integ = clampf(self.out - (p + d + ff), self.i_min, self.i_max)
            self.seeded = True

        u_raw = p + self.integ + d + ff
        u = clampf(u_raw, self.out_min, self.out_max)
        step = self.slew * dt
        u = clampf(u, self.out - step, self.out + step)

        tt = self.tt if self.tt > 1e-4 else 1e-4
        self.integ += (self.ki * err + (u - u_raw) / tt) * dt
        self.integ = clampf(self.integ, self.i_min, self.i_max)

        self.sat = abs(u - u_raw) > 1e-6
        self.out = u
        return u


class AxisStabilizer:
    """mode: 'full' | 'reactive' | 'level'"""

    def __init__(self, cfg, mode="full"):
        self.c, self.mode = cfg, mode
        fs = cfg["fs"]
        self.intent_lp = Biquad().lowpass(cfg["intent_hz"], fs)
        self.band_hi = Biquad().highpass(cfg["tremor_lo_hz"], fs)
        self.band_lo = Biquad().lowpass(cfg["tremor_hi_hz"], fs)
        self.wflc = Wflc(**cfg["wflc"])
        self.pid = Pid(**cfg["pid"])
        self.q_des = 0.0

    def reset(self, handle_angle, joint_angle):
        self.intent_lp.reset(handle_angle)
        self.band_hi.reset(0.0)
        self.band_lo.reset(0.0)
        self.wflc.reset()
        self.pid.reset(0.0)
        self.q_des = joint_angle

    def update(self, handle_angle, joint_angle, joint_rate,
               gravity_ff, base_accel, dt, enable):
        c = self.c
        intent = self.intent_lp.process(handle_angle)
        s = self.band_lo.process(self.band_hi.process(handle_angle - intent))

        model = self.wflc.update(s, dt)
        resid = s - model
        lat = c["latency_s"]

        if not enable:
            q_t = q_td = q_tdd = 0.0
        elif self.mode == "full":
            q_t = -(self.wflc.predict(lat) + c["residual_gain"] * resid)
            q_td = -self.wflc.predict_rate(lat)
            q_tdd = -self.wflc.predict_accel(lat)
        elif self.mode == "reactive":
            # "Detect the shake, oppose the shake": no prediction, no model.
            q_t, q_td, q_tdd = -s, 0.0, 0.0
        elif self.mode == "level":
            # "Keep it level in the world frame": also cancels intent.
            q_t, q_td, q_tdd = -handle_angle, 0.0, 0.0

        q_t = clampf(q_t, -c["travel_rad"], c["travel_rad"])
        step = c["max_joint_rate"] * dt
        q_prev = self.q_des
        self.q_des = clampf(q_t, q_prev - step, q_prev + step)
        if self.q_des != q_t:
            q_td, q_tdd = (self.q_des - q_prev) / dt, 0.0

        # World-frame acceleration feedforward. What the motor must supply is
        # set by the head's acceleration in the WORLD, not in the handle frame,
        # and (q_ddot + base_accel) is exactly that. For a direct-drive gimbal
        # this is ~0 in the tremor band, which is the mathematical statement of
        # the fact that a low-friction gimbal is already passively isolating.
        ff = gravity_ff + c["accel_ff_gain"] * (q_tdd + base_accel) if enable else 0.0
        return self.pid.update(self.q_des, q_td, joint_angle, joint_rate, ff, dt)


# ---------------------------------------------------------------------------
# Plant
# ---------------------------------------------------------------------------
J_JOINT = 6.0e-5      # [kg m^2]
KT = 0.060            # [N m per unit command]
B_VISC = 2.0e-4       # [N m s/rad]
TAU_COUL = 1.5e-3     # [N m]  bearing + seal Coulomb friction
TAU_COG = 6.0e-4      # [N m]  motor cogging amplitude
COG_CYCLES = 84.0     # cogging cycles per revolution
M_HEAD = 0.035        # [kg]   spoon head, fully loaded
L_HEAD = 0.055        # [m]
TRAVEL = 0.44         # [rad]  +/- 25 deg
L_ARM = 0.040         # [m]    tremor rotation radius at the wrist

ENC_LSB = TWO_PI / 16384.0
QD_NOISE = 0.05
GYRO_NOISE = 0.010
GYRO_BIAS = 0.020
ACC_NOISE = 0.15

FS = 500.0
DT = 1.0 / FS
N_SUB = 10
DT_P = DT / N_SUB

ACCEL_FF_GAIN = J_JOINT / KT             # [cmd per rad/s^2]
GRAV_FF_K = M_HEAD * G * L_HEAD / KT     # [cmd] to hold a loaded head level

PID_KW = dict(kp=19.0, ki=120.0, kd=0.22, b=1.0, c=1.0,
              out_min=-1.0, out_max=1.0, slew=200.0,
              i_min=-0.35, i_max=0.35, tracking_tc=0.02)


def joint_torque(u, q, qd, th, thdd):
    return (KT * u
            - B_VISC * qd
            - TAU_COUL * math.tanh(qd / 0.02)
            - TAU_COG * math.sin(COG_CYCLES * q)
            - M_HEAD * G * L_HEAD * math.cos(th + q)
            - J_JOINT * thdd)                       # base-acceleration coupling


def handle_angle(t, f0=5.0, f1=9.0, dur=30.0, a_tremor=0.070):
    voluntary = (0.35 * math.sin(TWO_PI * 0.25 * t) +
                 0.15 * math.sin(TWO_PI * 0.11 * t + 0.7))
    k = (f1 - f0) / dur
    phi = TWO_PI * (f0 * t + 0.5 * k * t * t)
    amp = a_tremor * (1.0 + 0.30 * math.sin(TWO_PI * 0.13 * t))
    tremor = amp * (math.sin(phi) + 0.25 * math.sin(2.0 * phi + 0.6))
    return voluntary + tremor, voluntary, tremor


def handle_derivs(t, **kw):
    h = 2.0e-4
    a0 = handle_angle(t - h, **kw)[0]
    a1, v1, tr1 = handle_angle(t, **kw)
    a2 = handle_angle(t + h, **kw)[0]
    return a1, v1, tr1, (a2 - a0) / (2 * h), (a2 - 2 * a1 + a0) / (h * h)


def true_freq(t, f0=5.0, f1=9.0, dur=30.0):
    return f0 + (f1 - f0) / dur * t


def make_cfg(latency_s=-0.003, mode="full", accel_ff=0.0,
             lo=0.8, hi=40.0, mu_w=8.0, mu_omega=200.0, residual_gain=0.25):
    return {
        "fs": FS, "intent_hz": 1.2,
        "tremor_lo_hz": lo, "tremor_hi_hz": hi,
        "latency_s": latency_s, "travel_rad": TRAVEL,
        "max_joint_rate": 12.0, "residual_gain": residual_gain,
        "accel_ff_gain": accel_ff,
        "wflc": dict(M=2, mu_w=mu_w, mu_omega=mu_omega, freq_harmonics=1,
                     f_init=6.0, f_min=3.0, f_max=14.0),
        "pid": dict(PID_KW),
    }


# ---------------------------------------------------------------------------
# PROBE 1 -- inner joint loop step response
# ---------------------------------------------------------------------------
def probe_inner_loop(step=math.radians(5.0), dur=0.5):
    pid = Pid(**PID_KW)
    pid.reset(0.0)
    q = qd = 0.0
    prev_cmd = 0.0
    n = int(dur * FS)
    hist = []
    for k in range(n):
        u_new = pid.update(step, 0.0, q, qd, GRAV_FF_K * math.cos(q), DT)
        u_applied, prev_cmd = prev_cmd, u_new
        for _ in range(N_SUB):
            qd += (joint_torque(u_applied, q, qd, 0.0, 0.0) / J_JOINT) * DT_P
            q += qd * DT_P
        hist.append(q)
    peak = max(hist)
    overshoot = 100.0 * (peak - step) / step
    settle = None
    for i in range(n - 1, -1, -1):
        if abs(hist[i] - step) > 0.02 * step:
            settle = (i + 1) / FS
            break
    final_err = math.degrees(hist[-1] - step)
    return overshoot, settle, final_err, max(abs(h) for h in hist) < 10 * step


# ---------------------------------------------------------------------------
# PROBE 2 -- WFLC lock and prediction accuracy, no plant in the loop
# ---------------------------------------------------------------------------
def probe_wflc(latency_s=-0.003, dur=30.0, warmup=6.0, lo=0.8, hi=40.0,
               mu_w=8.0, mu_omega=200.0):
    intent_lp = Biquad().lowpass(1.2, FS)
    bhi = Biquad().highpass(lo, FS)
    blo = Biquad().lowpass(hi, FS)
    w = Wflc(M=2, mu_w=mu_w, mu_omega=mu_omega)
    th0 = handle_angle(0.0)[0]
    intent_lp.reset(th0)

    n = int(dur * FS)
    i0 = int(warmup * FS)
    ferr, num, den = [], 0.0, 0.0
    for k in range(n):
        t = k * DT
        th = handle_angle(t)[0]
        s = blo.process(bhi.process(th - intent_lp.process(th)))
        w.update(s, DT)
        if k >= i0:
            pred = w.predict(latency_s)
            truth = handle_angle(t + latency_s)[2]      # true tremor at t+lat
            num += (pred - truth) ** 2
            den += truth ** 2
            ferr.append(abs(w.frequency_hz() - true_freq(t)))
    resid_ratio = math.sqrt(num / den)
    return resid_ratio, (sum(ferr) / len(ferr)), 20.0 * math.log10(resid_ratio)


# ---------------------------------------------------------------------------
# PROBE 3 -- full closed loop
# ---------------------------------------------------------------------------
def simulate(mode="full", latency_s=-0.003, accel_ff=0.0, dur=30.0, seed=12345,
             lo=0.8, hi=40.0, mu_w=8.0, mu_omega=200.0, residual_gain=0.25,
             extra_delay=0, sense_delay=0, pid_kw=None):
    rng = random.Random(seed)
    cfg = make_cfg(latency_s, mode, accel_ff, lo, hi, mu_w, mu_omega,
                   residual_gain)
    if pid_kw:
        cfg["pid"] = dict(PID_KW, **pid_kw)
    axis = AxisStabilizer(cfg, mode)
    tilt = TiltEstimator(0.30, 1.0)
    dgyro = Biquad().lowpass(30.0, FS)

    th0 = handle_derivs(0.0)[0]
    tilt.reset(th0)
    axis.reset(th0, 0.0)

    q = qd = 0.0
    fifo = [0.0] * (1 + extra_delay)     # 1 = the baseline one-tick delay
    sfifo = None
    prev_rate = 0.0
    n = int(dur * FS)
    rec = dict(th=[], ts=[], vol=[], f=[], u=[])
    n_sat = gate = 0

    for k in range(n):
        t = k * DT
        th, vol, _tr, thd, thdd = handle_derivs(t)

        gyro = thd + GYRO_BIAS + rng.gauss(0.0, GYRO_NOISE)
        a_lin = L_ARM * thdd
        acc_angle = th + a_lin / G + rng.gauss(0.0, ACC_NOISE) / G
        acc_mag = math.hypot(G, a_lin) + rng.gauss(0.0, ACC_NOISE)
        rate = tilt.update(acc_angle, acc_mag, gyro, DT)
        if tilt.gated:
            gate += 1
        th_hat = tilt.angle
        if sense_delay:                      # sensor-pipeline delay only:
            if sfifo is None:                # the joint servo loop is
                sfifo = [th_hat] * sense_delay   # left untouched
            sfifo.append(th_hat)
            th_hat = sfifo.pop(0)
        base_accel = dgyro.process((rate - prev_rate) / DT)
        prev_rate = rate

        q_meas = round(q / ENC_LSB) * ENC_LSB
        qd_meas = qd + rng.gauss(0.0, QD_NOISE)

        grav_ff = GRAV_FF_K * math.cos(th_hat + q_meas)
        u_new = axis.update(th_hat, q_meas, qd_meas, grav_ff, base_accel,
                            DT, True)
        fifo.append(u_new)
        u_applied = fifo.pop(0)
        if axis.pid.sat:
            n_sat += 1

        for _ in range(N_SUB):
            qd += (joint_torque(u_applied, q, qd, th, thdd) / J_JOINT) * DT_P
            q += qd * DT_P
            if q > TRAVEL:
                q, qd = TRAVEL, min(qd, 0.0)
            elif q < -TRAVEL:
                q, qd = -TRAVEL, max(qd, 0.0)

        rec["th"].append(th)
        rec["ts"].append(th + q)
        rec["vol"].append(vol)
        rec["f"].append(axis.wflc.frequency_hz())
        rec["u"].append(u_applied)

    rec["sat_frac"] = n_sat / n
    rec["gate_frac"] = gate / n
    return rec


def filtfilt(make, x):
    def once(sig):
        bq = make(); bq.reset(sig[0])
        return [bq.process(v) for v in sig]
    return list(reversed(once(list(reversed(once(x))))))


def bandpass_series(x, lo=3.0, hi=14.0):
    return filtfilt(lambda: Biquad().lowpass(hi, FS),
                    filtfilt(lambda: Biquad().highpass(lo, FS), x))


def rms(x):
    return math.sqrt(sum(v * v for v in x) / len(x)) if x else 0.0


def analyse(r, warmup=8.0, dur=30.0):
    i0, i1 = int(warmup * FS), int(dur * FS)
    th_bp = bandpass_series(r["th"])[i0:i1]
    ts_bp = bandpass_series(r["ts"])[i0:i1]
    ri, ro = rms(th_bp), rms(ts_bp)
    att = 20.0 * math.log10(ro / ri) if ri > 0 and ro > 0 else float("nan")
    ts_lp = filtfilt(lambda: Biquad().lowpass(2.0, FS), r["ts"])[i0:i1]
    vol_err = rms([a - b for a, b in zip(ts_lp, r["vol"][i0:i1])])
    ferr = [abs(r["f"][i] - true_freq(i / FS)) for i in range(i0, i1)]
    return dict(in_deg=math.degrees(ri), out_deg=math.degrees(ro), att_db=att,
                reduction_pct=100.0 * (1.0 - ro / ri),
                vol_err_deg=math.degrees(vol_err),
                f_err_hz=sum(ferr) / len(ferr),
                peak_u=max(abs(v) for v in r["u"][i0:i1]),
                sat_frac=r["sat_frac"], gate_frac=r["gate_frac"])


def main():
    print("=" * 76)
    print("PROBE 1 -- inner joint loop, 5 deg step with a loaded head")
    print("=" * 76)
    ov, ts, fe, stable = probe_inner_loop()
    print(f"  overshoot        {ov:6.1f} %")
    print(f"  settle to 2%     {ts*1000 if ts else 0:6.1f} ms")
    print(f"  steady-state err {fe:6.3f} deg")
    print(f"  bounded          {stable}")

    print()
    print("=" * 76)
    print("PROBE 2 -- WFLC lock + prediction, open loop (best case achievable)")
    print("=" * 76)
    print(f"{'predict-ahead [ms]':>20}{'resid/signal':>16}{'ideal dB':>12}"
          f"{'freq err [Hz]':>16}")
    best = (None, 9e9, None)
    for ms in range(-10, 9, 2):
        rr, fh, db = probe_wflc(ms / 1000.0)
        print(f"{ms:>20d}{rr:>16.3f}{db:>12.1f}{fh:>16.3f}")
        if db < best[1]:
            best = (ms, db, fh)
    print(f"\n  best open-loop predict-ahead: {best[0]} ms -> {best[1]:.1f} dB")

    # Report the SHIPPED default, not probe 2's open-loop optimum -- otherwise
    # this table silently measures a config that tuning.h does not use.
    lat = make_cfg()["latency_s"]
    print(f"  (open-loop optimum {best[0]:+d} ms; shipped default "
          f"{lat*1000:+.0f} ms -- closed loop below uses the shipped value)")

    print()
    print("=" * 76)
    print("PROBE 3 -- full closed loop, 30 s, 5.0->9.0 Hz chirped tremor")
    print("=" * 76)
    rows = []
    for mode, ff, label in [
            ("full", 0.0, "Full design (WFLC predict + PID)"),
            ("reactive", 0.0, "Reactive PID only (no prediction)"),
            ("level", 0.0, "Naive hold-level in world frame")]:
        rows.append((label, analyse(simulate(mode=mode, latency_s=lat,
                                             accel_ff=ff))))
    hdr = (f"{'controller':<38}{'in deg':>9}{'out deg':>9}{'dB':>8}"
           f"{'volerr deg':>12}")
    print(hdr)
    print("-" * len(hdr))
    for label, a in rows:
        print(f"{label:<38}{a['in_deg']:>9.3f}{a['out_deg']:>9.3f}"
              f"{a['att_db']:>8.1f}{a['vol_err_deg']:>12.2f}")

    a = rows[0][1]
    print(f"\n  tremor reduction        {a['reduction_pct']:.1f} %")
    print(f"  WFLC mean freq error    {a['f_err_hz']:.3f} Hz")
    print(f"  peak actuator command   {a['peak_u']:.3f} of 1.0")
    print(f"  samples PID saturated   {100*a['sat_frac']:.2f} %")
    print(f"  accel samples rejected  {100*a['gate_frac']:.1f} %")

    print()
    print("=" * 76)
    print("PROBE 4 -- closed-loop predict-ahead sweep (the phase budget)")
    print("=" * 76)
    print(f"{'predict-ahead [ms]':>20}{'attenuation dB':>18}{'reduction %':>15}")
    for ms in range(-10, 9, 2):
        a = analyse(simulate(mode="full", latency_s=ms / 1000.0))
        print(f"{ms:>20d}{a['att_db']:>18.1f}{a['reduction_pct']:>15.1f}")
    print("\n  Negative dB = attenuation. Positive dB = the device is adding")
    print("  energy at tremor frequency, i.e. making the shaking worse.")


if __name__ == "__main__":
    main()
