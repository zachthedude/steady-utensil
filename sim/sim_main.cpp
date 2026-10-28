// SPDX-License-Identifier: MIT
// sim/sim_main.cpp -- desktop closed-loop simulation of the stabilised utensil.
//
// C++ counterpart of validate/validate.py. Same plant, same disturbance, same
// metrics, so the two should agree to within float-vs-double rounding. Its job
// is to let you re-run the tuning sweeps against the actual shipping headers
// after you change them.
//
//   g++ -std=c++17 -O2 -I../include sim_main.cpp -o sim && ./sim
//
// Reports attenuation inside the 3-14 Hz band, plus how hard the actuator had
// to work to get it.
// ---------------------------------------------------------------------------

#include "steady/tuning.h"

#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

using namespace steady;

namespace {

// --- plant constants (mirror validate.py) --------------------------------
constexpr float kJ        = 6.0e-5f;    // [kg m^2]
constexpr float kKt       = 0.060f;     // [N m per unit command]
constexpr float kBVisc    = 2.0e-4f;    // [N m s/rad]
constexpr float kTauCoul  = 1.5e-3f;    // [N m] bearing Coulomb friction
constexpr float kTauCog   = 6.0e-4f;    // [N m] cogging amplitude
constexpr float kCogCyc   = 84.0f;      // cogging cycles per revolution
constexpr float kMHead    = 0.035f;     // [kg]
constexpr float kLHead    = 0.055f;     // [m]
constexpr float kTravel   = 0.44f;      // [rad]
constexpr float kLArm     = 0.040f;     // [m] tremor rotation radius at wrist

constexpr float kEncLsb   = kTwoPi / 16384.0f;   // 14-bit encoder
constexpr float kQdNoise  = 0.05f;      // [rad/s]
constexpr float kGyroNoise= 0.010f;     // [rad/s]
constexpr float kGyroBias = 0.020f;     // [rad/s]
constexpr float kAccNoise = 0.15f;      // [m/s^2]

constexpr float kFs   = 500.0f;
constexpr float kDt   = 1.0f / kFs;
constexpr int   kSub  = 10;
constexpr float kDtP  = kDt / kSub;

constexpr float kDur    = 30.0f;
constexpr float kWarmup = 8.0f;
constexpr float kF0 = 5.0f, kF1 = 9.0f;

// --- disturbance ---------------------------------------------------------
struct Handle { float angle, voluntary, tremor; };

Handle handleAt(double t) {
    const double vol = 0.35 * std::sin(kTwoPi * 0.25 * t) +
                       0.15 * std::sin(kTwoPi * 0.11 * t + 0.7);
    const double k   = (kF1 - kF0) / kDur;
    const double phi = kTwoPi * (kF0 * t + 0.5 * k * t * t);
    const double amp = 0.070 * (1.0 + 0.30 * std::sin(kTwoPi * 0.13 * t));
    const double tr  = amp * (std::sin(phi) + 0.25 * std::sin(2.0 * phi + 0.6));
    return { float(vol + tr), float(vol), float(tr) };
}

struct Derivs { float angle, voluntary, rate, accel; };

Derivs handleDerivs(double t) {
    constexpr double h = 2.0e-4;
    const double a0 = handleAt(t - h).angle;
    const Handle m  = handleAt(t);
    const double a2 = handleAt(t + h).angle;
    return { m.angle, m.voluntary,
             float((a2 - a0) / (2.0 * h)),
             float((a2 - 2.0 * m.angle + a0) / (h * h)) };
}

float jointTorque(float u, float q, float qd, float th, float thdd) {
    return kKt * u
         - kBVisc * qd
         - kTauCoul * std::tanh(qd / 0.02f)
         - kTauCog  * std::sin(kCogCyc * q)
         - kMHead * kG * kLHead * std::cos(th + q)
         - kJ * thdd;                       // base-acceleration coupling
}

// --- offline zero-phase analysis ----------------------------------------
enum class Kind { kLow, kHigh };

std::vector<float> pass(const std::vector<float>& x, Kind kind, float fc) {
    auto once = [&](const std::vector<float>& in) {
        Biquad b;
        if (kind == Kind::kLow) b.lowpass(fc, kFs); else b.highpass(fc, kFs);
        b.reset(in.front());
        std::vector<float> out(in.size());
        for (size_t i = 0; i < in.size(); ++i) out[i] = b.process(in[i]);
        return out;
    };
    std::vector<float> y = once(x);
    std::vector<float> r(y.rbegin(), y.rend());
    y = once(r);
    return std::vector<float>(y.rbegin(), y.rend());
}

std::vector<float> bandpass(const std::vector<float>& x) {
    return pass(pass(x, Kind::kHigh, 3.0f), Kind::kLow, 14.0f);
}

double rms(const std::vector<float>& x, size_t i0, size_t i1) {
    double s = 0.0;
    for (size_t i = i0; i < i1; ++i) s += double(x[i]) * x[i];
    return std::sqrt(s / double(i1 - i0));
}

// --- run -----------------------------------------------------------------
struct Result {
    double in_deg, out_deg, att_db, reduction, vol_err_deg;
    double peak_u, sat_pct, gate_pct, f_err_hz;
};

Result run(float latency_s, bool enabled = true) {
    AxisStabilizer::Config cfg = reference::axisConfig();
    cfg.latency_s = latency_s;

    AxisStabilizer axis;
    axis.configure(cfg);

    TiltEstimator tilt;
    tilt.configure(reference::kTiltCrossoverHz, reference::kTiltZeta);

    Biquad daccel;
    daccel.lowpass(reference::kBaseAccelLpHz, kFs);

    const Derivs d0 = handleDerivs(0.0);
    tilt.reset(d0.angle);
    axis.reset(d0.angle, 0.0f);

    std::mt19937 rng(12345u);
    std::normal_distribution<float> ng(0.0f, kGyroNoise), na(0.0f, kAccNoise),
                                    nq(0.0f, kQdNoise);

    const int n = int(kDur * kFs);
    std::vector<float> th(n), ts(n), vol(n);
    float q = 0.0f, qd = 0.0f, prev_cmd = 0.0f, prev_rate = 0.0f;
    double peak_u = 0.0, f_err = 0.0;
    int n_sat = 0, n_gate = 0;
    const size_t i0 = size_t(kWarmup * kFs), i1 = size_t(n);

    for (int k = 0; k < n; ++k) {
        const double t = k * double(kDt);
        const Derivs h = handleDerivs(t);

        const float gyro    = h.rate + kGyroBias + ng(rng);
        const float a_lin   = kLArm * h.accel;
        const float acc_ang = h.angle + a_lin / kG + na(rng) / kG;
        const float acc_mag = std::sqrt(kG * kG + a_lin * a_lin) + na(rng);

        const float rate = tilt.update(acc_ang, acc_mag, gyro, kDt);
        if (tilt.accelGated()) ++n_gate;
        const float base_accel = daccel.process((rate - prev_rate) / kDt);
        prev_rate = rate;

        const float q_meas  = std::round(q / kEncLsb) * kEncLsb;
        const float qd_meas = qd + nq(rng);
        const float grav_ff = reference::kLoadFeedforward *
                              std::cos(tilt.angle() + q_meas);

        AxisStabilizer::Telemetry tm;
        const float u_new = axis.update(tilt.angle(), q_meas, qd_meas, grav_ff,
                                        base_accel, kDt,
                                        /*enable=*/enabled, &tm);
        const float u = prev_cmd;          // one control period of transport delay
        prev_cmd = u_new;
        if (tm.saturated) ++n_sat;

        for (int i = 0; i < kSub; ++i) {
            qd += (jointTorque(u, q, qd, h.angle, h.accel) / kJ) * kDtP;
            q  += qd * kDtP;
            if (q >  kTravel) { q =  kTravel; qd = std::fmin(qd, 0.0f); }
            if (q < -kTravel) { q = -kTravel; qd = std::fmax(qd, 0.0f); }
        }

        th[k]  = h.angle;
        ts[k]  = h.angle + q;              // true utensil attitude in the world
        vol[k] = h.voluntary;
        if (size_t(k) >= i0) {
            peak_u = std::fmax(peak_u, std::fabs(double(u)));
            f_err += std::fabs(tm.tremor_hz - (kF0 + (kF1 - kF0) / kDur * t));
        }
    }

    const std::vector<float> th_bp = bandpass(th), ts_bp = bandpass(ts);
    const double ri = rms(th_bp, i0, i1), ro = rms(ts_bp, i0, i1);

    const std::vector<float> ts_lp = pass(ts, Kind::kLow, 2.0f);
    std::vector<float> verr(n);
    for (int k = 0; k < n; ++k) verr[k] = ts_lp[k] - vol[k];

    Result r;
    r.in_deg      = ri * 180.0 / kPi;
    r.out_deg     = ro * 180.0 / kPi;
    r.att_db      = 20.0 * std::log10(ro / ri);
    r.reduction   = 100.0 * (1.0 - ro / ri);
    r.vol_err_deg = rms(verr, i0, i1) * 180.0 / kPi;
    r.peak_u      = peak_u;
    r.sat_pct     = 100.0 * n_sat / n;
    r.gate_pct    = 100.0 * n_gate / n;
    r.f_err_hz    = f_err / double(i1 - i0);
    return r;
}

} // namespace

int main() {
    std::printf("Stabilised utensil -- closed-loop simulation\n");
    std::printf("30 s, 5.0->9.0 Hz chirped two-harmonic tremor over a 0.25 Hz reach\n");
    std::printf("Metrics over t = 8..30 s, attenuation measured in 3-14 Hz\n\n");

    const Result r = run(reference::axisConfig().latency_s);
    std::printf("  tremor at the handle    %6.3f deg rms\n", r.in_deg);
    std::printf("  tremor at the utensil   %6.3f deg rms\n", r.out_deg);
    std::printf("  attenuation             %6.1f dB  (%.1f%% reduction)\n",
                r.att_db, r.reduction);
    std::printf("  voluntary tracking err  %6.2f deg   <- must stay small or the\n",
                r.vol_err_deg);
    std::printf("                                        device fights the user\n");
    std::printf("  WFLC mean freq error    %6.3f Hz\n",  r.f_err_hz);
    std::printf("  peak actuator command   %6.3f of 1.0\n", r.peak_u);
    std::printf("  samples saturated       %6.2f %%\n",  r.sat_pct);
    std::printf("  accel samples rejected  %6.1f %%  (expected: gravity is\n",
                r.gate_pct);
    std::printf("                                        unusable during tremor)\n");

    std::printf("\nPhase-alignment sweep (latency_s), the parameter that decides it\n");
    std::printf("%18s%14s%14s\n", "latency_s [ms]", "att [dB]", "reduction [%]");
    for (int ms = -10; ms <= 6; ms += 2) {
        const Result s = run(ms / 1000.0f);
        std::printf("%18d%14.1f%14.1f\n", ms, s.att_db, s.reduction);
    }
    std::printf("\nNegative dB is attenuation. Positive dB means the device is\n");
    std::printf("adding energy at tremor frequency -- making the shaking worse.\n");
    return 0;
}
