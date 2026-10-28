# steady-utensil

**Project page:** https://zachthedude.github.io/steady-utensil/ — interactive look at why prediction beats reaction.

Tremor-cancelling control firmware for a 3D-printed self-stabilising eating
utensil. Two-axis gimbal, IMU on the handle, BLDC gimbal motors, PID inner
loops with an adaptive predictive feedforward.

---

## The design correction you should read first

The obvious specification — *"the IMU detects the shaking and the servos move
the spoon in the exact opposite direction, keeping the food perfectly level"* —
does not survive contact with the problem, for two reasons.

**1. "Level" is the wrong target.** Holding the utensil level in the world frame
means resisting the user every time they raise it toward their mouth. Measured
below: a hold-level controller produces **15.0° of RMS voluntary-tracking
error**. The device would feel like it was fighting you, because it would be.

Tremor is 4–12 Hz. Voluntary reaching is under about 2 Hz. The controller has to
*separate* those bands and cancel only the tremor, while tracking intent as the
setpoint. The kinematics then collapse to something simple: for a serial gimbal
where `θ_utensil = θ_handle + q`,

```
q_desired = θ_intent − θ_handle = −tremor(t)
```

The joint setpoint is exactly the negated tremor. The whole problem is
estimating `tremor(t)` accurately at the moment the actuator will get there.

**2. "Instantly" is not available.** IMU filtering, the control period and the
actuator rise time all contribute phase lag. At roughly 90° of accumulated lag,
cancellation inverts into amplification and the device makes the tremor *worse*.
A purely reactive loop — detect the shake, oppose the shake — measures
**−5.6 dB** here against **−17.5 dB** for the predictive design. Add 10 ms of
sensing delay and the uncompensated loop collapses to −5.2 dB, while the
predictive one holds −16.9 dB.

The fix is to stop reacting and start predicting. Tremor is narrowband, so a
two-harmonic adaptive oscillator (WFLC) locks onto it in under a second and can
be evaluated *ahead* of real time, cancelling pipeline lag instead of fighting
it. The PID is then left with a much easier job: tracking a setpoint that is
already correct.

---

## Measured results

30 s run. Tremor chirps 5.0 → 9.0 Hz with a 25% second harmonic and ±30%
amplitude modulation, superimposed on a 0.25 Hz reaching motion. Plant includes
transport delay, finite actuator authority, encoder quantisation, velocity-
estimate noise, gyro bias and noise, Coulomb friction, motor cogging, base-
acceleration coupling, and an accelerometer corrupted by the tremor's own
linear acceleration. Attenuation measured in the 3–14 Hz band over t = 8–30 s.

| Controller | Tremor out | Attenuation | Voluntary error |
|---|---|---|---|
| **Full design** (WFLC predict + PID) | 0.35° rms | **−17.5 dB (86.7% reduction)** | **0.20°** |
| Reactive PID only (no prediction) | 1.37° rms | −5.6 dB (47.8%) | 0.46° |
| Naive hold-level in world frame | 1.28° rms | −6.2 dB (51.3%) | **15.0°** ← fights the user |

Input tremor at the handle: 2.63° rms.

Supporting numbers for the full design: frequency-tracking error **0.108 Hz**,
peak actuator command **0.54 of 1.0**, **0%** of samples saturated, inner-loop
step response **0.011° static error with no overshoot**.

---

## Architecture

```
IMU 1 kHz ─► tilt observer (complementary PI, 0.3 Hz crossover, gated)
                │
                ├─► intent   = LP 1.2 Hz          ── what the user means to do
                └─► tremor   = (θ − intent), conditioned 0.8–40 Hz
                        │
                        ▼
                     WFLC ── locks frequency, amplitude, phase
                        │
                        ├─ predict(latency)      ─► joint setpoint  q_des
                        ├─ predictRate(latency)  ─► D-term setpoint rate
                        └─ predictAccel(latency) ─► accel feedforward (geared only)
                        │
                        ▼
        travel limit ─► slew limit ─► PID (D from measured rate) ─► motor
```

Files:

| Path | Contents |
|---|---|
| `include/steady/pid.h` | The PID. 2-DOF, back-calculation anti-windup, slew limiting, bumpless transfer |
| `include/steady/filters.h` | Biquads (RBJ, DF2T) and the gated tilt observer |
| `include/steady/wflc.h` | Weighted-Frequency Fourier Linear Combiner + predictors |
| `include/steady/stabilizer.h` | Kinematics, per-axis controller, safety state machine |
| `include/steady/tuning.h` | The reference configuration, with each value's justification |
| `firmware/main.cpp` | Real-time task structure, ISR/control split, HAL surface |
| `sim/sim_main.cpp` | Desktop closed-loop simulation against the shipping headers |
| `validate/validate.py` | Python transliteration + the probes that produced the numbers above |

---

## Five findings that cost real performance

Each of these was found by sweep, and each looks like a harmless choice until
you measure it.

**1. The conditioning low-pass is set by the tremor's *harmonic*, not the
tremor.** A 20 Hz low-pass looks generous for a 14 Hz signal — but the second
harmonic of a 9 Hz tremor is at 18 Hz, where that filter applies −22% gain and
−81° of phase. The estimator then faithfully reproduces the damage. Opening the
corner to 40 Hz took the open-loop residual from **−13 dB to −21 dB** on an
identical signal. Rule: `band_hi ≥ 3 × f_max`.

**2. Only the fundamental may drive the frequency estimate.** The WFLC frequency
gradient weights harmonic *r* by *r*, so a mere 25% second harmonic contributes
half the frequency drive of the fundamental it is supposed to be a detail of.
Including it destroys the lock — frequency error goes from **0.018 Hz to
3.7 Hz**, i.e. no lock at all. Harmonics still get reconstructed by the weight
vector; they just don't get a vote on what the frequency is.

**3. The accelerometer is useless in the tremor band.** A 6 Hz tremor at the
wrist generates ~7 m/s² of linear acceleration — comparable to gravity. The tilt
observer crossover therefore sits at **0.3 Hz**, not the usual 1–5 Hz, and gates
out accelerometer samples whose magnitude doesn't look like gravity. **31% of
samples get rejected during tremor**, which is the correct behaviour. Set this
crossover conventionally and the loop injects the tremor's own linear
acceleration back into its attitude estimate.

**4. The D-term setpoint weight must be 1, not 0.** Textbook practice is
derivative-on-measurement (`c = 0`) to avoid derivative kick. Here the setpoint
*itself* moves at tremor frequency, and its derivative is known analytically
from the model — so there is no kick to avoid, while `c = 0` would make the D
term fight every commanded motion to the tune of ~2.0 of command against a 1.0
limit. Instant saturation.

**5. A free gimbal is already passively isolating.** A head that stays
inertially still requires no net torque, so `q̈_des + θ̈_handle ≈ 0` in the tremor
band. On a direct-drive gimbal the acceleration feedforward should therefore be
*zero*; feeding `q̈_des` alone commands roughly double the required torque and
saturates the actuator. This is why `accel_ff_gain` defaults to 0 — and why it
must be set to `J/Kt` on a **geared** drive, where backdrive resistance destroys
that passive isolation. It also makes low bearing friction the single most
important mechanical requirement.

---

## Tuning on real hardware

1. **Measure the latency.** Toggle a GPIO on IMU data-ready, clear it after the
   motor driver write, look at it on a scope. Do not estimate it from
   datasheets.
2. **Measure `kLoadFeedforward`**: the command needed to hold a fully loaded
   head out horizontally. This absorbs magnet tolerance, driver scaling and
   gearbox efficiency in one number.
3. **Tune the inner loop first**, with the WFLC disabled, on step responses.
   `kp` sits on a broad plateau — 19 and 45 measured identically; 60 saturated
   96% of samples and 80 went unstable. Pick the low end for margin, quieter
   operation and lower current.
4. **Sweep `latency_s` last.** It is a *net phase alignment* — pipeline lag
   minus front-end phase lead — so it is legitimately negative (−3 ms on the
   reference build). The curve has one clear minimum and is worth about 10 dB:

   | `latency_s` | −6 ms | −4 ms | **−3 ms** | −2 ms | 0 ms | +4 ms |
   |---|---|---|---|---|---|---|
   | attenuation | −13.9 dB | −17.1 dB | **−17.5 dB** | −16.5 dB | −13.1 dB | −7.6 dB |

`latency_s` trades one-for-one against pipeline delay: with 10 ms of extra
sensing delay, moving the predict-ahead by the same 10 ms recovers **−16.9 dB**
versus **−5.2 dB** uncompensated. That tolerance is what lets you use cheaper
sensors and actuators.

Note this compensates *feedforward timing*, not loop stability. Delay added
inside the joint servo loop destabilises it, and no amount of prediction fixes
that — you must detune the PID instead.

---

## Safety

This mechanism operates about 100 mm from a person's eyes and mouth, so the
control law is bounded independently of whether it is working correctly:

- **Setpoint slew limit** (`max_joint_rate`, 12 rad/s ≈ 3.5× what a real tremor
  needs). No fault path — diverged estimator, stuck sensor, corrupt gain table —
  can become a fast swing.
- **Travel limit** enforced before the slew limit, backed by mechanical stops.
- **Latched fault state.** IMU timeout, non-finite attitude or power fault
  drives the joints to neutral and *stays* there until an explicit
  `clearFault()`, so an intermittent sensor can't make the utensil twitch in and
  out of service during a meal.
- **Fail to coasting, not braking.** A braked joint makes the utensil a rigid
  lever that transmits every tremor to the food; a coasting gimbal keeps its
  passive inertial isolation. Either way it degrades to an ordinary utensil,
  which is still perfectly usable. **It must never fail into motion.**
- Watchdog fed only from the control task, never from the ISR.

---

## Mechanical and 3D-printing notes

**Printed parts must never touch food.** FDM layer lines harbour bacteria and
cannot be reliably sanitised, and few printable polymers survive a dishwasher.
The utensil head must be off-the-shelf **stainless steel, removable and
dishwasher-safe**; the printed work is the housing and gimbal structure only.

- **Bearings are the performance-critical part.** Use miniature ball bearings,
  not printed journals. Coulomb friction is what destroys the passive isolation
  described in finding 5 — a printed plain bearing has roughly an order of
  magnitude more of it and will eat most of the attenuation.
- **Material:** PETG or SLS PA12 for the housing. PLA creeps and is brittle;
  ASA if UV exposure matters. TPU overmould for the grip sleeve.
- **Handle diameter 32–38 mm.** Larger diameters need less pinch force, which
  is the point for reduced grip strength or arthritis. Provide a universal-cuff
  variant (D-ring strap) for users who cannot grip at all.
- **Mass budget under ~200 g**, centre of gravity inside the grip. Extra mass
  damps tremor passively but fatigues the user — that trade-off should be
  deliberate.
- **Print orientation:** handle axis vertical, so hoop stress doesn't load the
  layer-adhesion axis. 4 perimeters, ~40% gyroid infill.
- **Head attachment:** magnetic with a mechanical detent that releases under
  overload. If the mechanism ever jams near the face, the head must let go.
- **Sealing:** IP54 minimum with a gasketed USB-C flap. Food and liquid get
  everywhere.
- **Accessibility of the device itself:** one large tactile, high-contrast
  power control; no small buttons; audible *and* haptic state feedback; a
  charging connector that can be mated one-handed without fine alignment.

---

## Verification status

The numbers in this README were produced by `validate/validate.py`, a
line-for-line transliteration of the C++ headers, run against the plant model
described above.

**The C++ has not been compiled** — there is no C++ toolchain on the machine
where this was written. `sim/sim_main.cpp` is provided so the same sweeps can be
re-run against the shipping headers; expect to fix compile errors on first
build.

The synthetic tremor is periodic by construction, which flatters the approach:
`residual_gain` measures as having no effect here precisely because there is no
non-periodic component for it to catch, and real tremor has one. Treat these
figures as evidence the control law is sound, not as predicted clinical
performance. Real validation needs recorded tremor from real users, and
anything intended for actual use needs clinical evaluation and the appropriate
medical-device process for your jurisdiction.
