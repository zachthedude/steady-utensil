// SPDX-License-Identifier: MIT
// steady/tuning.h -- the reference build's measured configuration.
//
// Every number here was found by sweep against the plant model in
// validate/validate.py, not chosen by feel. Where a value has a surprising
// justification it is stated inline. Re-run the sweeps after any mechanical
// change; several of these are properties of the hardware, not of the maths.
//
// Reference hardware assumed by these numbers:
//   direct-drive BLDC gimbal motor, ~0.06 N m per unit command
//   joint inertia (loaded head + gimbal ring)  6.0e-5 kg m^2
//   14-bit on-axis magnetic encoder, FOC driver supplying velocity
//   IMU sampled at 1 kHz, control loop at 500 Hz
//   head mass 35 g at 55 mm from the joint (a full spoon)
// ---------------------------------------------------------------------------
#pragma once

#include "stabilizer.h"

namespace steady {
namespace reference {

// --- mechanism constants -------------------------------------------------
constexpr float kJointInertia = 6.0e-5f;   ///< [kg m^2]
constexpr float kMotorKt      = 0.060f;    ///< [N m per unit command]
constexpr float kHeadMass     = 0.035f;    ///< [kg] loaded
constexpr float kHeadArm      = 0.055f;    ///< [m]
constexpr float kControlHz    = 500.0f;

/// Command needed to hold a fully loaded head out horizontally. Measure this
/// on the bench rather than trusting the arithmetic -- it absorbs gearbox
/// efficiency, magnet strength tolerance and driver scaling all at once.
constexpr float kLoadFeedforward = kHeadMass * kG * kHeadArm / kMotorKt;

inline AxisStabilizer::Config axisConfig() {
    AxisStabilizer::Config c;
    c.fs        = kControlHz;
    c.intent_hz = 1.2f;

    // Conditioning band. The high edge is set by the tremor's SECOND
    // HARMONIC, not by the tremor: 2 x 14 Hz = 28 Hz must pass undistorted.
    // Measured -13 dB at 20 Hz vs -21 dB at 40 Hz on the identical signal.
    c.band_lo_hz = 0.8f;
    c.band_hi_hz = 40.0f;

    // Net phase alignment: pipeline lag minus front-end lead. Negative here
    // because the 0.8 Hz high-pass leads at tremor frequencies. Swept:
    // -3 ms gives -17.7 dB, 0 ms gives -13.1 dB, +4 ms gives -7.6 dB.
    c.latency_s = -0.003f;

    c.travel_rad     = 0.44f;   // +/- 25 deg, set by the gimbal geometry
    c.max_joint_rate = 12.0f;   // safety: 3.5x the peak rate a real tremor needs
    c.residual_gain  = 0.25f;
    c.accel_ff_gain  = 0.0f;    // direct drive: the head self-isolates

    c.wflc.mu_w           = 8.0f;
    c.wflc.mu_omega       = 200.0f;
    c.wflc.f_init         = 6.0f;
    c.wflc.f_min          = 3.0f;
    c.wflc.f_max          = 14.0f;
    c.wflc.freq_harmonics = 1;   // do not raise; see wflc.h header note 1

    // Joint servo. kp is on a measured plateau: 19 and 45 perform identically
    // (-13 dB), 60 saturates 96% of samples and 80 goes unstable. 19 is chosen
    // for margin, quieter operation and lower current draw -- not because it
    // is the best-performing value.
    c.pid.kp = 19.0f;
    // Small on purpose: this term exists only to clear the static droop of a
    // loaded head against bearing stiction, not to shape the response. Swept:
    // 30 leaves 0.49 deg of standing error, 120 leaves 0.011 deg with zero
    // overshoot, 500 starts overshooting and costs 0.7 dB of attenuation.
    c.pid.ki = 120.0f;
    c.pid.kd = 0.22f;
    c.pid.b  = 1.0f;
    c.pid.c  = 1.0f;     // D acts on tracking-rate error, NOT on measurement.
                         // Unusual, and correct here: the setpoint itself moves
                         // at tremor frequency and its derivative is known
                         // analytically from the model, so there is no
                         // derivative kick to avoid -- while c = 0 would make
                         // the D term fight every commanded motion, costing
                         // ~2.0 of command against a 1.0 limit.

    c.pid_lim.out_min     = -1.0f;
    c.pid_lim.out_max     =  1.0f;
    c.pid_lim.slew        =  200.0f;   // safety cap, not normally binding
    c.pid_lim.i_min       = -0.35f;
    c.pid_lim.i_max       =  0.35f;
    c.pid_lim.tracking_tc =  0.02f;
    return c;
}

/// Tilt observer crossover. Deliberately far below the usual 1-5 Hz: during a
/// 6 Hz tremor the wrist generates ~7 m/s^2 of linear acceleration, so the
/// accelerometer is not a gravity reference anywhere near the tremor band.
/// Measured: 31% of accelerometer samples get gate-rejected during tremor,
/// which is the correct behaviour.
constexpr float kTiltCrossoverHz = 0.30f;
constexpr float kTiltZeta        = 1.0f;

/// Corner for the filtered gyro derivative used as base_accel. Only needed
/// when accel_ff_gain is non-zero.
constexpr float kBaseAccelLpHz = 30.0f;

} // namespace reference
} // namespace steady
