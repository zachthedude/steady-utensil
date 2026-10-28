// SPDX-License-Identifier: MIT
// steady/pid.h -- 2-DOF PID for the joint servos of a tremor-cancelling utensil.
//
// No heap, no exceptions, no RTTI, no <algorithm>. Safe to call from a timer ISR.
#pragma once

#include <cmath>

namespace steady {

constexpr float kPi  = 3.14159265358979f;
constexpr float kTwoPi = 6.28318530717959f;
constexpr float kG   = 9.80665f;   // [m/s^2]

inline constexpr float clampf(float v, float lo, float hi) {
    return v < lo ? lo : (v > hi ? hi : v);
}

inline float wrapPi(float a) {
    while (a >  kPi) a -= kTwoPi;
    while (a < -kPi) a += kTwoPi;
    return a;
}

// ---------------------------------------------------------------------------
// Two-degree-of-freedom PID.
//
//   u = Kp*(b*sp - pv) + Ki*Integral(sp - pv) + Kd*(c*sp_rate - rate) + u_ff
//
// Four choices here are specific to cancelling 4-12 Hz tremor and are not
// interchangeable with a generic PID:
//
//  1. D comes from a *measured* rate (motor encoder velocity or the gyro),
//     never from differentiating the error. A backward difference of a 12-bit
//     encoder at 500 Hz is mostly quantisation noise, and the low-pass you
//     would add to tame it burns phase margin in exactly the band that has to
//     stay phase-accurate.
//
//  2. Ki is deliberately tiny and hard-clamped. Integral action contributes
//     -90 deg of phase and buys nothing above ~1 Hz. Its only job is removing
//     static droop caused by the weight of a loaded spoon. Tuned like a
//     process-control loop, this term will make the utensil *amplify* tremor.
//
//  3. Anti-windup is back-calculation against the output actually applied, so
//     the integrator learns about the slew limiter too, not just the clamp.
//
//  4. Slew limiting is a safety requirement, not a tuning nicety. This
//     mechanism operates roughly 100 mm from a person's eyes and mouth, and a
//     saturated servo command must never become a fast swing.
// ---------------------------------------------------------------------------
class Pid {
public:
    struct Gains {
        float kp = 0.0f;    ///< [cmd/rad]
        float ki = 0.0f;    ///< [cmd/(rad*s)]
        float kd = 0.0f;    ///< [cmd/(rad/s)]
        float b  = 1.0f;    ///< setpoint weight on P, 0..1 (<1 softens steps)
        float c  = 0.0f;    ///< setpoint weight on D, 0..1 (0 = pure D-on-measurement)
    };

    struct Limits {
        float out_min     = -1.0f;    ///< [cmd]
        float out_max     =  1.0f;    ///< [cmd]
        float slew        =  1.0e9f;  ///< [cmd/s]
        float i_min       = -0.20f;   ///< hard clamp on the integral term itself
        float i_max       =  0.20f;
        float tracking_tc =  0.05f;   ///< [s] back-calculation constant Tt
    };

    struct Sample {
        float setpoint      = 0.0f;     ///< [rad]
        float setpoint_rate = 0.0f;     ///< [rad/s], 0 if unknown
        float measurement   = 0.0f;     ///< [rad]
        float rate          = 0.0f;     ///< [rad/s] measured, bias-corrected
        float feedforward   = 0.0f;     ///< [cmd]  open-loop term (gravity etc.)
        float dt            = 0.002f;   ///< [s]
    };

    Pid() = default;
    Pid(const Gains& g, const Limits& l) : gains_(g), lim_(l) {}

    void setGains (const Gains&  g) { gains_ = g; }
    void setLimits(const Limits& l) { lim_   = l; }

    /// Re-seed so that the next update() returns `output_now` with no transient.
    /// Call on enable, on gain change, and on any sensor-fault recovery.
    void reset(float output_now = 0.0f) {
        integ_     = 0.0f;
        out_       = clampf(output_now, lim_.out_min, lim_.out_max);
        seeded_    = false;
        saturated_ = false;
    }

    float update(const Sample& s) {
        const float dt = (s.dt > 1.0e-6f) ? s.dt : 1.0e-6f;

        const float err    = s.setpoint - s.measurement;              // for I
        const float err_p  = gains_.b * s.setpoint - s.measurement;   // for P
        const float rate_e = gains_.c * s.setpoint_rate - s.rate;     // for D

        const float p = gains_.kp * err_p;
        const float d = gains_.kd * rate_e;

        if (!seeded_) {
            // Bumpless transfer: pick I so this first tick reproduces out_.
            integ_  = clampf(out_ - (p + d + s.feedforward), lim_.i_min, lim_.i_max);
            seeded_ = true;
        }

        const float u_raw = p + integ_ + d + s.feedforward;

        // Magnitude limit first, then slew limit, so slew is the binding
        // constraint reported back to the integrator.
        float u = clampf(u_raw, lim_.out_min, lim_.out_max);
        const float max_step = lim_.slew * dt;
        u = clampf(u, out_ - max_step, out_ + max_step);

        // Back-calculation anti-windup:  dI/dt = Ki*e + (u_applied - u_raw)/Tt
        const float tt = (lim_.tracking_tc > 1.0e-4f) ? lim_.tracking_tc : 1.0e-4f;
        integ_ += (gains_.ki * err + (u - u_raw) / tt) * dt;
        integ_  = clampf(integ_, lim_.i_min, lim_.i_max);

        saturated_ = std::fabs(u - u_raw) > 1.0e-6f;
        out_       = u;
        return u;
    }

    float output()    const { return out_;   }
    float integral()  const { return integ_; }
    bool  saturated() const { return saturated_; }

private:
    Gains  gains_{};
    Limits lim_{};
    float  integ_ = 0.0f;
    float  out_   = 0.0f;
    bool   seeded_    = false;
    bool   saturated_ = false;
};

} // namespace steady
