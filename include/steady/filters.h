// SPDX-License-Identifier: MIT
// steady/filters.h -- biquads and the tilt observer.
#pragma once

#include <cmath>
#include "pid.h"   // clampf, wrapPi, kPi, kG

namespace steady {

// ---------------------------------------------------------------------------
// Direct-Form-II-transposed biquad, RBJ cookbook coefficients.
// DF2T is chosen over DF1 for its better round-off behaviour at the low
// normalised frequencies used here (1.2 Hz at fs = 500 Hz gives w0 = 0.015 rad).
// ---------------------------------------------------------------------------
class Biquad {
public:
    void lowpass(float fc, float fs, float q = 0.70710678f) {
        float cw, alpha; prewarp(fc, fs, q, cw, alpha);
        const float a0 = 1.0f + alpha;
        b0_ = ((1.0f - cw) * 0.5f) / a0;
        b1_ = (1.0f - cw) / a0;
        b2_ = b0_;
        a1_ = (-2.0f * cw) / a0;
        a2_ = (1.0f - alpha) / a0;
    }

    void highpass(float fc, float fs, float q = 0.70710678f) {
        float cw, alpha; prewarp(fc, fs, q, cw, alpha);
        const float a0 = 1.0f + alpha;
        b0_ = ((1.0f + cw) * 0.5f) / a0;
        b1_ = -(1.0f + cw) / a0;
        b2_ = b0_;
        a1_ = (-2.0f * cw) / a0;
        a2_ = (1.0f - alpha) / a0;
    }

    /// Constant 0 dB peak-gain bandpass.
    void bandpass(float fc, float fs, float q = 1.0f) {
        float cw, alpha; prewarp(fc, fs, q, cw, alpha);
        const float a0 = 1.0f + alpha;
        b0_ =  alpha / a0;
        b1_ =  0.0f;
        b2_ = -alpha / a0;
        a1_ = (-2.0f * cw) / a0;
        a2_ = (1.0f - alpha) / a0;
    }

    inline float process(float x) {
        const float y = b0_ * x + z1_;
        z1_ = b1_ * x - a1_ * y + z2_;
        z2_ = b2_ * x - a2_ * y;
        return y;
    }

    /// Seed the state so a constant input x produces its steady-state output
    /// immediately. Without this the intent low-pass takes seconds to charge on
    /// power-up and the utensil lurches while it settles.
    void reset(float x = 0.0f) {
        const float den = 1.0f + a1_ + a2_;
        const float g   = (std::fabs(den) > 1.0e-12f) ? (b0_ + b1_ + b2_) / den : 0.0f;
        const float y   = g * x;
        z1_ = y - b0_ * x;
        z2_ = b2_ * x - a2_ * y;
    }

private:
    static void prewarp(float fc, float fs, float q, float& cw, float& alpha) {
        const float w0 = kTwoPi * fc / fs;
        cw    = std::cos(w0);
        alpha = std::sin(w0) / (2.0f * q);
    }
    float b0_ = 1.0f, b1_ = 0.0f, b2_ = 0.0f, a1_ = 0.0f, a2_ = 0.0f;
    float z1_ = 0.0f, z2_ = 0.0f;
};

// ---------------------------------------------------------------------------
// Single-axis tilt observer (Mahony-style complementary PI).
//
// The crossover is deliberately LOW, around 0.3 Hz. During a 6 Hz tremor of
// only 5 mm amplitude the hand experiences (2*pi*6)^2 * 0.005 = 7.1 m/s^2 of
// linear acceleration, comparable to gravity itself. The accelerometer is
// therefore NOT a usable attitude reference anywhere near the tremor band; its
// only role is bounding gyro drift over seconds. Setting this crossover to the
// usual 1-5 Hz injects the tremor own linear acceleration back into the
// attitude estimate and the loop chases its tail.
// ---------------------------------------------------------------------------
class TiltEstimator {
public:
    void configure(float crossover_hz, float zeta = 1.0f) {
        const float wc = kTwoPi * crossover_hz;
        k_angle_ = 2.0f * zeta * wc;
        k_bias_  = wc * wc;
    }

    void reset(float angle) { angle_ = angle; bias_ = 0.0f; }

    /// @param accel_angle  tilt implied by the accelerometer [rad]
    /// @param accel_mag    magnitude of specific force [m/s^2], gates the fix
    /// @param gyro_rate    raw gyro [rad/s]
    /// @return             bias-corrected angular rate [rad/s]
    float update(float accel_angle, float accel_mag, float gyro_rate, float dt) {
        const float rate = gyro_rate - bias_;
        angle_ += rate * dt;

        // Only trust gravity when the specific force actually looks like
        // gravity alone. Under a strong tremor this gate rejects most samples,
        // which is correct: coast on the gyro instead.
        if (std::fabs(accel_mag - kG) < kAccelGate) {
            const float e = wrapPi(accel_angle - angle_);
            angle_ += k_angle_ * e * dt;
            bias_  -= k_bias_  * e * dt;
            bias_   = clampf(bias_, -kMaxBias, kMaxBias);
            gated_  = false;
        } else {
            gated_ = true;
        }
        angle_ = wrapPi(angle_);
        return rate;
    }

    float angle()      const { return angle_; }
    float bias()       const { return bias_;  }
    bool  accelGated() const { return gated_; }

private:
    static constexpr float kAccelGate = 1.5f;   // [m/s^2]
    static constexpr float kMaxBias   = 0.15f;  // [rad/s]
    float angle_ = 0.0f, bias_ = 0.0f;
    float k_angle_ = 1.0f, k_bias_ = 1.0f;
    bool  gated_ = false;
};

} // namespace steady
