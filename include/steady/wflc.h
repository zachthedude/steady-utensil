// SPDX-License-Identifier: MIT
// steady/wflc.h -- Weighted-Frequency Fourier Linear Combiner.
//
// Adaptive estimator for a quasi-periodic disturbance of unknown, slowly
// drifting frequency (Riviere / Rader / Thakor, IEEE TBME 1998, originally
// developed for microsurgical tremor cancellation).
//
// WHY THIS EXISTS ALONGSIDE THE PID
// ---------------------------------
// A feedback loop cannot cancel a 4-12 Hz disturbance without enormous
// bandwidth, because everything in the chain -- IMU filtering, the control
// period, the actuator rise time -- contributes phase lag. Around 90 deg of
// accumulated lag turns cancellation into amplification and the utensil makes
// the tremor WORSE than no device at all.
//
// The way out is to stop reacting and start predicting. Tremor is narrowband
// and close to sinusoidal, so a two-harmonic model locks onto it in well under
// a second. Once locked, predict() evaluates that model at (t + t_ahead),
// cancelling pipeline lag outright instead of fighting it. The PID is then
// left with the easy job: tracking a joint setpoint that is already correct.
//
// Measured on the reference plant (see validate/): the optimal t_ahead tracks
// added sensing delay one-for-one, holding -17 dB of attenuation across at
// least 10 ms of pipeline delay. Without prediction the same delay costs
// 9 dB and eventually inverts the sign of the correction.
//
// TWO NON-OBVIOUS DETAILS, both of which cost real performance if skipped:
//
//  1. The frequency gradient uses the FUNDAMENTAL ONLY (freq_harmonics = 1).
//     The gradient of harmonic r is weighted by r, so a mere 25% second
//     harmonic contributes half the frequency drive of the fundamental it is
//     supposed to be a detail of. Including it pulls omega clean off lock:
//     measured frequency error goes from 0.018 Hz to 3.7 Hz, i.e. no lock at
//     all. Harmonics are still reconstructed by the weight vector; they just
//     do not get a vote on what the frequency is.
//
//  2. Adaptation is in RATE form (multiplied by dt), so mu_w and mu_omega do
//     not have to be retuned when the control rate changes.
// ---------------------------------------------------------------------------
#pragma once

#include <cmath>
#include "pid.h"

namespace steady {

// M = number of harmonics reconstructed. 2 covers the mildly non-sinusoidal
// shape of essential tremor; M = 3 measured no better and costs sin/cos calls.
template <int M = 2>
class Wflc {
public:
    struct Config {
        float mu_w     = 8.0f;    ///< [1/s]     weight adaptation rate
        float mu_omega = 200.0f;  ///< [rad/s^2] frequency adaptation rate
        float f_init   = 6.0f;    ///< [Hz] initial guess, mid essential-tremor
        float f_min    = 3.0f;    ///< [Hz] lock floor: below is voluntary motion
        float f_max    = 14.0f;   ///< [Hz] lock ceiling: above is sensor noise
        /// Harmonics allowed to drive the frequency estimate. Keep this at 1.
        /// See note 1 in the file header before changing it.
        int   freq_harmonics = 1;
    };

    void configure(const Config& c) {
        cfg_ = c;
        mf_  = (c.freq_harmonics < 1) ? 1 : (c.freq_harmonics > M ? M
                                                                 : c.freq_harmonics);
        reset();
    }

    void reset() {
        for (int i = 0; i < 2 * M; ++i) w_[i] = 0.0f;
        phase_ = 0.0f;
        omega_ = kTwoPi * cfg_.f_init;
    }

    /// Feed one band-conditioned tremor sample. Returns the model estimate.
    float update(float s, float dt) {
        float x[2 * M];
        for (int r = 1; r <= M; ++r) {
            const float rp = static_cast<float>(r) * phase_;
            x[r - 1]     = std::sin(rp);
            x[M + r - 1] = std::cos(rp);
        }

        float y = 0.0f;
        for (int i = 0; i < 2 * M; ++i) y += w_[i] * x[i];
        const float e = s - y;

        // Frequency gradient d(y)/d(phase), fundamental only.
        float grad = 0.0f;
        float amp2 = 0.0f;
        for (int r = 1; r <= mf_; ++r) {
            const float a = w_[r - 1], b = w_[M + r - 1];
            grad += static_cast<float>(r) * (a * x[M + r - 1] - b * x[r - 1]);
            amp2 += a * a + b * b;
        }

        // Normalise by the model's own amplitude so lock speed does not depend
        // on how hard the user happens to be shaking. Without this, mu_omega
        // needs retuning per user and per meal.
        omega_ += 2.0f * cfg_.mu_omega * e * grad / (amp2 + kEps) * dt;
        omega_  = clampf(omega_, kTwoPi * cfg_.f_min, kTwoPi * cfg_.f_max);

        // LMS weight update. ||x||^2 == M identically (sin^2 + cos^2 per
        // harmonic), so the regressor is already normalised and plain LMS is
        // well conditioned.
        const float g = 2.0f * cfg_.mu_w * e * dt;
        for (int i = 0; i < 2 * M; ++i) w_[i] += g * x[i];

        phase_ = wrapPi(phase_ + omega_ * dt);
        return y;
    }

    /// Evaluate the locked model t_ahead seconds into the future.
    ///
    /// t_ahead is a NET PHASE ALIGNMENT, not a physical delay: it is the
    /// pipeline lag minus whatever phase lead the conditioning filters
    /// contribute. It is legitimately negative when the front end leads --
    /// on the reference build the optimum is -3 ms. Find it by sweep; the
    /// curve has a clear single minimum and is worth about 5 dB.
    float predict(float t_ahead) const { return eval(t_ahead, 0); }

    /// Analytic d/dt of predict(). Handing this to the PID as the D-term
    /// setpoint rate avoids differentiating a noisy setpoint.
    float predictRate(float t_ahead) const { return eval(t_ahead, 1); }

    /// Analytic d2/dt2 of predict(). Used for the world-frame acceleration
    /// feedforward on geared mechanisms; see stabilizer.h.
    float predictAccel(float t_ahead) const { return eval(t_ahead, 2); }

    float frequencyHz() const { return omega_ / kTwoPi; }

    /// Peak amplitude of the fundamental [rad]. Useful for telemetry, and for
    /// deciding whether stabilisation is worth spending battery on at all.
    float amplitude() const { return std::sqrt(w_[0] * w_[0] + w_[M] * w_[M]); }

private:
    float eval(float t_ahead, int order) const {
        const float ph = phase_ + omega_ * t_ahead;
        float out = 0.0f;
        for (int r = 1; r <= M; ++r) {
            const float rf = static_cast<float>(r);
            const float rp = rf * ph;
            const float s_ = std::sin(rp), c_ = std::cos(rp);
            const float a = w_[r - 1], b = w_[M + r - 1];
            if (order == 0)      out += a * s_ + b * c_;
            else if (order == 1) out += rf * omega_ * (a * c_ - b * s_);
            else                 out -= rf * rf * omega_ * omega_ * (a * s_ + b * c_);
        }
        return out;
    }

    static constexpr float kEps = 1.0e-8f;
    Config cfg_{};
    int    mf_ = 1;
    float  w_[2 * M] = {0.0f};
    float  omega_ = kTwoPi * 6.0f;
    float  phase_ = 0.0f;
};

} // namespace steady
