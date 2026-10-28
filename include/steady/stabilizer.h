// SPDX-License-Identifier: MIT
// steady/stabilizer.h -- per-axis and 2-DOF tremor-cancelling controller.
//
// KINEMATICS
// ----------
// Serial 2-DOF gimbal. The utensil head hangs off two orthogonal joints driven
// relative to the handle, so for each axis:
//
//     theta_utensil_world = theta_handle_world + q          (q = joint angle)
//
// The goal is NOT "hold the utensil level". Holding it level would fight the
// user every time they raise it toward their mouth. The goal is
//
//     theta_utensil_world = theta_intent
//
// where theta_intent is the low-frequency component of the handle attitude --
// what the person is actually trying to do. Substituting:
//
//     q_desired = theta_intent - theta_handle = -tremor(t)
//
// So the joint setpoint is exactly the negated tremor, and the whole problem
// reduces to estimating tremor(t) accurately at the moment the actuator will
// arrive there. That is what the WFLC does; the PID then only has to track q.
// ---------------------------------------------------------------------------
#pragma once

#include "pid.h"
#include "filters.h"
#include "wflc.h"

namespace steady {

class AxisStabilizer {
public:
    struct Config {
        float fs        = 500.0f;   ///< control rate [Hz]
        float intent_hz = 1.2f;     ///< voluntary-motion bandwidth

        /// Conditioning band fed to the estimator. These are NOT "the tremor
        /// band" -- band_hi_hz must clear the tremor HARMONICS or the front
        /// end mangles them and the estimator faithfully reproduces the
        /// damage. Rule: band_hi_hz >= 3 * wflc.f_max. Setting it to 20 Hz
        /// with f_max = 14 Hz measured -13 dB; opening it to 40 Hz measured
        /// -21 dB on the identical signal. It is the largest single win in
        /// the whole front end and it looks like a harmless filter choice.
        float band_lo_hz = 0.8f;
        float band_hi_hz = 40.0f;

        /// Net phase alignment handed to Wflc::predict(). Pipeline lag minus
        /// front-end phase lead, so it is legitimately NEGATIVE here. Sweep
        /// it on real hardware: single clear minimum, worth about 5 dB.
        float latency_s      = -0.003f;
        float travel_rad     = 0.44f;   ///< mechanical stop, +/- 25 deg
        float max_joint_rate = 12.0f;   ///< [rad/s] safety cap on setpoint slew

        /// Cancels the part of the tremor the periodic model does not explain.
        /// Applied with lag, hence sub-unity. Perfectly periodic test signals
        /// make this look useless; real tremor is not perfectly periodic.
        float residual_gain  = 0.25f;

        /// Inverse actuator gain, [cmd per rad/s^2], for the WORLD-FRAME
        /// acceleration feedforward -- see the note in update().
        ///
        /// Leave at 0 for a direct-drive gimbal: a low-friction gimbal is
        /// already passively isolating, because a head that stays inertially
        /// still needs no net torque. Set it to J/Kt for a GEARED drive,
        /// where backdrive resistance destroys that passive isolation.
        float accel_ff_gain  = 0.0f;
        Wflc<2>::Config wflc{};
        Pid::Gains      pid{};
        Pid::Limits     pid_lim{};
    };

    struct Telemetry {
        float intent     = 0.0f;   ///< [rad] estimated voluntary attitude
        float tremor     = 0.0f;   ///< [rad] band-passed disturbance, this sample
        float q_des      = 0.0f;   ///< [rad] joint setpoint after limiting
        float tremor_hz  = 0.0f;   ///< [Hz]  current WFLC lock
        float tremor_amp = 0.0f;   ///< [rad] fundamental amplitude
        float command    = 0.0f;   ///< [cmd] actuator output
        bool  saturated  = false;
    };

    void configure(const Config& c) {
        cfg_ = c;
        intent_lp_.lowpass(c.intent_hz, c.fs, 0.70710678f);
        band_hi_.highpass(c.band_lo_hz, c.fs, 0.70710678f);
        band_lo_.lowpass (c.band_hi_hz, c.fs, 0.70710678f);
        wflc_.configure(c.wflc);
        pid_.setGains(c.pid);
        pid_.setLimits(c.pid_lim);
    }

    /// Seed every filter from the current pose so nothing lurches on enable.
    void reset(float handle_angle, float joint_angle) {
        intent_lp_.reset(handle_angle);
        band_hi_.reset(0.0f);
        band_lo_.reset(0.0f);
        wflc_.reset();
        pid_.reset(0.0f);
        q_des_ = joint_angle;
    }

    /// @param handle_angle  fused handle attitude on this axis [rad]
    /// @param joint_angle   measured joint angle from the encoder [rad]
    /// @param joint_rate    measured joint velocity [rad/s]
    /// @param gravity_ff    open-loop command offsetting the utensil load [cmd]
    /// @param base_accel    handle angular acceleration [rad/s^2], from a
    ///                      filtered gyro derivative. Only used when
    ///                      accel_ff_gain is non-zero; pass 0 otherwise.
    /// @param enable        false: ramp the joint to neutral and hold there
    /// @return              actuator command, clamped to the PID output limits
    float update(float handle_angle,
                 float joint_angle,
                 float joint_rate,
                 float gravity_ff,
                 float base_accel,
                 float dt,
                 bool  enable,
                 Telemetry* tm = nullptr) {

        // 1. Split handle motion into intent and tremor.
        const float intent = intent_lp_.process(handle_angle);
        const float hp     = handle_angle - intent;            // complementary HP
        const float s      = band_lo_.process(band_hi_.process(hp));

        // 2. Lock the adaptive model onto the tremor.
        const float model = wflc_.update(s, dt);
        const float resid = s - model;

        // 3. Joint setpoint = negated tremor, evaluated at (t + latency) so the
        //    actuator ARRIVES in phase instead of chasing. The residual term
        //    mops up the non-periodic remainder; that part is applied with lag,
        //    hence the deliberately sub-unity gain.
        float q_target, q_target_rate, q_target_accel;
        if (enable) {
            q_target       = -(wflc_.predict(cfg_.latency_s) +
                               cfg_.residual_gain * resid);
            q_target_rate  = -wflc_.predictRate (cfg_.latency_s);
            q_target_accel = -wflc_.predictAccel(cfg_.latency_s);
        } else {
            q_target       = 0.0f;   // fail-safe pose: a plain, rigid utensil
            q_target_rate  = 0.0f;
            q_target_accel = 0.0f;
        }

        // 4. Limit travel, then slew. Setpoint slew is bounded to
        //    max_joint_rate so that no fault path -- a diverged estimator, a
        //    stuck sensor, a bad gain table -- can become a fast swing while
        //    the head is near the eyes or mouth.
        q_target = clampf(q_target, -cfg_.travel_rad, cfg_.travel_rad);
        const float max_step = cfg_.max_joint_rate * dt;
        const float q_prev   = q_des_;
        q_des_ = clampf(q_target, q_prev - max_step, q_prev + max_step);
        if (q_des_ != q_target) {          // limiter is binding: stop lying to
            q_target_rate  = (q_des_ - q_prev) / dt;   // the PID about where
            q_target_accel = 0.0f;                     // the setpoint is going
        }

        // 5. Inner loop: make the joint actually get there under load.
        Pid::Sample smp;
        smp.setpoint      = q_des_;
        smp.setpoint_rate = q_target_rate;
        smp.measurement   = joint_angle;
        smp.rate          = joint_rate;
        // World-frame acceleration feedforward. What the motor must supply is
        // set by the head's acceleration in the WORLD, not in the handle
        // frame, and (q_ddot + base_accel) is exactly that. In the tremor band
        // those two very nearly cancel -- which is the algebraic statement of
        // the fact that a free gimbal already isolates rotation. Feeding
        // q_ddot alone (without base_accel) commands roughly double the
        // required torque and pushes the actuator into saturation.
        smp.feedforward   = enable
                          ? gravity_ff +
                            cfg_.accel_ff_gain * (q_target_accel + base_accel)
                          : 0.0f;
        smp.dt            = dt;
        const float u = pid_.update(smp);

        if (tm) {
            tm->intent     = intent;
            tm->tremor     = s;
            tm->q_des      = q_des_;
            tm->tremor_hz  = wflc_.frequencyHz();
            tm->tremor_amp = wflc_.amplitude();
            tm->command    = u;
            tm->saturated  = pid_.saturated();
        }
        return u;
    }

    Pid&     pid()  { return pid_;  }
    Wflc<2>& wflc() { return wflc_; }

private:
    Config  cfg_{};
    Biquad  intent_lp_{}, band_hi_{}, band_lo_{};
    Wflc<2> wflc_{};
    Pid     pid_{};
    float   q_des_ = 0.0f;
};

// ---------------------------------------------------------------------------
// Two-axis front end plus the safety state machine.
//
// FAIL-SAFE PHILOSOPHY: every fault path drives the joints to neutral and
// holds them there. A stabilised utensil that fails becomes an ordinary
// utensil, which is still perfectly usable. It must never fail into motion.
// ---------------------------------------------------------------------------
class Stabilizer {
public:
    enum class State : unsigned char { kInit, kCalibrating, kActive, kFault };

    struct Inputs {
        float pitch_angle = 0.0f, roll_angle = 0.0f;  ///< fused handle attitude [rad]
        float q_pitch  = 0.0f, q_roll  = 0.0f;        ///< joint encoders [rad]
        float qd_pitch = 0.0f, qd_roll = 0.0f;        ///< joint velocities [rad/s]
        float base_accel_pitch = 0.0f;                ///< [rad/s^2] filtered
        float base_accel_roll  = 0.0f;                ///< gyro derivative
        float dt = 0.002f;
        bool  imu_healthy = true;
        bool  power_ok    = true;
    };

    struct Outputs {
        float cmd_pitch = 0.0f, cmd_roll = 0.0f;
        State state = State::kInit;
        AxisStabilizer::Telemetry pitch{}, roll{};
    };

    void configure(const AxisStabilizer::Config& pitch_cfg,
                   const AxisStabilizer::Config& roll_cfg,
                   float calibrate_s = 1.0f) {
        pitch_.configure(pitch_cfg);
        roll_.configure(roll_cfg);
        calibrate_s_ = calibrate_s;
        state_ = State::kInit;
    }

    /// Load compensation. Geometry: theta is the utensil pitch measured from
    /// horizontal, so the gravity moment arm goes as cos(theta). kff has units
    /// of [cmd] and is the command needed to hold a fully loaded head out
    /// horizontally; measure it once on the bench rather than deriving it.
    void setLoadFeedforward(float kff) { kff_ = kff; }

    Outputs update(const Inputs& in) {
        Outputs out;

        const bool faulted = !in.imu_healthy || !in.power_ok ||
                             !std::isfinite(in.pitch_angle) ||
                             !std::isfinite(in.roll_angle);

        switch (state_) {
        case State::kInit:
            pitch_.reset(in.pitch_angle, in.q_pitch);
            roll_.reset (in.roll_angle,  in.q_roll);
            timer_ = 0.0f;
            state_ = State::kCalibrating;
            break;

        case State::kCalibrating:
            // Hold neutral while the tilt observer settles and the gyro bias
            // estimate converges. Roughly one second with the device at rest.
            timer_ += in.dt;
            if (timer_ >= calibrate_s_) state_ = State::kActive;
            if (faulted) state_ = State::kFault;
            break;

        case State::kActive:
            if (faulted) state_ = State::kFault;
            break;

        case State::kFault:
            // Latched. Recovery requires an explicit clearFault(), so an
            // intermittent sensor cannot make the utensil twitch in and out of
            // service during a meal.
            break;
        }

        const bool  enable   = (state_ == State::kActive);
        const float ff_pitch = kff_ * std::cos(in.pitch_angle + in.q_pitch);

        out.cmd_pitch = pitch_.update(in.pitch_angle, in.q_pitch, in.qd_pitch,
                                      ff_pitch, in.base_accel_pitch,
                                      in.dt, enable, &out.pitch);
        out.cmd_roll  = roll_.update (in.roll_angle,  in.q_roll,  in.qd_roll,
                                      0.0f,     in.base_accel_roll,
                                      in.dt, enable, &out.roll);
        out.state = state_;
        return out;
    }

    void  clearFault()  { if (state_ == State::kFault) state_ = State::kInit; }
    State state() const { return state_; }

    AxisStabilizer& pitchAxis() { return pitch_; }
    AxisStabilizer& rollAxis()  { return roll_;  }

private:
    AxisStabilizer pitch_{}, roll_{};
    State state_ = State::kInit;
    float timer_ = 0.0f;
    float calibrate_s_ = 1.0f;
    float kff_ = 0.0f;
};

} // namespace steady
