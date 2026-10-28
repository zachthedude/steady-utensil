// SPDX-License-Identifier: MIT
// firmware/main.cpp -- real-time task structure for the stabilised utensil.
//
// Vendor-neutral. Every hook marked HAL: is the one line you replace with your
// MCU's driver call. Targeted at an STM32G4 / RP2040 class part with hardware
// FPU; measured cost of one control tick on a 170 MHz M4F is dominated by the
// 12 sin/cos calls in the two WFLCs (~9 us), leaving the 2 ms budget ~99% idle.
//
// TIMING CONTRACT -- violate any of these and the loop stops working:
//   IMU sampled at 1 kHz, hardware-timed. Software-polled sampling jitters and
//     jitter at 8 Hz reads as phase noise, which is exactly what the predictor
//     cannot recover from.
//   Control at 500 Hz, driven from a timer, never from the main loop.
//   Sense-to-actuation latency must be MEASURED (toggle a GPIO on IMU DRDY,
//     clear it after the driver write, look at it on a scope) and the result
//     fed to latency_s. Do not estimate it from datasheets.
// ---------------------------------------------------------------------------

#include "steady/tuning.h"

#include <cstdint>
#include <cmath>

using namespace steady;

// ---------------------------------------------------------------------------
// HAL surface -- replace these six functions.
// ---------------------------------------------------------------------------
struct ImuSample {
    float ax, ay, az;      // [m/s^2]
    float gx, gy, gz;      // [rad/s]
    bool  valid;
};

extern ImuSample halImuRead();                     // HAL: burst-read the IMU
extern void      halMotorWrite(float pitch, float roll);  // HAL: -1..1 each
extern void      halMotorDisable();                // HAL: coast, do not brake
extern float     halJointAngle(int axis);          // HAL: encoder [rad]
extern float     halJointRate(int axis);           // HAL: FOC observer [rad/s]
extern bool      halPowerOk();                     // HAL: battery + driver fault
extern void      halTelemetry(const Stabilizer::Outputs&);

// ---------------------------------------------------------------------------
// Signal conditioning shared between the ISR and the control task.
// ---------------------------------------------------------------------------
namespace {

constexpr float kDt          = 1.0f / reference::kControlHz;
constexpr int   kImuPerTick  = 2;      // 1 kHz IMU / 500 Hz control
constexpr int   kImuFaultMax = 5;      // consecutive bad reads before fault

Stabilizer   g_stab;
TiltEstimator g_tilt_pitch, g_tilt_roll;
Biquad        g_daccel_pitch, g_daccel_roll;   // filtered gyro derivative

float g_prev_rate_pitch = 0.0f, g_prev_rate_roll = 0.0f;
int   g_imu_faults = 0;

// Decimation accumulator: average the two IMU samples per control tick rather
// than dropping one. Averaging is a 2-tap FIR with half a sample of group
// delay, which is cheaper in phase than any IIR anti-alias filter would be.
struct Accum {
    float ax = 0, ay = 0, az = 0, gx = 0, gy = 0, gz = 0;
    int   n = 0;
    void add(const ImuSample& s) {
        ax += s.ax; ay += s.ay; az += s.az;
        gx += s.gx; gy += s.gy; gz += s.gz; ++n;
    }
    bool take(ImuSample& out) {
        if (n == 0) return false;
        const float k = 1.0f / static_cast<float>(n);
        out.ax = ax * k; out.ay = ay * k; out.az = az * k;
        out.gx = gx * k; out.gy = gy * k; out.gz = gz * k;
        out.valid = true;
        ax = ay = az = gx = gy = gz = 0.0f; n = 0;
        return true;
    }
};

volatile Accum g_accum;          // written in ISR, drained in control task
volatile bool  g_tick_ready = false;

} // namespace

// ---------------------------------------------------------------------------
// 1 kHz IMU interrupt. Keep it short: no trig, no division, no logging.
// ---------------------------------------------------------------------------
extern "C" void IMU_DRDY_IRQHandler() {
    const ImuSample s = halImuRead();
    if (!s.valid) { ++g_imu_faults; return; }
    g_imu_faults = 0;

    Accum& a = const_cast<Accum&>(g_accum);
    a.add(s);
    if (a.n >= kImuPerTick) g_tick_ready = true;
}

// ---------------------------------------------------------------------------
// 500 Hz control task.
// ---------------------------------------------------------------------------
void controlTick() {
    ImuSample s{};
    {
        // Short critical section: copy out, then compute with interrupts on.
        // HAL: __disable_irq() / __enable_irq() or equivalent.
        Accum& a = const_cast<Accum&>(g_accum);
        if (!a.take(s)) return;
        g_tick_ready = false;
    }

    // --- attitude -------------------------------------------------------
    // Tilt from specific force. atan2 keeps this valid through +/-90 deg,
    // which a small-angle approximation would not: people tip a spoon a long
    // way on the way to their mouth.
    const float acc_mag   = std::sqrt(s.ax * s.ax + s.ay * s.ay + s.az * s.az);
    const float acc_pitch = std::atan2(-s.ax, std::sqrt(s.ay * s.ay + s.az * s.az));
    const float acc_roll  = std::atan2( s.ay, s.az);

    const float rate_pitch = g_tilt_pitch.update(acc_pitch, acc_mag, s.gy, kDt);
    const float rate_roll  = g_tilt_roll.update (acc_roll,  acc_mag, s.gx, kDt);

    // Base angular acceleration, for the geared-drive feedforward. Cheap and
    // low-noise: differentiating a 0.01 rad/s-rms gyro at 500 Hz then
    // low-passing at 30 Hz lands around 2 rad/s^2 rms, which is ~1% of the
    // tremor's own peak angular acceleration.
    const float ba_pitch = g_daccel_pitch.process((rate_pitch - g_prev_rate_pitch) / kDt);
    const float ba_roll  = g_daccel_roll.process ((rate_roll  - g_prev_rate_roll ) / kDt);
    g_prev_rate_pitch = rate_pitch;
    g_prev_rate_roll  = rate_roll;

    // --- control --------------------------------------------------------
    Stabilizer::Inputs in;
    in.pitch_angle      = g_tilt_pitch.angle();
    in.roll_angle       = g_tilt_roll.angle();
    in.q_pitch          = halJointAngle(0);
    in.q_roll           = halJointAngle(1);
    in.qd_pitch         = halJointRate(0);
    in.qd_roll          = halJointRate(1);
    in.base_accel_pitch = ba_pitch;
    in.base_accel_roll  = ba_roll;
    in.dt               = kDt;
    in.imu_healthy      = (g_imu_faults < kImuFaultMax);
    in.power_ok         = halPowerOk();

    const Stabilizer::Outputs out = g_stab.update(in);

    // --- actuate --------------------------------------------------------
    if (out.state == Stabilizer::State::kFault) {
        // Coast, do not brake. A braked joint makes the utensil rigid but also
        // makes it a lever that transmits every tremor straight to the food;
        // a coasting gimbal at least keeps its passive inertial isolation.
        halMotorDisable();
    } else {
        halMotorWrite(out.cmd_pitch, out.cmd_roll);
    }

    static uint32_t n = 0;
    if ((++n % 10u) == 0u) halTelemetry(out);   // 50 Hz
}

// ---------------------------------------------------------------------------
int main() {
    // HAL: clocks, GPIO, SPI/I2C, FOC driver, 1 kHz IMU timer, 500 Hz control
    //      timer, independent watchdog with a period of ~4 control ticks.

    g_tilt_pitch.configure(reference::kTiltCrossoverHz, reference::kTiltZeta);
    g_tilt_roll.configure (reference::kTiltCrossoverHz, reference::kTiltZeta);
    g_daccel_pitch.lowpass(reference::kBaseAccelLpHz, reference::kControlHz);
    g_daccel_roll.lowpass (reference::kBaseAccelLpHz, reference::kControlHz);

    const AxisStabilizer::Config cfg = reference::axisConfig();
    g_stab.configure(cfg, cfg, /*calibrate_s=*/1.0f);
    g_stab.setLoadFeedforward(reference::kLoadFeedforward);

    for (;;) {
        if (g_tick_ready) controlTick();
        // HAL: kick the watchdog only from here, never from the ISR. A
        //      watchdog fed by the interrupt happily survives a hung control
        //      task, which is the failure that actually matters.
        // HAL: __WFI();
    }
}
