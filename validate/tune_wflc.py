#!/usr/bin/env python3
"""Isolate the tremor estimator and find gains that actually lock.

The failure mode being hunted: if the LMS weights adapt fast they absorb the
phase error by rotating the sin/cos pair, the residual collapses, and the
frequency loop never sees an error to descend. The model then FITS the signal
without KNOWING its frequency -- predict(0) looks fine and predict(t) is
garbage, which is exactly useless for latency compensation.
"""

import math
from validate import Biquad, Wflc, FS, DT, handle_angle, true_freq

DUR = 30.0
WARMUP = 6.0
N = int(DUR * FS)
I0 = int(WARMUP * FS)
ADV_MS = (0, 2, 4, 6, 8, 10, 12, 14, 16, 20)

# Precompute the plant-free truth once; the sweep only re-runs the estimator.
TH = [handle_angle(k * DT)[0] for k in range(N)]
TRUTH = {a: [handle_angle(k * DT + a / 1000.0)[2] for k in range(N)]
         for a in ADV_MS}
TR0 = [handle_angle(k * DT)[2] for k in range(N)]
DEN = sum(v * v for v in TR0[I0:])
FTRUE = [true_freq(k * DT) for k in range(N)]


def run(mu_w, mu_omega, M=2, front="bandpass", lo=3.0, hi=14.0):
    intent = Biquad().lowpass(1.2, FS)
    bhi = Biquad().highpass(lo, FS)
    blo = Biquad().lowpass(hi, FS)
    w = Wflc(M=M, mu_w=mu_w, mu_omega=mu_omega)
    intent.reset(TH[0])

    ferr = 0.0
    num = {a: 0.0 for a in ADV_MS}
    for k in range(N):
        hp = TH[k] - intent.process(TH[k])
        s = blo.process(bhi.process(hp)) if front == "bandpass" else hp
        w.update(s, DT)
        if k >= I0:
            ferr += abs(w.frequency_hz() - FTRUE[k])
            for a in ADV_MS:
                num[a] += (w.predict(a / 1000.0) - TRUTH[a][k]) ** 2
    ferr /= (N - I0)
    best_a, best_db = None, 9e9
    for a in ADV_MS:
        db = 20.0 * math.log10(math.sqrt(num[a] / DEN))
        if db < best_db:
            best_a, best_db = a, db
    return ferr, best_a, best_db


def main():
    print("Does the frequency loop lock? true tremor chirps 5.0 -> 9.0 Hz")
    print("'best dB' = residual after ideal cancellation, lower is better\n")
    print(f"{'front':>9}{'M':>3}{'mu_w':>7}{'mu_om':>8}"
          f"{'f err Hz':>10}{'best ms':>9}{'best dB':>9}")
    print("-" * 55)
    best = None
    for front in ("bandpass", "hponly"):
        for M in (1, 2):
            for mu_w in (0.3, 1.0, 3.0, 8.0):
                for mu_om in (100.0, 400.0, 1500.0, 5000.0):
                    fe, a, db = run(mu_w, mu_om, M, front)
                    mark = ""
                    if best is None or db < best[0]:
                        best = (db, front, M, mu_w, mu_om, fe, a)
                        mark = "  <-"
                    print(f"{front:>9}{M:>3}{mu_w:>7.1f}{mu_om:>8.0f}"
                          f"{fe:>10.3f}{a:>9d}{db:>9.1f}{mark}")
    db, front, M, mu_w, mu_om, fe, a = best
    print(f"\nBEST  front={front}  M={M}  mu_w={mu_w}  mu_omega={mu_om}")
    print(f"      freq err {fe:.3f} Hz, predict-ahead {a} ms, residual {db:.1f} dB")


if __name__ == "__main__":
    main()
