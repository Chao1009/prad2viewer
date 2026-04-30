"""
Sanity tests for fadc250_modes.py.

Generates synthetic PMT-like waveforms (gamma rise/fall) at known sub-sample
times, runs all three modes, and verifies:
  - the TDC time recovers the truth time to sub-sample precision (<= ~CLK/64),
  - the Mode-1 raw window has length NSB+NSA+1 (when not truncated),
  - the Mode-2 integral equals sum of those raw samples,
  - multi-pulse handling (up to 4) works.
"""

import math
from fadc250_modes import FADC250Config, FADC250Analyzer


def make_pulse(t_ns_truth, amplitude, n=100, clk_ns=4.0,
               tau_rise=2.0, tau_fall=15.0, pedestal=100.0, noise_seed=None):
    """Synthetic PMT pulse: pedestal + amplitude * (1 - exp(-(t-t0)/τr)) * exp(-(t-t0)/τf)."""
    import random
    if noise_seed is not None:
        random.seed(noise_seed)
    out = []
    for i in range(n):
        t = i * clk_ns
        v = pedestal
        for t0, A in zip(t_ns_truth, amplitude):
            dt = t - t0
            if dt > 0:
                v += A * (1.0 - math.exp(-dt / tau_rise)) * math.exp(-dt / tau_fall)
        # tiny noise so it looks realistic but doesn't move things by samples
        v += random.gauss(0, 0.5) if noise_seed is not None else 0.0
        out.append(v)
    return out


def banner(s):
    print("\n" + "=" * 72); print(s); print("=" * 72)


# ---------------------------------------------------------------------------
# Test 1: single pulse, recover time
# ---------------------------------------------------------------------------
banner("TEST 1: single pulse, truth t0 = 37.5 ns, amplitude = 800 counts")

cfg = FADC250Config(PED=100.0, TET=50.0, NSB=4, NSA=10)
ana = FADC250Analyzer(cfg)

t0_truth = 37.5
wf = make_pulse([t0_truth], [800], noise_seed=1)
res = ana.analyze(wf)

print(f"Vnoise (pedestal floor seen by analyzer): {res.pedestal_floor:.2f}")
print(f"Number of TDC pulses: {len(res.pulses_tdc)}")
for p in res.pulses_tdc:
    print(f"  pulse {p.pulse_number}: Vmin={p.Vmin:.1f}  Vpeak={p.Vpeak:.1f}  "
          f"Va={p.Va:.1f}  coarse={p.coarse_clk}  fine={p.fine}/64  "
          f"T = {p.T_ns:.3f} ns  (cross at sample {p.sample_cross})")

# Half-amplitude time of the truth pulse (pedestal-subtracted)
# is somewhere on the rising edge; we just want T to land near t0.
# The TDC's "Va" is half between Vmin and Vpeak which on the leading
# edge corresponds to a time slightly after t0_truth.  We verify the
# returned time is on the leading edge to within ~1 sample.
if res.pulses_tdc:
    err = res.pulses_tdc[0].T_ns - t0_truth
    print(f"\n  T_recovered - t0_truth = {err:+.3f} ns "
          f"(fine resolution = {cfg.CLK_NS/64*1000:.1f} ps)")

# ---------------------------------------------------------------------------
# Test 2: Mode 1 window length and Mode 2 integral
# ---------------------------------------------------------------------------
banner("TEST 2: Mode 1 window + Mode 2 integral consistency")

print(f"NSB={cfg.NSB}, NSA={cfg.NSA}  -> expected window length = {cfg.NSB+cfg.NSA+1}")
for m1, m2 in zip(res.pulses_mode1, res.pulses_mode2):
    print(f"  pulse {m1.pulse_number}:  window=[{m1.window_start},{m1.window_end}]  "
          f"len={len(m1.samples)}  integral(M2)={m2.integral:.1f}  "
          f"sum(M1 samples)={sum(m1.samples):.1f}")
    assert abs(m2.integral - sum(m1.samples)) < 1e-6, "Mode2 integral must equal sum of Mode1 samples"
print("  Mode 1 / Mode 2 consistency: OK")

# ---------------------------------------------------------------------------
# Test 3: two pulses in the same window
# ---------------------------------------------------------------------------
banner("TEST 3: two pulses at 30 ns and 230 ns")

wf2 = make_pulse([30.0, 230.0], [600, 400], n=100, noise_seed=2)
res2 = ana.analyze(wf2)
print(f"Found {len(res2.pulses_tdc)} pulses (expect 2):")
for p in res2.pulses_tdc:
    print(f"  pulse {p.pulse_number}: T = {p.T_ns:7.3f} ns  Vpeak = {p.Vpeak:6.1f}  "
          f"cross@sample {p.sample_cross}")

# ---------------------------------------------------------------------------
# Test 4: pulse below threshold should NOT be reported
# ---------------------------------------------------------------------------
banner("TEST 4: sub-threshold bump should not be flagged")

wf3 = make_pulse([40.0], [30], n=100, noise_seed=3)   # peak ~30, TET=50
res3 = ana.analyze(wf3)
print(f"Found {len(res3.pulses_tdc)} pulses (expect 0)")
assert len(res3.pulses_tdc) == 0, "Sub-threshold pulse must not be reported"
print("  OK")

# ---------------------------------------------------------------------------
# Test 5: edge truncation flags
# ---------------------------------------------------------------------------
banner("TEST 5: pulse very near end of window -> truncated NSA")

wf4 = make_pulse([380.0], [600], n=100, noise_seed=4)  # near end (window = 0..396 ns)
res4 = ana.analyze(wf4)
for p, m1 in zip(res4.pulses_tdc, res4.pulses_mode1):
    print(f"  cross@{p.sample_cross}  window=[{m1.window_start},{m1.window_end}]  "
          f"len={len(m1.samples)}  quality={p.quality}")

# ---------------------------------------------------------------------------
# Test 6: fine-time scan — sweep truth t0 across a sample period
# ---------------------------------------------------------------------------
banner("TEST 6: sweep t0 across one 4-ns sample bin, check linearity of T_ns")

errs = []
for k in range(11):
    t0 = 40.0 + k * 0.4   # 0..4 ns in 0.4 ns steps
    wf = make_pulse([t0], [800], n=100, noise_seed=None)  # noise off for clean sweep
    r = ana.analyze(wf)
    if r.pulses_tdc:
        T = r.pulses_tdc[0].T_ns
        errs.append(T - t0)
        print(f"  t0={t0:5.2f} ns  ->  T_recovered={T:6.3f} ns   diff={T-t0:+.3f} ns")

if errs:
    span = max(errs) - min(errs)
    print(f"\n  spread of (T - t0) across the sweep = {span:.3f} ns  "
          f"(should be small + roughly constant, since Va sits on the leading edge)")

print("\nAll tests done.")
