//=============================================================================
// test_fadc_fw.cpp — hand-traced regression tests for Fadc250FwAnalyzer.
//
// Every expected value below was derived by hand from the manual algorithm
// (see docs/clas_fadc/FADC250UsersManual.pdf and FADC250_algorithms.md), then
// asserted against the implementation.  Test coverage:
//
//   T1: single trapezoid pulse, integer Vpeak, no truncation
//   T2: sub-threshold bump → 0 pulses
//   T3: two pulses in same window
//   T4: NSA truncation at end of window
//   T5: NSB truncation at start of window
//   T6: peak landing on the last sample (boundary peak)
//   T7: Va exactly == samples[k] → fine carries into coarse
//
// Build target: prad2dec_test_fadc_fw  (CTest: prad2dec.fadc_fw).
//=============================================================================

#include "Fadc250FwAnalyzer.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

using namespace fdec;

// ---------------------------------------------------------------------------
// Tiny test harness.  Stays self-contained so prad2dec keeps zero external
// test deps.  Each EXPECT_* records a failure but lets the test continue,
// so we see every divergence in one run.
// ---------------------------------------------------------------------------
namespace {
int g_total    = 0;
int g_failures = 0;

#define EXPECT_EQ(actual, expected) do {                                    \
    ++g_total;                                                              \
    auto _a = (actual); auto _e = (expected);                               \
    if (!(_a == _e)) {                                                      \
        ++g_failures;                                                       \
        std::fprintf(stderr,                                                \
            "  FAIL  %s:%d  %s == %s  (got %lld, want %lld)\n",             \
            __FILE__, __LINE__, #actual, #expected,                         \
            static_cast<long long>(_a), static_cast<long long>(_e));        \
    }                                                                       \
} while (0)

#define EXPECT_NEAR(actual, expected, tol) do {                             \
    ++g_total;                                                              \
    double _a = (double)(actual); double _e = (double)(expected);           \
    double _t = (double)(tol);                                              \
    if (std::fabs(_a - _e) > _t) {                                          \
        ++g_failures;                                                       \
        std::fprintf(stderr,                                                \
            "  FAIL  %s:%d  %s ≈ %s  (got %.6f, want %.6f, tol %.6g)\n",    \
            __FILE__, __LINE__, #actual, #expected, _a, _e, _t);            \
    }                                                                       \
} while (0)

void banner(const char *name)
{
    std::printf("\n[T] %s\n", name);
}

// Dispatch helper: run the analyzer on a vector with the given config + PED.
DaqWaveResult run(const std::vector<uint16_t> &raw, float PED,
                  const evc::DaqConfig::Fadc250FwConfig &cfg)
{
    Fadc250FwAnalyzer ana(cfg);
    DaqWaveResult res;
    ana.Analyze(raw.data(), static_cast<int>(raw.size()), PED, res);
    return res;
}

} // namespace

// ---------------------------------------------------------------------------
// T1 — Single trapezoid pulse, integer Vpeak, no truncation.
//
// Hand trace (PED=100, TET=50, NSB=4, NSA=10):
//   raw      = [100,100,100,100, 110,130,170,250,200,130,110,
//               100,100,100,100,100,100,100,100,100]
//   pedsub s = [  0,  0,  0,  0,  10, 30, 70,150,100, 30, 10,
//                 0,  0,  0,  0,  0,  0,  0,  0,  0]
//   Vnoise = mean(s[0..3]) = 0.  Vmin = 0.
//   Search starts i=4.  s[4]=10>0 → i_start=4.
//   Walk to peak: s[7]=150 is the last before s[8]=100<s[7] → i_peak=7,
//                  Vpeak=150.
//   Vpeak > TET ✓.  Tcross: first sample > TET on rising edge → s[6]=70 → cross=6.
//   Va = 0 + (150-0)/2 = 75.
//   Bracket: s[6]=70 < 75, s[7]=150 ≥ 75 → k=7.
//     Vba = 70, Vaa = 150, coarse = 6.
//   fine = round((75-70)/(150-70) * 64) = round(4.0) = 4.
//   time_units = 6*64 + 4 = 388.  time_ns = 388 * 4/64 = 24.25 ns.
//   Window = [6-4, 6+10] = [2, 16] (no truncation in n=20).
//   Σ s[2..16] = 0+0+10+30+70+150+100+30+10+0+0+0+0+0+0 = 400.
// ---------------------------------------------------------------------------
void test_T1_single_pulse()
{
    banner("T1 — single trapezoid pulse, integer Vpeak");

    std::vector<uint16_t> raw = {
        100,100,100,100, 110,130,170,250,200,130,110,
        100,100,100,100,100,100,100,100,100
    };
    evc::DaqConfig::Fadc250FwConfig cfg;   // TET=50, NSB=4, NSA=10, MAX_PULSES=4

    auto r = run(raw, 100.0f, cfg);

    EXPECT_EQ(r.npeaks, 1);
    EXPECT_NEAR(r.vnoise, 0.0f, 1e-6f);
    if (r.npeaks >= 1) {
        const auto &p = r.peaks[0];
        EXPECT_EQ(p.pulse_id, 0);
        EXPECT_NEAR(p.vmin, 0.0f, 1e-6f);
        EXPECT_NEAR(p.vpeak, 150.0f, 1e-6f);
        EXPECT_NEAR(p.va, 75.0f, 1e-6f);
        EXPECT_EQ(p.coarse, 6);
        EXPECT_EQ(p.fine, 4);
        EXPECT_EQ(p.time_units, 388);
        EXPECT_NEAR(p.time_ns, 24.25, 1e-6);
        EXPECT_EQ(p.cross_sample, 6);
        EXPECT_EQ(p.peak_sample, 7);
        EXPECT_NEAR(p.integral, 400.0f, 1e-4f);
        EXPECT_EQ(p.window_lo, 2);
        EXPECT_EQ(p.window_hi, 16);
        EXPECT_EQ(p.quality, Q_GOOD);
    }
}

// ---------------------------------------------------------------------------
// T2 — Sub-threshold bump.
//
// Vpeak after PED-subtract = 30 ≤ TET = 50 → no pulse reported.
// ---------------------------------------------------------------------------
void test_T2_sub_threshold()
{
    banner("T2 — sub-threshold bump → no pulses");

    std::vector<uint16_t> raw = {
        100,100,100,100, 105,110,130,110,105,100,100,100,100,100
    };
    evc::DaqConfig::Fadc250FwConfig cfg;   // TET=50

    auto r = run(raw, 100.0f, cfg);

    EXPECT_EQ(r.npeaks, 0);
}

// ---------------------------------------------------------------------------
// T3 — Two pulses in the same window.
//
// First pulse identical to T1.  Second pulse starts at sample 25.
//   pedsub s[25..30] = 10, 40, 100, 160, 80, 20
//   i_peak = 28, Vpeak = 160.
//   cross: s[27]=100 first > TET=50 → cross = 27.
//   Va = 80.  Bracket: s[27]=100 ≥ 80 → k=27.
//     Vba = s[26]=40, Vaa = s[27]=100, coarse = 26.
//   fine = round((80-40)/(100-40) * 64) = round(40/60 * 64) = round(42.666…) = 43.
//   time_units = 26*64 + 43 = 1707.  time_ns = 1707 * 4/64 = 106.6875 ns.
//   Window = [23, 37].  Σ = 0+0+10+40+100+160+80+20+5+0+0+0+0+0+0 = 415.
// ---------------------------------------------------------------------------
void test_T3_two_pulses()
{
    banner("T3 — two pulses in same window");

    std::vector<uint16_t> raw = {
        100,100,100,100, 110,130,170,250,200,130,110,
        100,100,100,100,100,100,100,100,100,100,100,100,100,100,
        110,140,200,260,180,120,105,
        100,100,100,100,100,100,100,100
    };
    evc::DaqConfig::Fadc250FwConfig cfg;

    auto r = run(raw, 100.0f, cfg);

    EXPECT_EQ(r.npeaks, 2);
    if (r.npeaks >= 2) {
        // First pulse — same as T1 (asserts only the differentiators).
        const auto &p0 = r.peaks[0];
        EXPECT_EQ(p0.pulse_id, 0);
        EXPECT_EQ(p0.cross_sample, 6);
        EXPECT_NEAR(p0.vpeak, 150.0f, 1e-6f);
        EXPECT_NEAR(p0.time_ns, 24.25, 1e-6);

        const auto &p1 = r.peaks[1];
        EXPECT_EQ(p1.pulse_id, 1);
        EXPECT_NEAR(p1.vpeak, 160.0f, 1e-6f);
        EXPECT_NEAR(p1.va, 80.0f, 1e-6f);
        EXPECT_EQ(p1.coarse, 26);
        EXPECT_EQ(p1.fine, 43);
        EXPECT_EQ(p1.time_units, 26 * 64 + 43);
        EXPECT_NEAR(p1.time_ns, 106.6875, 1e-6);
        EXPECT_EQ(p1.cross_sample, 27);
        EXPECT_NEAR(p1.integral, 415.0f, 1e-4f);
        EXPECT_EQ(p1.window_lo, 23);
        EXPECT_EQ(p1.window_hi, 37);
        EXPECT_EQ(p1.quality, Q_GOOD);
    }
}

// ---------------------------------------------------------------------------
// T4 — NSA truncation at end of window.
//
// n=14, pulse with cross=8, NSA=10 → cross+NSA=18 > 13 → truncated to 13.
//   pedsub s = [0,0,0,0,0,0,0,0, 70,150,100,30,10,0]
//   i_peak = 9, Vpeak = 150.  cross = 8.
//   Va = 75.  Bracket: s[8]=70 < 75, s[9]=150 ≥ 75 → k=9.  coarse=8, fine=4.
//   time_units = 8*64+4 = 516.  time_ns = 32.25.
//   Window = [4, 18] → clamped [4, 13].  Σ s[4..13] = 70+150+100+30+10 = 360.
// ---------------------------------------------------------------------------
void test_T4_nsa_truncated()
{
    banner("T4 — NSA truncation");

    std::vector<uint16_t> raw = {
        100,100,100,100,100,100,100,100, 170,250,200,130,110,100
    };
    evc::DaqConfig::Fadc250FwConfig cfg;   // NSA=10, NSB=4

    auto r = run(raw, 100.0f, cfg);

    EXPECT_EQ(r.npeaks, 1);
    if (r.npeaks >= 1) {
        const auto &p = r.peaks[0];
        EXPECT_EQ(p.cross_sample, 8);
        EXPECT_NEAR(p.vpeak, 150.0f, 1e-6f);
        EXPECT_EQ(p.coarse, 8);
        EXPECT_EQ(p.fine, 4);
        EXPECT_NEAR(p.time_ns, 32.25, 1e-6);
        EXPECT_EQ(p.window_lo, 4);
        EXPECT_EQ(p.window_hi, 13);
        EXPECT_NEAR(p.integral, 360.0f, 1e-4f);
        // Quality should at least set NSA-truncated.  Boundary flag also
        // legal (peak at sample 9, n=14 — i_peak ≠ n-1 here, so just NSA).
        EXPECT_EQ((int)(p.quality & Q_NSA_TRUNCATED), (int)Q_NSA_TRUNCATED);
        EXPECT_EQ((int)(p.quality & Q_NSB_TRUNCATED), 0);
        EXPECT_EQ((int)(p.quality & Q_PEAK_AT_BOUNDARY), 0);
    }
}

// ---------------------------------------------------------------------------
// T5 — NSB truncation at start of window.
//
// NSB = 6, pulse with cross at sample 4.  cross-NSB = -2 → clamp to 0.
//   pedsub s = [0,0,0,0, 70,150,100,30,10,0,0,…]
//   i_peak = 5, Vpeak = 150.  cross = 4.
//   Va = 75.  Bracket: s[4]=70<75, s[5]=150≥75 → k=5.  coarse=4, fine=4.
//   time_units = 4*64+4 = 260.  time_ns = 16.25.
//   Window = [-2, 14] → clamped [0, 14].  Σ s[0..14] = 70+150+100+30+10 = 360.
// ---------------------------------------------------------------------------
void test_T5_nsb_truncated()
{
    banner("T5 — NSB truncation");

    std::vector<uint16_t> raw = {
        100,100,100,100, 170,250,200,130,110,100,100,100,100,100,100,
        100,100,100,100,100
    };
    evc::DaqConfig::Fadc250FwConfig cfg;
    cfg.NSB = 6;
    cfg.NSA = 10;

    auto r = run(raw, 100.0f, cfg);

    EXPECT_EQ(r.npeaks, 1);
    if (r.npeaks >= 1) {
        const auto &p = r.peaks[0];
        EXPECT_EQ(p.cross_sample, 4);
        EXPECT_NEAR(p.vpeak, 150.0f, 1e-6f);
        EXPECT_EQ(p.coarse, 4);
        EXPECT_EQ(p.fine, 4);
        EXPECT_NEAR(p.time_ns, 16.25, 1e-6);
        EXPECT_EQ(p.window_lo, 0);
        EXPECT_EQ(p.window_hi, 14);
        EXPECT_NEAR(p.integral, 360.0f, 1e-4f);
        EXPECT_EQ((int)(p.quality & Q_NSB_TRUNCATED), (int)Q_NSB_TRUNCATED);
        EXPECT_EQ((int)(p.quality & Q_NSA_TRUNCATED), 0);
    }
}

// ---------------------------------------------------------------------------
// T6 — Boundary peak (i_peak == n-1).
//
// Monotonically rising waveform — peak walks to last sample.  Quality
// gets Q_PEAK_AT_BOUNDARY ∪ Q_NSA_TRUNCATED.
//
//   pedsub s = [0,0,0,0, 10,30,70,100,130,160]   (n=10)
//   i_peak = 9, Vpeak = 160.  cross: s[6]=70 first > 50 → cross = 6.
//   Va = 80.  Bracket: s[6]=70<80, s[7]=100≥80 → k=7.  coarse=6.
//   fine = round((80-70)/(100-70) * 64) = round(10/30*64) = round(21.333…) = 21.
//   time_units = 6*64+21 = 405.  time_ns = 405*0.0625 = 25.3125 ns.
//   Window = [2, 10] → clamped [2, 9].  Σ = 0+0+10+30+70+100+130+160 = 500.
// ---------------------------------------------------------------------------
void test_T6_boundary_peak()
{
    banner("T6 — boundary peak");

    std::vector<uint16_t> raw = {
        100,100,100,100, 110,130,170,200,230,260
    };
    evc::DaqConfig::Fadc250FwConfig cfg;
    cfg.NSB = 4;
    cfg.NSA = 4;

    auto r = run(raw, 100.0f, cfg);

    EXPECT_EQ(r.npeaks, 1);
    if (r.npeaks >= 1) {
        const auto &p = r.peaks[0];
        EXPECT_EQ(p.cross_sample, 6);
        EXPECT_NEAR(p.vpeak, 160.0f, 1e-6f);
        EXPECT_NEAR(p.va, 80.0f, 1e-6f);
        EXPECT_EQ(p.coarse, 6);
        EXPECT_EQ(p.fine, 21);
        EXPECT_EQ(p.time_units, 6 * 64 + 21);
        EXPECT_NEAR(p.time_ns, 25.3125, 1e-6);
        EXPECT_EQ(p.window_lo, 2);
        EXPECT_EQ(p.window_hi, 9);
        EXPECT_NEAR(p.integral, 500.0f, 1e-4f);
        EXPECT_EQ((int)(p.quality & Q_PEAK_AT_BOUNDARY), (int)Q_PEAK_AT_BOUNDARY);
        EXPECT_EQ((int)(p.quality & Q_NSA_TRUNCATED),    (int)Q_NSA_TRUNCATED);
    }
}

// ---------------------------------------------------------------------------
// T7 — Va exactly equal to samples[k] → fine carries into coarse.
//
//   pedsub s = [0,0,0,0, 20,40,80,160,80,40,20,0,…]
//   i_peak = 7, Vpeak = 160.  cross = 6.
//   Va = 80.  Bracket: s[6]=80 ≥ 80 → k=6.  Vba=s[5]=40, Vaa=s[6]=80.
//     coarse pre-carry = 5.
//   f = (80-40)/(80-40) * 64 = 64.0.  Carry: ++coarse → 6, fine = 0.
//   time_units = 6*64 + 0 = 384.  time_ns = 24.0 ns (lands exactly on s[6]).
// ---------------------------------------------------------------------------
void test_T7_va_on_sample()
{
    banner("T7 — Va exactly == samples[k]");

    std::vector<uint16_t> raw = {
        100,100,100,100, 120,140,180,260,180,140,120,100,100,100,100,100,
        100,100,100,100
    };
    evc::DaqConfig::Fadc250FwConfig cfg;

    auto r = run(raw, 100.0f, cfg);

    EXPECT_EQ(r.npeaks, 1);
    if (r.npeaks >= 1) {
        const auto &p = r.peaks[0];
        EXPECT_NEAR(p.vpeak, 160.0f, 1e-6f);
        EXPECT_NEAR(p.va, 80.0f, 1e-6f);
        EXPECT_EQ(p.coarse, 6);     // carried
        EXPECT_EQ(p.fine, 0);
        EXPECT_EQ(p.time_units, 384);
        EXPECT_NEAR(p.time_ns, 24.0, 1e-6);
        EXPECT_EQ(p.cross_sample, 6);
    }
}

// ---------------------------------------------------------------------------
int main()
{
    test_T1_single_pulse();
    test_T2_sub_threshold();
    test_T3_two_pulses();
    test_T4_nsa_truncated();
    test_T5_nsb_truncated();
    test_T6_boundary_peak();
    test_T7_va_on_sample();

    std::printf("\n[Summary] %d/%d assertions passed (%d failed)\n",
                g_total - g_failures, g_total, g_failures);
    return (g_failures == 0) ? 0 : 1;
}
