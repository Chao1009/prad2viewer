#pragma once
//=============================================================================
// gain_factor.h — load per-module LMS gain factors from the database
//
// File naming convention:  <dir>/prad_XXXXXX_LMS.dat
// File format (whitespace-delimited, 7 columns):
//   Name  lms_peak  lms_sigma  lms_chi2/ndf  g1  g2  g3
//
// Only W (PWO) and G (PbGlass) module lines are read; LMS header lines and
// any other prefixes are ignored.
//
// Selection rule (same as RunInfoConfig):
//   run_num >= 0 -> file whose run number is the largest that is <= run_num
//   run_num <  0 -> file with the largest run number ("latest")
//
// Usage:
//   auto tbl = prad2::LoadGainFactors("/path/to/database/gain_factor", run);
//   float g1_W1   = tbl.w[1].g[0];
//   float g2_G123 = tbl.g[123].g[1];
//
//   auto corr = prad2::ComputeGainCorrection(db_dir + "/gain_factor", cur_run, ref_run);
//   float c_avg  = corr.w[id].avg;      // W模块id，三个LMS均值
//   float c_lms2 = corr.g[id].corr[1];  // G模块id，使用第2路LMS的correction
//   new_adc2mev = old_adc2mev * corr.w[id].avg

//=============================================================================

#include <algorithm>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <regex>
#include <string>

namespace prad2 {

// Single module's three gain factors (g1, g2, g3 from the LMS/alpha fit).
// Zero-initialised by default so missing entries are safe to use.
struct GainFactor {
    float g[3] = {0.f, 0.f, 0.f};  // g[0]=g1, g[1]=g2, g[2]=g3
};

// Full table for one run.  Arrays indexed by module numeric ID; index 0 is
// unused (modules are 1-based in the dat file).
struct GainFactorTable {
    static constexpr int MAX_W = 1157;  // W1 .. W1156
    static constexpr int MAX_G = 901;   // G1 .. G900

    GainFactor w[MAX_W];  // PWO crystal modules
    GainFactor g[MAX_G];  // PbGlass modules
    int  run_number = -1;
    bool loaded     = false;
};

// Returns the path to the best-matching gain factor file for run_num.
// Selection rule:
//   run_num >= 0 -> file whose run number is the largest that is <= run_num;
//                   if no such file exists, falls back to the nearest available file.
//   run_num <  0 -> file with the largest run number ("latest")
// Returns an empty string if the directory is empty or inaccessible.
inline std::string FindGainFactorFile(const std::string &dir, int run_num)
{
    const std::regex pat(R"(prad_(\d{6})_LMS\.dat)");

    int    best_run   = -1;
    std::string best_path;
    // fallback: nearest file when no file satisfies run <= run_num
    int    near_run   = -1;
    std::string near_path;

    std::error_code ec;
    for (auto &entry : std::filesystem::directory_iterator(dir, ec)) {
        if (!entry.is_regular_file()) continue;
        std::smatch m;
        std::string fname = entry.path().filename().string();
        if (!std::regex_match(fname, m, pat)) continue;

        int rn = std::stoi(m[1].str());
        if (run_num < 0) {
            if (rn > best_run) { best_run = rn; best_path = entry.path().string(); }
        } else if (rn <= run_num && rn > best_run) {
            best_run = rn; best_path = entry.path().string();
        }
        // Track nearest file regardless of direction (for fallback)
        if (near_run < 0 || std::abs(rn - run_num) < std::abs(near_run - run_num))
        { near_run = rn; near_path = entry.path().string(); }
    }
    if (ec) {
        std::cerr << "Warning: cannot iterate gain_factor dir " << dir
                  << ": " << ec.message() << "\n";
    }

    if (run_num >= 0 && best_path.empty() && !near_path.empty()) {
        std::cerr << "Warning: no gain factor file with run <= " << run_num
                  << " in " << dir << "; using nearest run " << near_run << " instead.\n";
        return near_path;
    }
    return best_path;
}

// Load gain factors for the given run number.
// On any failure returns a default-constructed (unloaded) table.
inline GainFactorTable LoadGainFactors(const std::string &dir, int run_num)
{
    GainFactorTable tbl;

    std::string path = FindGainFactorFile(dir, run_num);
    if (path.empty()) {
        std::cerr << "Warning: no gain factor file found in " << dir
                  << " for run " << run_num << "\n";
        return tbl;
    }

    std::ifstream f(path);
    if (!f) {
        std::cerr << "Warning: cannot open gain factor file " << path << "\n";
        return tbl;
    }

    // Parse the embedded run number from the file name.
    {
        std::regex pat(R"(prad_(\d{6})_LMS\.dat)");
        std::smatch m;
        std::string fname = std::filesystem::path(path).filename().string();
        if (std::regex_search(fname, m, pat))
            tbl.run_number = std::stoi(m[1].str());
    }

    std::string name;
    float col2, col3, col4, g1, g2, g3;
    while (f >> name >> col2 >> col3 >> col4 >> g1 >> g2 >> g3) {
        if (name.empty()) continue;

        if (name[0] == 'W') {
            int id = std::stoi(name.substr(1));
            if (id >= 1 && id < GainFactorTable::MAX_W) {
                tbl.w[id].g[0] = g1;
                tbl.w[id].g[1] = g2;
                tbl.w[id].g[2] = g3;
            }
        } else if (name[0] == 'G') {
            int id = std::stoi(name.substr(1));
            if (id >= 1 && id < GainFactorTable::MAX_G) {
                tbl.g[id].g[0] = g1;
                tbl.g[id].g[1] = g2;
                tbl.g[id].g[2] = g3;
            }
        }
        // LMS lines and anything else are silently skipped.
    }

    tbl.loaded = true;
    std::cerr << "GainFactor: loaded run " << tbl.run_number
              << " from " << path << "\n";
    return tbl;
}

// Per-module gain correction factor: correction[id] = g_ref / g_current.
// Applying this to the current ADC->MeV scale compensates for gain drift.
// A value of 1.0 means no correction needed; > 1.0 means gain dropped.
struct GainCorrTable {
    static constexpr int MAX_W = GainFactorTable::MAX_W;
    static constexpr int MAX_G = GainFactorTable::MAX_G;

    // correction[id][j] = ref.g[j] / cur.g[j]  (j = 0,1,2 for g1,g2,g3)
    // avg[id]           = mean of the valid (non-zero) per-LMS corrections
    struct Entry {
        float corr[3] = {1.f, 1.f, 1.f};  // per-LMS correction
        float avg      = 1.f;              // average over valid LMS channels
    };

    Entry w[MAX_W];
    Entry g[MAX_G];

    int ref_run = -1;
    int cur_run = -1;
};

// Compute the average of non-zero elements in a 3-element array.
// Returns 1.0 if no valid elements (safe fallback — no correction applied).
namespace detail {
inline float avgValid3(const float v[3])
{
    float sum = 0.f;
    int   n   = 0;
    for (int j = 0; j < 3; ++j) {
        if (v[j] != 0.f) { sum += v[j]; ++n; }
    }
    return n > 0 ? sum / static_cast<float>(n) : 1.f;
}
} // namespace detail

// Build a GainCorrTable from two already-loaded GainFactorTables.
// ref_tbl  — reference run (denominator in g_ref / g_cur, should be stable)
// cur_tbl  — current run (numerator target to be corrected)
// If either table is not loaded, returns a default (all-ones) table.
inline GainCorrTable ComputeGainCorrection(const GainFactorTable &ref_tbl,
                                           const GainFactorTable &cur_tbl)
{
    GainCorrTable corr;
    corr.ref_run = ref_tbl.run_number;
    corr.cur_run = cur_tbl.run_number;

    if (!ref_tbl.loaded || !cur_tbl.loaded) {
        std::cerr << "Warning: GainCorrTable: one or both tables not loaded, "
                     "returning identity correction.\n";
        return corr;
    }

    auto fill = [](GainCorrTable::Entry &e,
                   const GainFactor &ref, const GainFactor &cur)
    {
        float raw[3];
        for (int j = 0; j < 3; ++j) {
            if (ref.g[j] > 0.f && cur.g[j] > 0.f)
                raw[j] = ref.g[j] / cur.g[j];
            else
                raw[j] = 0.f;  // mark invalid (both ref and cur must be valid)
            e.corr[j] = (raw[j] > 0.f) ? raw[j] : 1.f;
        }
        e.avg = detail::avgValid3(raw);
    };

    for (int i = 1; i < GainCorrTable::MAX_W; ++i)
        fill(corr.w[i], ref_tbl.w[i], cur_tbl.w[i]);
    for (int i = 1; i < GainCorrTable::MAX_G; ++i)
        fill(corr.g[i], ref_tbl.g[i], cur_tbl.g[i]);

    std::cerr << "GainCorr  : ref=" << corr.ref_run
              << " cur=" << corr.cur_run << "\n";
    return corr;
}

// Convenience overload: load both tables from the database directory and
// compute the correction in one call.
inline GainCorrTable ComputeGainCorrection(const std::string &dir,
                                           int cur_run, int ref_run)
{
    auto cur_tbl = LoadGainFactors(dir, cur_run);
    auto ref_tbl = LoadGainFactors(dir, ref_run);
    return ComputeGainCorrection(ref_tbl, cur_tbl);
}

} // namespace prad2

