#pragma once
//=============================================================================
// ConfigSetup.h — analysis-side helpers around RunInfo configs.
//
// The RunConfig struct, LoadRunConfig() and WriteRunConfig() live in
// prad2det/include/RunInfoConfig.h so they can be reused by the viewer,
// Python bindings and ROOT scripts. This header keeps:
//   - the analysis-only gRunConfig global + backward-compat aliases
//   - BuildLabTransforms() — bridge from RunConfig to DetectorTransform
//   - TransformDetData(MollerData, ...) for explicit calibration shifts
//   - get_run_str() / get_run_int() filename parsers
//
// "Detector-frame → lab-frame" is owned by prad2det's DetectorTransform —
// build the per-detector poses once via BuildLabTransforms(), then call
// xform.toLab(x, y[, z]) per hit.  The earlier RotateDetData/TransformDetData
// helpers duplicated that math and have been retired.
//=============================================================================

#include "DetectorTransform.h"
#include "PhysicsTools.h"
#include "RunInfoConfig.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cstdio>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

namespace analysis {

namespace fs = std::filesystem;

// Re-export the shared type so existing analysis code that says
// `analysis::RunConfig` keeps compiling without source changes.
using RunConfig = ::prad2::RunConfig;
using ::prad2::LoadRunConfig;
using ::prad2::WriteRunConfig;

// Global geometry config for single-run tools.
// Multi-run code should capture LoadRunConfig()'s return value into a
// local RunConfig instead of relying on this global.
inline RunConfig gRunConfig;

// HyCal + 4-GEM DetectorTransform bundle, ready for per-hit `toLab` calls.
// Build once per run (rotation matrices precomputed inside) and reuse.
struct LabTransforms {
    DetectorTransform hycal;
    std::array<DetectorTransform, 4> gem;
};

inline LabTransforms BuildLabTransforms(const RunConfig &geo = gRunConfig)
{
    LabTransforms t;
    t.hycal.set(geo.hycal_x, geo.hycal_y, geo.hycal_z,
                geo.hycal_tilt_x, geo.hycal_tilt_y, geo.hycal_tilt_z);
    for (int d = 0; d < 4; ++d) {
        t.gem[d].set(geo.gem_x[d], geo.gem_y[d], geo.gem_z[d],
                     geo.gem_tilt_x[d], geo.gem_tilt_y[d], geo.gem_tilt_z[d]);
    }
    return t;
}

// Apply a DetectorTransform to a hit in-place (lab = R*[x,y,z] + [tx,ty,tz]).
template <typename Hit>
inline void ApplyToLab(const DetectorTransform &xform, Hit &h)
{
    float lx, ly, lz;
    xform.toLab(h.x, h.y, h.z, lx, ly, lz);
    h.x = lx; h.y = ly; h.z = lz;
}

// MollerData is a translation-only calibration shift — used by det_calib to
// apply per-detector alignment offsets to fitted Moller pairs.  Not a
// detector-frame transform, so it stays as plain arithmetic.
inline void TransformDetData(MollerData &mollers, float detX, float detY, float ZfromTarget)
{
    for (auto &moller : mollers) {
        moller.first.x  += detX;
        moller.first.y  += detY;
        moller.first.z  += ZfromTarget;
        moller.second.x += detX;
        moller.second.y += detY;
        moller.second.z += ZfromTarget;
    }
}

// --- run number utilities ---------------------------------------------------
// Extract the run number embedded in a file name of the form
// ".../prad_<digits>...". Returns "unknown" / -1 on failure.
inline std::string get_run_str(const std::string &file_name)
{
    std::string fname = fs::path(file_name).filename().string();
    auto ppos = fname.find("prad_");
    if (ppos != std::string::npos) {
        size_t s = ppos + 5;
        size_t e = s;
        while (e < fname.size() && std::isdigit((unsigned char)fname[e])) e++;
        if (e > s) return std::to_string(std::stoul(fname.substr(s, e - s)));
    }
    std::cerr << "Warning: cannot extract run number from file name " << file_name << ", using 'unknown'.\n";
    return "unknown";
}

inline int get_run_int(const std::string &file_name)
{
    std::string run_str = get_run_str(file_name);
    if (run_str == "unknown") return -1;
    return std::stoi(run_str);
}

} // namespace analysis
