#pragma once
//=============================================================================
// ConfigSetup.h — analysis-side helpers around RunInfo configs.
//
// The RunConfig struct, LoadRunConfig() and WriteRunConfig() live in
// prad2det/include/RunInfoConfig.h so they can be reused by the viewer,
// Python bindings and ROOT scripts. This header keeps:
//   - the analysis-only gRunConfig global + backward-compat aliases
//   - TransformDetData() / RotateDetData() overloads (need PhysicsTools)
//   - get_run_str() / get_run_int() filename parsers
//=============================================================================

#include "PhysicsTools.h"
#include "RunInfoConfig.h"

#include <algorithm>
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

// Transform detector-frame coordinates to the target/beam-centered frame.
//
// Explicit-offset overloads (float beamX, float beamY, float ZfromTarget):
//   Always available; callers supply the numbers directly.
//
// RunConfig overloads (const RunConfig &geo = gRunConfig):
//   Use the geometry loaded by LoadRunConfig().
//   Default argument = gRunConfig, so existing single-arg calls like
//     TransformDetData(hc_hits);
//   still compile and use the global config unchanged.
//   Multi-run code can pass an explicit RunConfig:
//     TransformDetData(hc_hits, geo1);

// -- single-hit primitives (used internally by the vector overloads) ---------
inline void TransformDetData(HCHit &h, float detX, float detY, float ZfromTarget)
{
    h.x += detX;
    h.y += detY;
    h.z += ZfromTarget;
}
inline void TransformDetData(HCHit &h, const RunConfig &geo = gRunConfig)
{
    h.x += geo.hycal_x;
    h.y += geo.hycal_y;
    h.z += geo.hycal_z;
}
inline void TransformDetData(GEMHit &h, float detX, float detY, float ZfromTarget)
{
    h.x += detX;
    h.y += detY;
    h.z += ZfromTarget;
}
inline void TransformDetData(GEMHit &h, const RunConfig &geo = gRunConfig)
{
    int det_id = h.det_id;
    if (det_id >= 0 && det_id < 4) {
        h.x += geo.gem_x[det_id];
        h.y += geo.gem_y[det_id];
        h.z += geo.gem_z[det_id];
    }
}

// Apply successive rotations Rz → Ry → Rx (extrinsic, small-angle convention).
// Each angle is in degrees.  Only non-zero axes incur any computation cost.
inline void RotateDetData(HCHit &h, float x_deg, float y_deg, float z_deg)
{
    constexpr float kDeg2Rad = 3.14159265f / 180.f;
    float x = h.x, y = h.y, z = h.z;

    if (z_deg != 0.f) {
        float c = std::cos(z_deg * kDeg2Rad), s = std::sin(z_deg * kDeg2Rad);
        float nx = x * c - y * s;
        float ny = x * s + y * c;
        x = nx; y = ny;
    }
    if (y_deg != 0.f) {
        float c = std::cos(y_deg * kDeg2Rad), s = std::sin(y_deg * kDeg2Rad);
        float nx =  x * c + z * s;
        float nz = -x * s + z * c;
        x = nx; z = nz;
    }
    if (x_deg != 0.f) {
        float c = std::cos(x_deg * kDeg2Rad), s = std::sin(x_deg * kDeg2Rad);
        float ny = y * c - z * s;
        float nz = y * s + z * c;
        y = ny; z = nz;
    }
    h.x = x; h.y = y; h.z = z;
}
inline void RotateDetData(HCHit &h, const RunConfig &geo = gRunConfig)
{
    RotateDetData(h, geo.hycal_tilt_x, geo.hycal_tilt_y, geo.hycal_tilt_z);
}

inline void RotateDetData(GEMHit &h, float x_deg, float y_deg, float z_deg)
{
    constexpr float kDeg2Rad = 3.14159265f / 180.f;
    float x = h.x, y = h.y, z = h.z;

    if (z_deg != 0.f) {
        float c = std::cos(z_deg * kDeg2Rad), s = std::sin(z_deg * kDeg2Rad);
        float nx = x * c - y * s;
        float ny = x * s + y * c;
        x = nx; y = ny;
    }
    if (y_deg != 0.f) {
        float c = std::cos(y_deg * kDeg2Rad), s = std::sin(y_deg * kDeg2Rad);
        float nx =  x * c + z * s;
        float nz = -x * s + z * c;
        x = nx; z = nz;
    }
    if (x_deg != 0.f) {
        float c = std::cos(x_deg * kDeg2Rad), s = std::sin(x_deg * kDeg2Rad);
        float ny = y * c - z * s;
        float nz = y * s + z * c;
        y = ny; z = nz;
    }
    h.x = x; h.y = y; h.z = z;
}
inline void RotateDetData(GEMHit &h, const RunConfig &geo = gRunConfig)
{
    int det_id = h.det_id;
    if (det_id >= 0 && det_id < 4) {
        RotateDetData(h, geo.gem_tilt_x[det_id], geo.gem_tilt_y[det_id], geo.gem_tilt_z[det_id]);
    }
}

// -- HCHit vector ------------------------------------------------------------
inline void RotateDetData(std::vector<HCHit> &hc_hits,
                          float x_deg, float y_deg, float z_deg)
{
    for (auto &h : hc_hits) RotateDetData(h, x_deg, y_deg, z_deg);
}

// -- GEMHit vector -----------------------------------------------------------
inline void RotateDetData(std::vector<GEMHit> &gem_hits,
                          float x_deg, float y_deg, float z_deg)
{
    for (auto &h : gem_hits) RotateDetData(h, x_deg, y_deg, z_deg);
}

// -- RunConfig overloads (use tilting angles stored in config) -------------
inline void RotateDetData(std::vector<HCHit> &hc_hits,
                          const RunConfig &geo = gRunConfig)
{
    RotateDetData(hc_hits, geo.hycal_tilt_x, geo.hycal_tilt_y, geo.hycal_tilt_z);
}

inline void RotateDetData(std::vector<GEMHit> &gem_hits,
                          const RunConfig &geo = gRunConfig)
{
    for (auto &h : gem_hits) {
        int det_id = h.det_id;
        if (det_id >= 0 && det_id < 4) {
            RotateDetData(h, geo.gem_tilt_x[det_id],
                             geo.gem_tilt_y[det_id],
                             geo.gem_tilt_z[det_id]);
        }
    }
}

// -- HCHit vector ------------------------------------------------------------
inline void TransformDetData(std::vector<HCHit> &hc_hits, float detX, float detY, float ZfromTarget)
{
    for (auto &h : hc_hits) TransformDetData(h, detX, detY, ZfromTarget);
}

inline void TransformDetData(std::vector<HCHit> &hc_hits, const RunConfig &geo = gRunConfig)
{
    TransformDetData(hc_hits, geo.hycal_x, geo.hycal_y, geo.hycal_z);
}

// -- GEMHit vector -----------------------------------------------------------
inline void TransformDetData(std::vector<GEMHit> &gem_hits, float detX, float detY, float ZfromTarget)
{
    for (auto &h : gem_hits) TransformDetData(h, detX, detY, ZfromTarget);
}

// Each GEM hit is transformed using its own detector id.
inline void TransformDetData(std::vector<GEMHit> &gem_hits, const RunConfig &geo = gRunConfig)
{
    for (auto &h : gem_hits) {
        int det_id = h.det_id;
        if (det_id >= 0 && det_id < 4) {
            TransformDetData(h, geo.gem_x[det_id], geo.gem_y[det_id], geo.gem_z[det_id]);
        } else {
            std::cerr << "Warning: Invalid GEM det_id " << det_id << " for coordinate transformation\n";
        }
    }
}

// -- MollerData --------------------------------------------------------------
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
