#pragma once
//=============================================================================
// ConfigSetup.h — detector geometry configuration for PRad2
//
// Loads beam energy and HyCal/GEM coordinates from a JSON config file,
// and provides helpers to transform hit coordinates from the detector
// frame to the target/beam-centered frame. Also parses run numbers from
// input file names.
// Depends on PhysicsTools (HCHit/GEMHit/MollerData) and nlohmann::json.
//=============================================================================

#include "PhysicsTools.h"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace analysis {

// --- detector geometry configuration struct ---------------------------------
// Holds all run-specific detector geometry and beam parameters.
// Using a struct allows multi-run processing without shared mutable state:
//
//   auto geo1 = LoadTransformConfig(path, run1);
//   auto geo2 = LoadTransformConfig(path, run2);
//   TransformDetData(hits1, geo1);
//   TransformDetData(hits2, geo2);
//
struct GeoConfig {
    float Ebeam   = 0.f;
    float hycal_z = 6222.5f;
    float hycal_x = 0.f;
    float hycal_y = 0.f;
    float gem_z[4] = {5423.71f, 5384.00f, 5823.71f, 5784.00f};
    float gem_x[4] = {0.f, 0.f, 0.f, 0.f};
    float gem_y[4] = {0.f, 0.f, 0.f, 0.f};
};

// Global geometry config for single-run tools.
// Multi-run code should capture LoadTransformConfig()'s return value into a
// local GeoConfig instead of relying on this global.
inline GeoConfig gGeoConfig;

// Backward-compatible aliases pointing into gGeoConfig (zero overhead).
// Existing code that reads/writes Ebeam_, hycal_x_, gem_z_[] etc. continues
// to work without any modification.
inline float        &Ebeam_   = gGeoConfig.Ebeam;
inline float        &hycal_z_ = gGeoConfig.hycal_z;
inline float        &hycal_x_ = gGeoConfig.hycal_x;
inline float        &hycal_y_ = gGeoConfig.hycal_y;
inline float *const  gem_z_   = gGeoConfig.gem_z;
inline float *const  gem_x_   = gGeoConfig.gem_x;
inline float *const  gem_y_   = gGeoConfig.gem_y;

// --- config file loading ----------------------------------------------------
// Returns a GeoConfig populated from the best-matching entry in the JSON file.
// Selects the entry whose run_number is the largest value <= run_num.
// If run_num < 0 (unknown), uses the entry with the largest run_number.
//
// Single-run tools:  gGeoConfig = LoadTransformConfig(path, run);
// Multi-run tools:   auto geo1 = LoadTransformConfig(path, run1);
//                    auto geo2 = LoadTransformConfig(path, run2);
inline GeoConfig LoadTransformConfig(const std::string &transform_config, int run_num)
{
    GeoConfig result;   // start from defaults defined in GeoConfig

    std::ifstream cfg_f(transform_config);
    if (!cfg_f) {
        std::cerr << "Warning: cannot open config file " << transform_config << ", using defaults.\n";
        std::cerr << "Warning: Ebeam not set, may cause something wrong\n";
        return result;
    }
    auto cfg = nlohmann::json::parse(cfg_f, nullptr, false, true);
    if (cfg.is_discarded()) {
        std::cerr << "Warning: failed to parse " << transform_config << ", using defaults.\n";
        std::cerr << "Warning: Ebeam not set, may cause something wrong\n";
        return result;
    }
    if (!cfg.contains("configurations") || !cfg["configurations"].is_array()) {
        std::cerr << "Warning: " << transform_config << " has no \"configurations\" array, using defaults.\n";
        return result;
    }

    // find the best-matching entry
    const nlohmann::json *best = nullptr;
    int best_run = -1;

    if (run_num < 0) {
        std::cerr << "Warning: unknown run number, using the entry with the largest run_number.\n";
    }
    for (const auto &entry : cfg["configurations"]) {
        if (!entry.contains("run_number")) continue;
        int rn = entry["run_number"].get<int>();
        if (run_num < 0) {
            if (rn > best_run) { best = &entry; best_run = rn; }
        } else {
            if (rn <= run_num && rn > best_run) { best = &entry; best_run = rn; }
        }
    }

    if (best == nullptr) {
        std::cerr << "Warning: no matching configuration found in " << transform_config
                  << " for run " << run_num << ", using defaults.\n";
        return result;
    }

    const auto &c = *best;
    if (c.contains("Ebeam")) result.Ebeam = c["Ebeam"].get<float>();
    if (c.contains("hycal")) {
        const auto &h = c["hycal"];
        if (h.contains("z")) result.hycal_z = h["z"].get<float>();
        if (h.contains("x")) result.hycal_x = h["x"].get<float>();
        if (h.contains("y")) result.hycal_y = h["y"].get<float>();
    }
    if (c.contains("gem")) {
        const auto &g = c["gem"];
        auto load_gem_array = [&](const char *key, float (&dst)[4]) {
            if (!g.contains(key)) return;
            const auto &arr = g[key];
            if (!arr.is_array() || arr.size() != 4) {
                std::cerr << "Warning: gem." << key << " must be an array of 4 numbers, ignoring.\n";
                return;
            }
            for (int d = 0; d < 4; d++) dst[d] = arr[d].get<float>();
        };
        load_gem_array("z", result.gem_z);
        load_gem_array("x", result.gem_x);
        load_gem_array("y", result.gem_y);
    }
    std::cerr << "Loaded detector coordinates config (run_number=" << best_run
              << ") from: " << transform_config << "\n";
    return result;
}

// --- config file writing ----------------------------------------------------
// Appends a new entry (run_number + GeoConfig) to the "configurations" array
// in the given JSON file. If the file does not exist, it is created from
// scratch. If an entry with the same run_number already exists, it is
// overwritten in-place. The updated JSON is written back atomically via a
// temporary file to avoid corruption on failure.
inline bool WriteTransformConfig(const std::string &transform_config, int run_num,
                                 const GeoConfig &geo)
{
    // --- load existing file (or start empty) --------------------------------
    nlohmann::json cfg;
    {
        std::ifstream cfg_f(transform_config);
        if (cfg_f) {
            cfg = nlohmann::json::parse(cfg_f, nullptr, false, true);
            if (cfg.is_discarded()) {
                std::cerr << "Warning: failed to parse " << transform_config
                          << ", will overwrite with new data.\n";
                cfg = nlohmann::json::object();
            }
        } else {
            cfg = nlohmann::json::object();
        }
    }

    // ensure top-level structure
    if (!cfg.contains("information"))
        cfg["information"] = "z positions(mm) from target center, x/y (mm) from beam center";
    if (!cfg.contains("units"))
        cfg["units"] = "Ebeam is in MeV, z in mm, x/y in mm";
    if (!cfg.contains("configurations") || !cfg["configurations"].is_array())
        cfg["configurations"] = nlohmann::json::array();

    // --- build new entry ----------------------------------------------------
    nlohmann::json entry;
    entry["run_number"] = run_num;
    entry["Ebeam"]      = geo.Ebeam;
    entry["hycal"]["z"] = geo.hycal_z;
    entry["hycal"]["x"] = geo.hycal_x;
    entry["hycal"]["y"] = geo.hycal_y;
    entry["gem"]["z"]   = nlohmann::json::array({geo.gem_z[0], geo.gem_z[1], geo.gem_z[2], geo.gem_z[3]});
    entry["gem"]["x"]   = nlohmann::json::array({geo.gem_x[0], geo.gem_x[1], geo.gem_x[2], geo.gem_x[3]});
    entry["gem"]["y"]   = nlohmann::json::array({geo.gem_y[0], geo.gem_y[1], geo.gem_y[2], geo.gem_y[3]});

    // --- replace existing entry or append -----------------------------------
    auto &arr = cfg["configurations"];
    bool replaced = false;
    for (auto &e : arr) {
        if (e.contains("run_number") && e["run_number"].get<int>() == run_num) {
            e = entry;
            replaced = true;
            break;
        }
    }
    if (!replaced) arr.push_back(entry);

    // sort by run_number for readability
    std::sort(arr.begin(), arr.end(), [](const nlohmann::json &a, const nlohmann::json &b) {
        int ra = a.contains("run_number") ? a["run_number"].get<int>() : -1;
        int rb = b.contains("run_number") ? b["run_number"].get<int>() : -1;
        return ra < rb;
    });

    // --- write back atomically via a temporary file -------------------------
    std::string tmp_path = transform_config + ".tmp";
    {
        std::ofstream out(tmp_path);
        if (!out) {
            std::cerr << "Error: cannot write to " << tmp_path << "\n";
            return false;
        }
        out << cfg.dump(4) << "\n";
    }
    if (std::rename(tmp_path.c_str(), transform_config.c_str()) != 0) {
        std::cerr << "Error: failed to rename " << tmp_path << " -> " << transform_config << "\n";
        return false;
    }
    std::cerr << (replaced ? "Updated" : "Appended") << " run_number=" << run_num
              << " in " << transform_config << "\n";
    return true;
}


// Transform detector-frame coordinates to the target/beam-centered frame.
//
// Explicit-offset overloads (float beamX, float beamY, float ZfromTarget):
//   Always available; callers supply the numbers directly.
//
// GeoConfig overloads (const GeoConfig &geo = gGeoConfig):
//   Use the geometry loaded by LoadTransformConfig().
//   Default argument = gGeoConfig, so existing single-arg calls like
//     TransformDetData(hc_hits);
//   still compile and use the global config unchanged.
//   Multi-run code can pass an explicit GeoConfig:
//     TransformDetData(hc_hits, geo1);

// -- single-hit primitives (used internally by the vector overloads) ---------
inline void TransformDetData(HCHit &h, float beamX, float beamY, float ZfromTarget)
{
    h.x -= beamX;
    h.y -= beamY;
    h.z += ZfromTarget;
}

inline void TransformDetData(GEMHit &h, float beamX, float beamY, float ZfromTarget)
{
    h.x -= beamX;
    h.y -= beamY;
    h.z += ZfromTarget;
}

// -- HCHit vector ------------------------------------------------------------
inline void TransformDetData(std::vector<HCHit> &hc_hits, float beamX, float beamY, float ZfromTarget)
{
    for (auto &h : hc_hits) TransformDetData(h, beamX, beamY, ZfromTarget);
}

inline void TransformDetData(std::vector<HCHit> &hc_hits, const GeoConfig &geo = gGeoConfig)
{
    TransformDetData(hc_hits, geo.hycal_x, geo.hycal_y, geo.hycal_z);
}

// -- GEMHit vector -----------------------------------------------------------
inline void TransformDetData(std::vector<GEMHit> &gem_hits, float beamX, float beamY, float ZfromTarget)
{
    for (auto &h : gem_hits) TransformDetData(h, beamX, beamY, ZfromTarget);
}

// Each GEM hit is transformed using its own detector id.
inline void TransformDetData(std::vector<GEMHit> &gem_hits, const GeoConfig &geo = gGeoConfig)
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
inline void TransformDetData(MollerData &mollers, float beamX, float beamY, float ZfromTarget)
{
    for (auto &moller : mollers) {
        moller.first.x  -= beamX;
        moller.first.y  -= beamY;
        moller.first.z  += ZfromTarget;
        moller.second.x -= beamX;
        moller.second.y -= beamY;
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
