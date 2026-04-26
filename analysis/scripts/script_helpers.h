#pragma once
//============================================================================
// script_helpers.h — small utility functions shared between analysis ACLiC
// scripts (gem_hycal_matching.C, plot_hits_at_hycal.C, …).
//
// Why a header instead of `static` helpers per-script:
//   Cling shares its dictionary scope across all ACLiC-loaded .C files in
//   the same ROOT session.  `static` / anonymous-namespace helpers in two
//   scripts therefore collide with `redefinition of …` errors at the
//   second `.L`.  Marking the helpers `inline` here gives them weak
//   external linkage so the dict-payload merge accepts them.
//
//   Each script just `#include "script_helpers.h"` and the symbols are
//   shared.  Add a new helper here whenever a second script needs it;
//   keep one-script-only helpers private to that script.
//============================================================================

#include "DaqConfig.h"

#include <nlohmann/json.hpp>

#include <cstdlib>
#include <fstream>
#include <map>
#include <regex>
#include <string>

// Resolve a possibly-relative database path to an absolute one using
// PRAD2_DATABASE_DIR.  Empty / already-absolute paths pass through.
inline std::string resolve_db_path(const std::string &p)
{
    if (p.empty()) return p;
    if (p[0] == '/' || p[0] == '\\') return p;
    if (p.size() >= 2 && p[1] == ':') return p;       // Windows drive letter
    const char *db = std::getenv("PRAD2_DATABASE_DIR");
    if (!db) return p;
    return std::string(db) + "/" + p;
}

// Sniff the run number out of a path like "prad_NNNNNN.evio.*".  Returns
// -1 if no plausible match is found, in which case LoadRunConfig falls
// back to the largest known runinfo entry (with a warning).
inline int extract_run_number_from_path(const std::string &path)
{
    static const std::regex pat(R"((?:prad|run)_0*(\d+))",
                                std::regex_constants::icase);
    std::smatch m;
    if (std::regex_search(path, m, pat)) {
        try { return std::stoi(m[1].str()); } catch (...) {}
    }
    return -1;
}

// Read database/config.json (under PRAD2_DATABASE_DIR or ./database) and
// return the resolved runinfo path, or "" if the pointer is missing /
// malformed.
inline std::string discover_runinfo_path()
{
    const char *db = std::getenv("PRAD2_DATABASE_DIR");
    std::string db_dir = db ? db : "database";
    std::ifstream f(db_dir + "/config.json");
    if (!f) return {};
    auto j = nlohmann::json::parse(f, nullptr, false, true);
    if (j.is_discarded() || !j.contains("runinfo") || !j["runinfo"].is_string())
        return {};
    return resolve_db_path(j["runinfo"].get<std::string>());
}

// EVIO bank-tag → logical-crate index (every ROC type).  Mirrors the
// roc_to_crate map in src/app_state_init.cpp / analysis/src/Replay.cpp.
inline std::map<int, int> build_full_crate_remap(const evc::DaqConfig &cfg)
{
    std::map<int, int> remap;
    for (const auto &re : cfg.roc_tags)
        remap[(int)re.tag] = re.crate;
    return remap;
}

// Same shape but only GEM ROCs — for GemSystem::LoadPedestals().
inline std::map<int, int> build_gem_crate_remap(const evc::DaqConfig &cfg)
{
    std::map<int, int> remap;
    for (const auto &re : cfg.roc_tags)
        if (re.type == "gem") remap[(int)re.tag] = re.crate;
    return remap;
}

// Strip the extension off a path so "out.pdf" becomes "out".  Used by
// scripts that derive a sibling .root output from a user-supplied
// canvas filename.  Leaves the directory alone.
inline std::string strip_extension(const std::string &p)
{
    auto dot = p.find_last_of('.');
    auto slash = p.find_last_of("/\\");
    if (dot == std::string::npos) return p;
    if (slash != std::string::npos && dot < slash) return p;
    return p.substr(0, dot);
}
