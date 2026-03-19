#pragma once
//=============================================================================
// load_daq_config.h — load DaqConfig from JSON file
//
// Utility for applications. Requires nlohmann/json.
// Not part of prad2dec library (which has no JSON dependency).
//=============================================================================

#include "DaqConfig.h"
#include <nlohmann/json.hpp>
#include <fstream>
#include <iostream>
#include <string>

namespace evc
{

// parse hex string like "0xFF50" or plain integer to uint32_t
inline uint32_t parse_hex(const nlohmann::json &j)
{
    if (j.is_number()) return j.get<uint32_t>();
    std::string s = j.get<std::string>();
    return static_cast<uint32_t>(std::stoul(s, nullptr, 0));
}

inline bool load_daq_config(const std::string &path, DaqConfig &cfg)
{
    std::ifstream f(path);
    if (!f.is_open()) {
        std::cerr << "load_daq_config: cannot open " << path << std::endl;
        return false;
    }

    nlohmann::json j;
    try { j = nlohmann::json::parse(f, nullptr, true, true); } // allow comments
    catch (const nlohmann::json::parse_error &e) {
        std::cerr << "load_daq_config: parse error: " << e.what() << std::endl;
        return false;
    }

    // event tags
    if (j.contains("event_tags")) {
        auto &et = j["event_tags"];
        if (et.contains("physics")) {
            cfg.physics_tags.clear();
            for (auto &v : et["physics"])
                cfg.physics_tags.push_back(parse_hex(v));
        }
        if (et.contains("prestart"))     cfg.prestart_tag    = parse_hex(et["prestart"]);
        if (et.contains("go"))           cfg.go_tag          = parse_hex(et["go"]);
        if (et.contains("end"))          cfg.end_tag         = parse_hex(et["end"]);
        if (et.contains("sync"))         cfg.sync_tag        = parse_hex(et["sync"]);
        if (et.contains("epics"))        cfg.epics_tag       = parse_hex(et["epics"]);
    }

    // bank tags
    if (j.contains("bank_tags")) {
        auto &bt = j["bank_tags"];
        if (bt.contains("fadc_composite")) cfg.fadc_composite_tag = parse_hex(bt["fadc_composite"]);
        if (bt.contains("ti_data"))        cfg.ti_bank_tag        = parse_hex(bt["ti_data"]);
        if (bt.contains("run_info"))       cfg.run_info_tag       = parse_hex(bt["run_info"]);
        if (bt.contains("daq_config"))     cfg.daq_config_tag     = parse_hex(bt["daq_config"]);
        if (bt.contains("epics_data"))     cfg.epics_bank_tag     = parse_hex(bt["epics_data"]);
    }

    // TI format
    if (j.contains("ti_format")) {
        auto &ti = j["ti_format"];
        if (ti.contains("time_low_word"))   cfg.ti_time_low_word   = ti["time_low_word"].get<int>();
        if (ti.contains("time_high_word"))  cfg.ti_time_high_word  = ti["time_high_word"].get<int>();
        if (ti.contains("time_high_mask"))  cfg.ti_time_high_mask  = parse_hex(ti["time_high_mask"]);
        if (ti.contains("time_high_shift")) cfg.ti_time_high_shift = ti["time_high_shift"].get<int>();
    }

    // run info format
    if (j.contains("run_info")) {
        auto &ri = j["run_info"];
        if (ri.contains("run_number_word"))     cfg.ri_run_number_word   = ri["run_number_word"].get<int>();
        if (ri.contains("event_count_word"))    cfg.ri_event_count_word  = ri["event_count_word"].get<int>();
        if (ri.contains("unix_timestamp_word")) cfg.ri_unix_time_word    = ri["unix_timestamp_word"].get<int>();
    }

    // ROC tags
    if (j.contains("roc_tags")) {
        cfg.roc_tags.clear();
        for (auto &entry : j["roc_tags"]) {
            DaqConfig::RocEntry re;
            re.tag  = parse_hex(entry["tag"]);
            re.name = entry.value("name", "");
            cfg.roc_tags.push_back(re);
        }
    }

    // TI master
    if (j.contains("ti_master_tag"))
        cfg.ti_master_tag = parse_hex(j["ti_master_tag"]);

    return true;
}

} // namespace evc
