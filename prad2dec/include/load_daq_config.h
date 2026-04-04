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

    // ADC format
    if (j.contains("adc_format"))
        cfg.adc_format = j["adc_format"].get<std::string>();

    // zero-suppression threshold (in sigma)
    if (j.contains("sparsify_sigma"))
        cfg.sparsify_sigma = j["sparsify_sigma"].get<float>();

    // bank tags
    if (j.contains("bank_tags")) {
        auto &bt = j["bank_tags"];
        if (bt.contains("fadc_composite")) cfg.fadc_composite_tag = parse_hex(bt["fadc_composite"]);
        if (bt.contains("adc1881m"))       cfg.adc1881m_bank_tag  = parse_hex(bt["adc1881m"]);
        if (bt.contains("ti_data"))        cfg.ti_bank_tag        = parse_hex(bt["ti_data"]);
        if (bt.contains("trigger_bank"))   cfg.trigger_bank_tag   = parse_hex(bt["trigger_bank"]);
        if (bt.contains("run_info"))       cfg.run_info_tag       = parse_hex(bt["run_info"]);
        if (bt.contains("daq_config"))     cfg.daq_config_tag     = parse_hex(bt["daq_config"]);
        if (bt.contains("epics_data"))     cfg.epics_bank_tag     = parse_hex(bt["epics_data"]);
        if (bt.contains("ssp_raw")) {
            cfg.ssp_bank_tags.clear();
            auto &v = bt["ssp_raw"];
            if (v.is_array()) {
                for (auto &item : v)
                    cfg.ssp_bank_tags.push_back(parse_hex(item));
            } else {
                cfg.ssp_bank_tags.push_back(parse_hex(v));
            }
        }
        if (bt.contains("fadc_raw"))       cfg.fadc_raw_tag       = parse_hex(bt["fadc_raw"]);
    }

    // TI format
    if (j.contains("ti_format")) {
        auto &ti = j["ti_format"];
        if (ti.contains("trigger_word"))        cfg.ti_trigger_word        = ti["trigger_word"].get<int>();
        if (ti.contains("time_low_word"))       cfg.ti_time_low_word       = ti["time_low_word"].get<int>();
        if (ti.contains("time_high_word"))      cfg.ti_time_high_word      = ti["time_high_word"].get<int>();
        if (ti.contains("time_high_mask"))      cfg.ti_time_high_mask      = parse_hex(ti["time_high_mask"]);
        if (ti.contains("time_high_shift"))     cfg.ti_time_high_shift     = ti["time_high_shift"].get<int>();
        if (ti.contains("trigger_type_word"))   cfg.ti_trigger_type_word   = ti["trigger_type_word"].get<int>();
        if (ti.contains("trigger_type_shift"))  cfg.ti_trigger_type_shift  = ti["trigger_type_shift"].get<int>();
        if (ti.contains("trigger_type_mask"))   cfg.ti_trigger_type_mask   = parse_hex(ti["trigger_type_mask"]);
    }

    // trigger bank format (0xC000)
    if (j.contains("trigger_bank")) {
        auto &tb = j["trigger_bank"];
        if (tb.contains("event_number_word")) cfg.trig_event_number_word = tb["event_number_word"].get<int>();
        if (tb.contains("event_type_word"))   cfg.trig_event_type_word   = tb["event_type_word"].get<int>();
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
            re.tag   = parse_hex(entry["tag"]);
            re.name  = entry.value("name", "");
            re.crate = entry.value("crate", -1);
            re.type  = entry.value("type", "");
            cfg.roc_tags.push_back(re);
        }
    }

    // TI master
    if (j.contains("ti_master_tag"))
        cfg.ti_master_tag = parse_hex(j["ti_master_tag"]);

    // pedestal file (loaded separately if specified)
    // Handled by the caller via load_pedestals()

    return true;
}

// Load per-channel pedestals from JSON file.
// Format: [{"crate":6,"slot":23,"channel":0,"mean":297.878,"rms":2.6972}, ...]
inline bool load_pedestals(const std::string &path, DaqConfig &cfg)
{
    std::ifstream f(path);
    if (!f.is_open()) {
        std::cerr << "load_pedestals: cannot open " << path << std::endl;
        return false;
    }

    nlohmann::json j;
    try { j = nlohmann::json::parse(f, nullptr, true, true); }
    catch (const nlohmann::json::parse_error &e) {
        std::cerr << "load_pedestals: parse error: " << e.what() << std::endl;
        return false;
    }

    cfg.pedestals.clear();
    for (auto &entry : j) {
        int crate   = entry.value("crate", 0);
        int slot    = entry.value("slot", 0);
        int channel = entry.value("channel", 0);
        float mean  = entry.value("mean", 0.f);
        float rms   = entry.value("rms", 0.f);
        cfg.pedestals[DaqConfig::pack_daq_key(crate, slot, channel)] = {mean, rms};
    }
    return true;
}

} // namespace evc
