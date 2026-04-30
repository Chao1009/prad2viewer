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
        if (et.contains("monitoring")) {
            cfg.monitoring_tags.clear();
            for (auto &v : et["monitoring"])
                cfg.monitoring_tags.push_back(parse_hex(v));
        }
        if (et.contains("physics_base")) cfg.physics_base    = parse_hex(et["physics_base"]);
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
        if (bt.contains("tdc"))            cfg.tdc_bank_tag       = parse_hex(bt["tdc"]);
    }

    // DSC2 scaler bank (livetime measurement).  Section is optional; absent
    // or with bank_tag < 0 leaves the measurement disabled.
    if (j.contains("dsc_scaler")) {
        auto &ds = j["dsc_scaler"];
        if (ds.contains("bank_tag")) {
            auto &v = ds["bank_tag"];
            if (v.is_number_integer()) cfg.dsc_scaler.bank_tag = v.get<int>();
            else if (v.is_string()) {
                std::string s = v.get<std::string>();
                if (s.empty()) cfg.dsc_scaler.bank_tag = -1;
                else {
                    try { cfg.dsc_scaler.bank_tag = (int)std::stoul(s, nullptr, 0); }
                    catch (...) { cfg.dsc_scaler.bank_tag = -1; }
                }
            }
        }
        if (ds.contains("slot"))    cfg.dsc_scaler.slot    = ds["slot"].get<int>();
        if (ds.contains("channel")) cfg.dsc_scaler.channel = ds["channel"].get<int>();
        if (ds.contains("source")) {
            std::string src = ds["source"].get<std::string>();
            if      (src.empty())  { /* leave default */ }
            else if (src == "ref") cfg.dsc_scaler.source = DaqConfig::DscScaler::Source::Ref;
            else if (src == "trg") cfg.dsc_scaler.source = DaqConfig::DscScaler::Source::Trg;
            else if (src == "tdc") cfg.dsc_scaler.source = DaqConfig::DscScaler::Source::Tdc;
            else std::cerr << "load_daq_config: dsc_scaler.source '" << src
                           << "' unknown — expected 'ref' | 'trg' | 'tdc'\n";
        }
    }

    // FADC250 firmware-emulation parameters (TDC + Mode 1/2 reproduction).
    // Optional — absent leaves defaults; intended for DAQ signal studies and
    // for matching firmware-reported pulse values when running offline.
    // Field names follow the Hall-D V3 firmware API (faV3HallDSetProcMode);
    // NPEAK is accepted as an alias for MAX_PULSES.
    if (j.contains("fadc250_firmware")) {
        auto &ff = j["fadc250_firmware"];
        if (ff.contains("TET"))        cfg.fadc250_fw.TET        = ff["TET"].get<float>();
        if (ff.contains("NSB"))        cfg.fadc250_fw.NSB        = ff["NSB"].get<int>();
        if (ff.contains("NSA"))        cfg.fadc250_fw.NSA        = ff["NSA"].get<int>();
        if (ff.contains("MAX_PULSES")) cfg.fadc250_fw.MAX_PULSES = ff["MAX_PULSES"].get<int>();
        else if (ff.contains("NPEAK")) cfg.fadc250_fw.MAX_PULSES = ff["NPEAK"].get<int>();
        if (ff.contains("NSAT"))       cfg.fadc250_fw.NSAT       = ff["NSAT"].get<int>();
        if (ff.contains("NPED"))       cfg.fadc250_fw.NPED       = ff["NPED"].get<int>();
        if (ff.contains("MAXPED"))     cfg.fadc250_fw.MAXPED     = ff["MAXPED"].get<int>();
        if (ff.contains("CLK_NS"))     cfg.fadc250_fw.CLK_NS     = ff["CLK_NS"].get<float>();
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

    // sync / control-event absolute-time banks
    if (j.contains("sync_format")) {
        auto &sf = j["sync_format"];
        if (sf.contains("head_bank")) {
            auto &hb = sf["head_bank"];
            if (hb.contains("tag"))             cfg.sync_head_tag             = parse_hex(hb["tag"]);
            if (hb.contains("run_number_word")) cfg.sync_head_run_number_word = hb["run_number_word"].get<int>();
            if (hb.contains("counter_word"))    cfg.sync_head_counter_word    = hb["counter_word"].get<int>();
            if (hb.contains("unix_time_word"))  cfg.sync_head_unix_time_word  = hb["unix_time_word"].get<int>();
            if (hb.contains("event_tag_word"))  cfg.sync_head_event_tag_word  = hb["event_tag_word"].get<int>();
        }
        if (sf.contains("control_event")) {
            auto &ce = sf["control_event"];
            if (ce.contains("unix_time_word"))  cfg.sync_control_unix_time_word  = ce["unix_time_word"].get<int>();
            if (ce.contains("run_number_word")) cfg.sync_control_run_number_word = ce["run_number_word"].get<int>();
            if (ce.contains("run_type_word"))   cfg.sync_control_run_type_word   = ce["run_type_word"].get<int>();
        }
    }

    // bank structure: tag → { module, product, type }
    // This is the new authoritative source used by EvChannel's lazy accessors
    // to dispatch decoders by data product.  Absent entries fall back to the
    // legacy hard-coded dispatch in DecodeEvent — see EvChannel.cpp.
    if (j.contains("bank_structure")) {
        auto &bs = j["bank_structure"];
        if (bs.contains("data_banks")) {
            cfg.data_banks.clear();
            for (auto it = bs["data_banks"].begin(); it != bs["data_banks"].end(); ++it) {
                uint32_t tag;
                try { tag = static_cast<uint32_t>(std::stoul(it.key(), nullptr, 0)); }
                catch (...) {
                    std::cerr << "load_daq_config: bank_structure.data_banks key '"
                              << it.key() << "' is not a valid integer; skipping\n";
                    continue;
                }
                DaqConfig::DataBankInfo info;
                info.module  = it.value().value("module",  "");
                info.product = it.value().value("product", "");
                info.type    = it.value().value("type",    "");
                cfg.data_banks[tag] = std::move(info);
            }
        }
    }

    // companion-file pointers — application layer resolves them against the
    // database directory and decides which to actually load.  pedestal_file
    // is consumed by load_pedestals().
    if (j.contains("modules_file"))       cfg.modules_file       = j["modules_file"].get<std::string>();
    if (j.contains("hycal_daq_map_file")) cfg.hycal_daq_map_file = j["hycal_daq_map_file"].get<std::string>();
    if (j.contains("gem_daq_map_file"))   cfg.gem_daq_map_file   = j["gem_daq_map_file"].get<std::string>();
    if (j.contains("pedestal_file"))      cfg.pedestal_file      = j["pedestal_file"].get<std::string>();

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
