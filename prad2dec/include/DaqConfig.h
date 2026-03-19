#pragma once
//=============================================================================
// DaqConfig.h — configurable DAQ bank tags and event type identification
//
// All tags are configurable to accommodate DAQ format changes.
// Defaults match PRad-II data (prad_023109.evio format).
//
// This is a plain struct — no JSON dependency. Loading from JSON is handled
// by the application layer (see load_daq_config.h).
//=============================================================================

#include <cstdint>
#include <vector>
#include <string>

namespace evc
{

struct DaqConfig
{
    // --- event type identification (top-level bank tag ranges) ---------------

    // physics event tags (JLab CODA convention)
    // Single-event mode: 0xB1 or 0xFE (depends on CODA event writer)
    // Built-trigger mode: 0xFF50-0xFF8F (num = event count)
    std::vector<uint32_t> physics_tags = {0x00B1, 0x00FE};

    // control event tags (CODA2/JLab legacy — confirmed in PRad-II data)
    // CODA3 uses 0xFFD0-0xFFD4, recognized via is_control() range check
    uint32_t prestart_tag     = 0x11;
    uint32_t go_tag           = 0x12;
    uint32_t end_tag          = 0x14;

    // sync event tag
    uint32_t sync_tag         = 0xC1;

    // EPICS slow control event tag
    uint32_t epics_tag        = 0x1F;

    // --- bank tags within physics events ------------------------------------

    // FADC250 composite data bank tag
    uint32_t fadc_composite_tag = 0xE101;

    // Trigger Interface (TI) data bank tag (present in every ROC data block)
    uint32_t ti_bank_tag        = 0xE10A;

    // Run info bank (in TI master crate only)
    uint32_t run_info_tag       = 0xE10F;

    // DAQ configuration readback string bank
    uint32_t daq_config_tag     = 0xE10E;

    // EPICS data bank tag (within EPICS events)
    uint32_t epics_bank_tag     = 0xE114;

    // --- TI data format (fallback for single-event / non-CODA3 data) --------
    // TI bank layout: word[0]=header, word[1]=trigger#, word[2]=ts_low, word[3]=ts_high
    int ti_trigger_word   = 1;
    int ti_time_low_word  = 2;      // lower 32 bits of 48-bit timestamp
    int ti_time_high_word = 3;      // upper bits of timestamp (shifted)
    uint32_t ti_time_high_mask  = 0xFFFF0000;
    int      ti_time_high_shift = 16;   // right-shift before combining

    // --- trigger type extraction from TI bank --------------------------------
    // PRad-I: trigger bits in TI data[2] >> 24 (old format with 2 extra headers)
    // PRad-II: likely TI word[0] >> 24 (header encodes event/trigger type)
    // Set to -1 to disable trigger type extraction.
    int ti_trigger_type_word  = 0;
    int ti_trigger_type_shift = 24;
    uint32_t ti_trigger_type_mask = 0xFF;

    // --- run info bank format (tag 0xE10F, in TI master crate) --------------
    int ri_run_number_word     = 1;
    int ri_event_count_word    = 2;
    int ri_unix_time_word      = 3;

    // --- ROC identification -------------------------------------------------
    struct RocEntry {
        uint32_t    tag;
        std::string name;
    };
    std::vector<RocEntry> roc_tags;

    // TI master crate tag (contains run info bank)
    uint32_t ti_master_tag = 0x27;

    // --- helpers ------------------------------------------------------------
    bool is_physics(uint32_t tag) const
    {
        // single-event tags
        for (auto t : physics_tags)
            if (tag == t) return true;
        // built-trigger range (0xFF50-0xFF8F: PEB, SEB, streaming)
        return (tag >= 0xFF50 && tag <= 0xFF8F);
    }

    bool is_control(uint32_t tag) const
    {
        // CODA3 control event range (0xFFD0-0xFFD4)
        if (tag >= 0xFFD0 && tag <= 0xFFD4) return true;
        // configured tags (may be legacy CODA2: 0x11, 0x12, 0x14)
        return tag == prestart_tag || tag == go_tag || tag == end_tag;
    }

    bool is_sync(uint32_t tag) const
    {
        return tag == sync_tag || tag == 0xFFD0;
    }

    bool is_epics(uint32_t tag) const { return tag == epics_tag; }

    // CODA trigger bank identification (spec pages 21, 26, 31)
    // Built trigger bank: 0xFF20-0xFF2F (created by Event Builder)
    // Raw trigger bank:   0xFF10-0xFF1F (from ROC, before event building)
    static bool is_built_trigger_bank(uint32_t tag) { return tag >= 0xFF20 && tag <= 0xFF2F; }
    static bool is_raw_trigger_bank(uint32_t tag)   { return tag >= 0xFF10 && tag <= 0xFF1F; }
    static bool is_trigger_bank(uint32_t tag)       { return tag >= 0xFF10 && tag <= 0xFF4F; }

    // Trigger bank tag encodes what data is present (page 26):
    //   bit 0: has timestamps
    //   bit 1: has run number & run type
    //   bit 2: NO run-specific data (inverted)
    static bool trigger_bank_has_timestamps(uint32_t tag) { return (tag & 0x01) != 0; }
    static bool trigger_bank_has_run_info(uint32_t tag)   { return (tag & 0x02) != 0; }
};

// --- event type enum --------------------------------------------------------
enum class EventType : uint8_t {
    Unknown   = 0,
    Physics   = 1,
    Sync      = 2,
    Epics     = 3,
    Prestart  = 4,
    Go        = 5,
    End       = 6,
    Control   = 7,
};

inline EventType classify_event(uint32_t tag, const DaqConfig &cfg)
{
    // CODA3 control events (0xFFD0-0xFFD4)
    if (tag == 0xFFD1 || tag == cfg.prestart_tag) return EventType::Prestart;
    if (tag == 0xFFD2 || tag == cfg.go_tag)       return EventType::Go;
    if (tag == 0xFFD4 || tag == cfg.end_tag)      return EventType::End;
    if (cfg.is_sync(tag))                         return EventType::Sync;
    if (cfg.is_epics(tag))                        return EventType::Epics;
    if (cfg.is_physics(tag))                      return EventType::Physics;
    if (cfg.is_control(tag))                      return EventType::Control;
    return EventType::Unknown;
}

} // namespace evc
