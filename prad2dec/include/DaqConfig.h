#pragma once
//=============================================================================
// DaqConfig.h — configurable DAQ bank tags and event type identification
//
// All tags are configurable to accommodate DAQ format changes.
// Default values follow PRad-I conventions; update for PRad-II when confirmed.
//
// This is a plain struct — no JSON dependency. Loading from JSON is handled
// by the application layer (see load_daq_config() in DaqConfig.cpp).
//=============================================================================

#include <cstdint>
#include <vector>
#include <string>

namespace evc
{

struct DaqConfig
{
    // --- event type identification (top-level bank tag ranges) ---------------

    // physics event tag range (CODA built events)
    uint32_t physics_tag_min  = 0xFF50;
    uint32_t physics_tag_max  = 0xFF8F;

    // control event tags
    uint32_t prestart_tag     = 0x11;
    uint32_t go_tag           = 0x12;
    uint32_t end_tag          = 0x20;

    // sync event tag
    uint32_t sync_tag         = 0xC1;

    // EPICS slow control event tag
    uint32_t epics_tag        = 0x1F;

    // --- bank tags within physics events ------------------------------------

    // FADC250 composite data bank tag
    uint32_t fadc_composite_tag = 0xE101;

    // Trigger Interface (TI) data bank tag
    uint32_t ti_bank_tag        = 0xE10A;

    // EPICS data bank tag (within EPICS events)
    uint32_t epics_bank_tag     = 0xE114;

    // --- TI data format -----------------------------------------------------
    // Word offsets within TI bank data for timestamp extraction
    int ti_time_low_word  = 4;    // 32-bit low word of 48-bit timestamp
    int ti_time_high_word = 5;    // upper 16 bits of 48-bit timestamp
    int ti_time_high_mask = 0xFFFF;

    // --- ROC identification -------------------------------------------------
    // List of known ROC bank tags (populated from config)
    // If empty, all parent banks of FADC composite banks are accepted
    struct RocEntry {
        uint32_t    tag;
        std::string name;    // e.g. "HyCal_0", "HyCal_1"
    };
    std::vector<RocEntry> roc_tags;

    // --- helpers ------------------------------------------------------------
    bool is_physics(uint32_t tag) const
    {
        return tag >= physics_tag_min && tag <= physics_tag_max;
    }

    bool is_control(uint32_t tag) const
    {
        return tag == prestart_tag || tag == go_tag || tag == end_tag;
    }

    bool is_sync(uint32_t tag) const { return tag == sync_tag; }
    bool is_epics(uint32_t tag) const { return tag == epics_tag; }
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
    Control   = 7,   // other control events
};

inline EventType classify_event(uint32_t tag, const DaqConfig &cfg)
{
    if (cfg.is_physics(tag)) return EventType::Physics;
    if (cfg.is_sync(tag))    return EventType::Sync;
    if (cfg.is_epics(tag))   return EventType::Epics;
    if (tag == cfg.prestart_tag) return EventType::Prestart;
    if (tag == cfg.go_tag)       return EventType::Go;
    if (tag == cfg.end_tag)      return EventType::End;
    if (cfg.is_control(tag))     return EventType::Control;
    return EventType::Unknown;
}

} // namespace evc
