#pragma once
//=============================================================================
// DaqConfig.h — configurable DAQ bank tags and event type identification
//
// All tags are configurable to accommodate DAQ format changes.
// No defaults — a DAQ config JSON must be loaded before use.
//
// This is a plain struct — no JSON dependency. Loading from JSON is handled
// by the application layer (see load_daq_config.h).
//=============================================================================

#include <cstdint>
#include <vector>
#include <string>
#include <unordered_map>

namespace evc
{

struct DaqConfig
{
    // --- event type identification (top-level bank tag ranges) ---------------

    // physics event tags (JLab CODA convention)
    // Single-event mode: e.g. 0xB1, 0xB9, 0xFE (depends on CODA event writer)
    // Built-trigger mode: 0xFF50-0xFF8F (num = event count)
    std::vector<uint32_t> physics_tags;

    // monitoring event tags — physics events with TI data but no waveforms
    // (e.g. 100Hz pulser, scaler triggers). Classified as physics by is_physics()
    // but skipped in file viewer navigation and event iteration.
    std::vector<uint32_t> monitoring_tags;

    // control event tags (CODA2/JLab legacy — confirmed in PRad-II data)
    // CODA3 uses 0xFFD0-0xFFD4, recognized via is_control() range check
    uint32_t prestart_tag;
    uint32_t go_tag;
    uint32_t end_tag;

    // sync event tag
    uint32_t sync_tag;

    // EPICS slow control event tag
    uint32_t epics_tag;

    // --- ADC format selection ------------------------------------------------
    // "fadc250"  — FADC250 composite format c,i,l,N(c,Ns) (PRad-II)
    // "adc1881m" — Fastbus ADC1881M raw words (PRad)
    std::string adc_format;

    // Zero-suppression threshold for ADC1881M (in units of pedestal sigma).
    // Channels with (raw - ped_mean) < sparsify_sigma * ped_rms are suppressed.
    float sparsify_sigma;

    // --- bank tags within physics events ------------------------------------

    // FADC250 composite data bank tag (used when adc_format == "fadc250")
    uint32_t fadc_composite_tag;

    // ADC1881M raw data bank tag (used when adc_format == "adc1881m")
    uint32_t adc1881m_bank_tag;

    // Trigger Interface (TI) data bank tag (present in every ROC data block)
    uint32_t ti_bank_tag;

    // JLab event number/type bank (depth 1, single-event mode)
    uint32_t trigger_bank_tag;

    // Run info bank (in TI master crate only)
    uint32_t run_info_tag;

    // DAQ configuration readback string bank
    uint32_t daq_config_tag;

    // EPICS data bank tag (within EPICS events)
    uint32_t epics_bank_tag;

    // SSP/MPD raw data bank tags (GEM readout)
    // Multiple tags: 0xE10C (SSP trigger in TI master), 0x0DEA (VTP/MPD GEM data)
    std::vector<uint32_t> ssp_bank_tags;

    // FADC250 hardware-format raw data bank tag (0xE109, used when rol2 is skipped)
    uint32_t fadc_raw_tag = 0;

    // --- TI data format (fallback for single-event / non-CODA3 data) --------
    // TI bank layout: word[0]=header, word[1]=trigger#, word[2]=ts_low, word[3]=ts_high
    int ti_trigger_word;
    int ti_time_low_word;       // lower 32 bits of 48-bit timestamp
    int ti_time_high_word;      // upper bits of timestamp (shifted)
    uint32_t ti_time_high_mask;
    int      ti_time_high_shift;    // right-shift before combining

    // --- trigger bits extraction -----------------------------------------------
    // Two-layer extraction:
    //   1. Baseline (always): d[0] bits[31:24] = TI event_type = TS input pattern.
    //      tiLoadTriggerTable(3) encodes which TS inputs fired.
    //      event_tag = 0x80 + event_type (e.g. 0xB0 = monitoring, 0xA9 = physics).
    //   2. TI master override: d[trigger_type_word] with configured mask/shift.
    //      With tiSetFPInputReadout(1): d[5] = 32-bit FP trigger inputs.
    //      Bits 16-31 = v1495 triggers, bit 16 = LMS.
    int ti_trigger_type_word;
    int ti_trigger_type_shift;
    uint32_t ti_trigger_type_mask;

    // --- JLab trigger bank format (tag 0xC000, single-event mode) -----------
    // 3 words: event_number, event_tag, reserved
    int trig_event_number_word;
    int trig_event_type_word;

    // --- run info bank format (tag 0xE10F, in TI master crate) --------------
    int ri_run_number_word;
    int ri_event_count_word;
    int ri_unix_time_word;

    // --- ROC identification -------------------------------------------------
    struct RocEntry {
        uint32_t    tag;
        std::string name;
        int         crate = -1;
        std::string type;   // "fadc" (default), "gem", etc.
    };
    std::vector<RocEntry> roc_tags;

    // TI master crate tag (contains run info bank)
    uint32_t ti_master_tag;

    // --- diagnostics -----------------------------------------------------------
    bool verbose_decode = false;   // log unmatched bank tags in DecodeEvent

    // --- per-channel pedestals (ADC1881M) ------------------------------------
    struct PedEntry { float mean = 0.f; float rms = 0.f; };

    // pedestal lookup: packed key (crate<<32 | slot<<16 | channel) → PedEntry
    std::unordered_map<uint64_t, PedEntry> pedestals;

    static uint64_t pack_daq_key(int crate, int slot, int ch)
    {
        return (static_cast<uint64_t>(crate) << 32) |
               (static_cast<uint64_t>(slot)  << 16) |
               static_cast<uint64_t>(ch);
    }

    const PedEntry *get_pedestal(int crate, int slot, int ch) const
    {
        auto it = pedestals.find(pack_daq_key(crate, slot, ch));
        return (it != pedestals.end()) ? &it->second : nullptr;
    }

    // --- helpers ------------------------------------------------------------
    bool is_physics(uint32_t tag) const
    {
        // built-in trigger range for physics (0xFF50-0xFF8F, 0x00A0-0x00BF)
        if ((tag >= 0x00A0 && tag <= 0x00BF) || (tag >= 0xFF50 && tag <= 0xFF8F))
            return true;
        // single-event tags
        for (auto t : physics_tags)
            if (tag == t) return true;
        return false;
    }

    bool is_monitoring(uint32_t tag) const
    {
        for (auto t : monitoring_tags)
            if (tag == t) return true;
        return false;
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

    bool is_ssp_bank(uint32_t tag) const
    {
        for (auto t : ssp_bank_tags)
            if (t == tag) return true;
        return false;
    }

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
