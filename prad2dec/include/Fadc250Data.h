#pragma once
//=============================================================================
// Fadc250Data.h — pre-allocated, flat event data for high-performance reading
//
// No heap allocation. Fixed-size arrays sized for worst case.
// Slot/channel indexed by hardware number for O(1) access.
//=============================================================================

#include <cstdint>
#include <cstring>

namespace fdec
{

// --- capacity limits (adjust to match your hardware) ------------------------
static constexpr int MAX_SAMPLES  = 200;   // samples per channel per event
static constexpr int MAX_CHANNELS = 64;    // channels per slot (16 for FADC250, 64 for ADC1881M)
static constexpr int MAX_SLOTS    = 32;    // slot IDs 0..31 (VME 3-20, Fastbus 0-25)
static constexpr int MAX_ROCS     = 10;    // number of ROC crates
static constexpr int MAX_PEAKS    = 8;     // max peaks per channel waveform

// --- per-channel data -------------------------------------------------------
struct ChannelData {
    int      nsamples;                     // 0 = channel not present
    uint16_t samples[MAX_SAMPLES];
};

// --- per-slot data ----------------------------------------------------------
struct SlotData {
    bool     present;                      // was this slot in the data?
    int32_t  trigger;                      // event number (from composite 'i')
    int64_t  timestamp;                    // 48-bit timestamp (from composite 'l')
    int      nchannels;                    // count of active channels
    uint64_t channel_mask;                 // bitmask: bit i set = channel i present
    ChannelData channels[MAX_CHANNELS];    // indexed by channel number

    void clear()
    {
        present = false;
        trigger = 0;
        timestamp = 0;
        nchannels = 0;
        channel_mask = 0;
        // only need to zero nsamples to mark channels empty
        for (int i = 0; i < MAX_CHANNELS; ++i)
            channels[i].nsamples = 0;
    }
};

// --- per-ROC data -----------------------------------------------------------
struct RocData {
    bool     present;
    uint32_t tag;                          // ROC bank tag (0x80, 0x84, etc.)
    int      nslots;                       // count of active slots
    SlotData slots[MAX_SLOTS];             // indexed by slot number

    void clear()
    {
        present = false;
        tag = 0;
        nslots = 0;
        for (int i = 0; i < MAX_SLOTS; ++i)
            slots[i].clear();
    }
};

// --- event-level information (extracted from TI bank + trigger bank) ---------
struct EventInfo {
    uint8_t  type;              // evc::EventType cast to uint8_t

    // --- two independent trigger fields (see database/trigger_bit.json) ------
    //
    // trigger_type: WHICH trigger caused this event (single, from TI event header)
    //   TI d[0] bits 31:24 = TS trigger table output (tiLoadTriggerTable(3))
    //   event_tag = 0x80 + trigger_type
    //   e.g. 0x29 = SSP raw sum, 0x30 = 100Hz pulser, 0x39 = LMS
    //
    // trigger_bits: WHAT signals were active at trigger time (multiple, from TI master d[5])
    //   32-bit FP input snapshot — multiple bits can fire simultaneously
    //   e.g. 0x01000100 = LMS (bit24) + SSP TRGBIT0 (bit8) both active
    //
    uint8_t  trigger_type;      // TI event_type — which trigger (single)
    uint32_t trigger_bits;      // FP trigger inputs — what fired (multi-bit)
    uint32_t event_tag;         // top-level bank tag (raw)
    int32_t  event_number;      // from trigger bank (0xC000) or TI
    int32_t  trigger_number;    // from TI data bank
    uint64_t timestamp;         // 48-bit TI timestamp (250MHz ticks)
    uint32_t run_number;        // from run info bank (0xE10F), 0 if absent
    uint32_t unix_time;         // from run info bank, 0 if absent

    void clear()
    {
        type = 0;
        trigger_type = 0;
        trigger_bits = 0;
        event_tag = 0;
        event_number = 0;
        trigger_number = 0;
        timestamp = 0;
        run_number = 0;
        unix_time = 0;
    }
};

// --- full event data --------------------------------------------------------
struct EventData {
    EventInfo info;                         // event-level metadata
    int      nrocs;                        // count of active ROCs
    int      roc_index[MAX_ROCS];          // maps i -> ROC index in rocs[]
    RocData  rocs[MAX_ROCS];

    void clear()
    {
        info.clear();
        nrocs = 0;
        for (int i = 0; i < MAX_ROCS; ++i)
            rocs[i].clear();
    }

    // find ROC by tag, returns nullptr if not found
    const RocData *findRoc(uint32_t tag) const
    {
        for (int i = 0; i < nrocs; ++i)
            if (rocs[i].tag == tag) return &rocs[i];
        return nullptr;
    }

    // convenience: iterate active slots in a ROC
    // usage: for (int s = 0; s < MAX_SLOTS; ++s) if (roc.slots[s].present) { ... }
};

// --- analysis helpers (to be filled later) ----------------------------------
struct Peak {
    float height, integral, time;
    int   pos, left, right;
    bool  overflow;
};

struct Pedestal {
    float mean, rms;
};

// --- firmware-faithful (Mode 1/2/3) analysis result -------------------------
//
// Quality bitmask — uint8_t so we have headroom past the firmware's 2-bit
// field.  Multiple flags can compose: e.g. peak at boundary AND NSA truncated.
//
// Q_GOOD              = 0       no flags
// Q_PEAK_AT_BOUNDARY  = 1 << 0  peak landed on the last sample
// Q_NSB_TRUNCATED     = 1 << 1  cross - NSB clipped at 0
// Q_NSA_TRUNCATED     = 1 << 2  cross + NSA clipped at N-1
// Q_VA_OUT_OF_RANGE   = 1 << 3  Va not bracketed by samples on the leading edge
//                                (very fast rise, or numerical edge case)
constexpr uint8_t Q_GOOD             = 0;
constexpr uint8_t Q_PEAK_AT_BOUNDARY = 1u << 0;
constexpr uint8_t Q_NSB_TRUNCATED    = 1u << 1;
constexpr uint8_t Q_NSA_TRUNCATED    = 1u << 2;
constexpr uint8_t Q_VA_OUT_OF_RANGE  = 1u << 3;

// One firmware-mode pulse (Mode 1 + Mode 2 + Mode 3 combined).
// All ADC values are pedestal-subtracted (s' = max(0, s - PED)).
struct DaqPeak {
    int      pulse_id;        // 0..MAX_PULSES-1
    float    vmin;            // = vnoise (manual step 1c)
    float    vpeak;           // peak ADC value (last sample before strict decrease)
    float    va;              // mid value = vmin + (vpeak - vmin) / 2
    int      coarse;          // 4-ns clock index of Vba (10-bit field in firmware)
    int      fine;            // 0..63 (6-bit firmware field)
    int      time_units;      // coarse * 64 + fine, LSB = 1/(CLK*64) = 62.5 ps
    float    time_ns;         // time_units * (CLK_NS / 64)
    int      cross_sample;    // first leading-edge sample > TET (Tcross)
    int      peak_sample;     // sample index where vpeak was found (i_peak)
    float    integral;        // Σ s' over [cross-NSB, cross+NSA] (Mode 2)
    int      window_lo;       // clamped Mode 1 window start (inclusive)
    int      window_hi;       // clamped Mode 1 window end   (inclusive)
    uint8_t  quality;         // bitmask of Q_* flags above
};

struct DaqWaveResult {
    float    vnoise;             // mean of first 4 pedestal-subtracted samples
    int      npeaks;
    DaqPeak  peaks[MAX_PEAKS];   // MAX_PEAKS=8 ≥ MAX_PULSES=4
    void clear() { vnoise = 0; npeaks = 0; }
};

} // namespace fdec
