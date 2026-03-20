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
static constexpr int MAX_CHANNELS = 16;    // channels per FADC250 slot
static constexpr int MAX_SLOTS    = 22;    // slot IDs 0..21 (VME slots 3-20 typical)
static constexpr int MAX_ROCS     = 10;    // number of ROC crates

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
    uint32_t channel_mask;                 // bitmask: bit i set = channel i present
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
    uint32_t trigger_bits;      // 32 FP trigger input bits from TI master (word[5])
    uint32_t event_tag;         // top-level bank tag (raw)
    int32_t  event_number;      // from trigger bank (0xC000) or TI
    int32_t  trigger_number;    // from TI data bank
    uint64_t timestamp;         // 48-bit TI timestamp (250MHz ticks)
    uint32_t run_number;        // from run info bank (0xE10F), 0 if absent
    uint32_t unix_time;         // from run info bank, 0 if absent

    void clear()
    {
        type = 0;
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

} // namespace fdec
