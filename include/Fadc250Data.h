#pragma once
//=============================================================================
// Fadc250Data.h — decoded FADC250 data structures
//=============================================================================

#include <vector>
#include <cstdint>
#include <utility>

namespace fdec
{

struct Peak {
    double height = 0, integral = 0, time = 0;
    uint32_t pos = 0, left = 0, right = 0;
    bool overflow = false;
};

struct Pedestal {
    double mean = 0, err = 0;
};

// raw waveform for one channel
struct ChannelData {
    uint8_t  channel = 0;
    Pedestal ped;
    std::vector<Peak>     peaks;
    std::vector<uint16_t> samples;
};

// one slot from composite format "c,i,l,N(c,Ns)"
struct SlotData {
    uint8_t  slot = 0;
    int32_t  trigger = 0;
    int64_t  timestamp = 0;
    std::vector<ChannelData> channels;
};

} // namespace fdec
