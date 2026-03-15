#pragma once
//=============================================================================
// Fadc250Decoder.h — decode FADC250 data from composite evio banks
//
// Composite format "c,i,l,N(c,Ns)" — packed byte stream.
// The byte order depends on the evio library version and source endianness.
// Decode() auto-detects by trying both LE and BE for the first slot.
//=============================================================================

#include "Fadc250Data.h"
#include <vector>
#include <cstdint>
#include <cstring>

namespace fdec
{

class Fadc250Decoder
{
public:
    // Decode composite payload.
    // data/nbytes = inner bank payload (after tagseg + inner bank header).
    // Byte order is auto-detected from the first slot's fields.
    static std::vector<SlotData> Decode(const uint8_t *data, size_t nbytes);

    // Decode with explicit byte order (false = native/LE, true = swapped/BE)
    static std::vector<SlotData> Decode(const uint8_t *data, size_t nbytes, bool swap);

private:
    // read helpers for native and byte-swapped
    static uint16_t rd16(const uint8_t *p, bool swap);
    static uint32_t rd32(const uint8_t *p, bool swap);
    static uint64_t rd64(const uint8_t *p, bool swap);
};

} // namespace fdec
