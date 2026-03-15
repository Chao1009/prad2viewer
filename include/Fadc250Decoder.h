#pragma once
//=============================================================================
// Fadc250Decoder.h — decode FADC250 data from composite evio banks
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
    // Decode composite payload: format "c,i,l,N(c,Ns)"
    // data/nbytes = payload of the inner bank inside the composite envelope.
    // Returns one SlotData per FADC250 slot found in the payload.
    static std::vector<SlotData> Decode(const uint8_t *data, size_t nbytes);
};

} // namespace fdec
