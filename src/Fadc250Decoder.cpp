#include "Fadc250Decoder.h"
#include <iostream>

using namespace fdec;

// --- byte-order aware reads -------------------------------------------------
// swap=false: native (LE on x86).  swap=true: byte-swapped (BE source, word-swapped by evio).
// When evio word-swaps a BE source, each 4-byte word gets reversed.
// For composite data this means:
//   c (1 byte): moves from word offset 0 to offset 3
//   i (4 bytes): correctly swapped as a whole word
//   l (8 bytes): each 4-byte half is swapped, but the halves may be in swapped order
//   s (2 bytes): within a word, bytes are reversed
//
// When swap=true, we undo the word-swap by reading each element in big-endian.

uint16_t Fadc250Decoder::rd16(const uint8_t *p, bool swap)
{
    if (!swap) { uint16_t v; std::memcpy(&v, p, 2); return v; }
    return (uint16_t(p[0]) << 8) | p[1];
}

uint32_t Fadc250Decoder::rd32(const uint8_t *p, bool swap)
{
    if (!swap) { uint32_t v; std::memcpy(&v, p, 4); return v; }
    return (uint32_t(p[0]) << 24) | (uint32_t(p[1]) << 16) | (uint32_t(p[2]) << 8) | p[3];
}

uint64_t Fadc250Decoder::rd64(const uint8_t *p, bool swap)
{
    if (!swap) { uint64_t v; std::memcpy(&v, p, 8); return v; }
    return (uint64_t(rd32(p, true)) << 32) | rd32(p + 4, true);
}

// --- auto-detect: try native first, if nchan > 16 try swapped --------------
std::vector<SlotData> Fadc250Decoder::Decode(const uint8_t *data, size_t nbytes)
{
    if (!data || nbytes < 17) return {};

    // try native (LE) first
    uint32_t nchan_le;
    std::memcpy(&nchan_le, data + 13, 4);           // offset 13 = after c(1)+i(4)+l(8)
    if (nchan_le > 0 && nchan_le <= 16)
        return Decode(data, nbytes, false);

    // try swapped (BE)
    uint32_t nchan_be = rd32(data + 13, true);
    if (nchan_be > 0 && nchan_be <= 16)
        return Decode(data, nbytes, true);

    std::cerr << "Fadc250Decoder: cannot determine byte order (nchan_le="
              << nchan_le << " nchan_be=" << nchan_be << ")\n";
    return {};
}

// --- decode with known byte order -------------------------------------------
std::vector<SlotData> Fadc250Decoder::Decode(const uint8_t *data, size_t nbytes, bool swap)
{
    std::vector<SlotData> slots;
    if (!data || !nbytes) return slots;

    size_t pos = 0;
    while (pos + 17 <= nbytes) {
        SlotData s;
        s.slot      = data[pos]; pos += 1;
        s.trigger   = static_cast<int32_t>(rd32(data + pos, swap)); pos += 4;
        s.timestamp = static_cast<int64_t>(rd64(data + pos, swap)); pos += 8;

        uint32_t nchan = rd32(data + pos, swap); pos += 4;
        if (nchan > 16) {
            std::cerr << "Fadc250Decoder: bad nchan=" << nchan
                      << " slot=" << (int)s.slot << " at pos=" << pos - 4 << "\n";
            break;
        }

        s.channels.resize(nchan);
        for (uint32_t i = 0; i < nchan; ++i) {
            if (pos + 5 > nbytes) break;
            s.channels[i].channel = data[pos]; pos += 1;

            uint32_t nsamp = rd32(data + pos, swap); pos += 4;
            if (nsamp > 4096) {
                std::cerr << "Fadc250Decoder: bad nsamp=" << nsamp << "\n";
                return slots;
            }

            size_t bytes = size_t(nsamp) * 2;
            if (pos + bytes > nbytes) { nsamp = (nbytes - pos) / 2; bytes = nsamp * 2; }

            s.channels[i].samples.resize(nsamp);
            if (!swap) {
                std::memcpy(s.channels[i].samples.data(), data + pos, bytes);
            } else {
                for (uint32_t j = 0; j < nsamp; ++j)
                    s.channels[i].samples[j] = rd16(data + pos + j * 2, true);
            }
            pos += bytes;
        }
        slots.push_back(std::move(s));
    }
    return slots;
}
