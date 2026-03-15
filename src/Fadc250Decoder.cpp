#include "Fadc250Decoder.h"

using namespace fdec;

// Composite format "c,i,l,N(c,Ns)" — packed, native-endian.
// Multiple slots back-to-back until end of payload.

std::vector<SlotData> Fadc250Decoder::Decode(const uint8_t *data, size_t nbytes)
{
    std::vector<SlotData> slots;
    if (!data || !nbytes) return slots;

    size_t pos = 0;
    while (pos + 17 <= nbytes) {       // minimum: c(1)+i(4)+l(8)+N(4)
        SlotData s;
        s.slot = data[pos]; pos += 1;
        std::memcpy(&s.trigger,   data + pos, 4); pos += 4;
        std::memcpy(&s.timestamp, data + pos, 8); pos += 8;

        uint32_t nchan;
        std::memcpy(&nchan, data + pos, 4); pos += 4;
        s.channels.resize(nchan);

        for (uint32_t i = 0; i < nchan; ++i) {
            if (pos + 5 > nbytes) break;       // c(1)+N(4)
            s.channels[i].channel = data[pos]; pos += 1;

            uint32_t nsamp;
            std::memcpy(&nsamp, data + pos, 4); pos += 4;

            size_t bytes = static_cast<size_t>(nsamp) * 2;
            if (pos + bytes > nbytes) { nsamp = (nbytes - pos) / 2; bytes = nsamp * 2; }

            s.channels[i].samples.resize(nsamp);
            std::memcpy(s.channels[i].samples.data(), data + pos, bytes);
            pos += bytes;
        }
        slots.push_back(std::move(s));
    }
    return slots;
}
