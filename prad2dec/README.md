# prad2dec

Static library for reading CODA evio data and analyzing FADC250 waveforms. Used by the PRad2 event viewer and online monitor, but also usable standalone.

## Components

**EvChannel** — Reads evio files via the C evio library. Calls `evRead` to get one buffer at a time, then `Scan()` builds a flat tree of bank nodes. `DecodeEvent()` extracts FADC250 composite data into pre-allocated arrays.

**EtChannel** — Subclass of EvChannel that reads from a running ET system instead of a file. Handles connection, station creation, event filtering, and byte-swap detection. The same `Scan()` / `DecodeEvent()` interface works on the received data.

**Fadc250Decoder** — Decodes the composite bank payload (format `c,i,l,N(c,Ns)`) into `SlotData` / `ChannelData` structs. Slots and channels are indexed by hardware number for direct access. No heap allocation — everything goes into the caller-provided `EventData`.

**WaveAnalyzer** — Waveform analysis on raw ADC samples:

- Triangular-kernel smoothing
- Iterative outlier-rejection pedestal (mean ± RMS)
- Local-maxima peak search with flat-plateau handling (overflow detection)
- Integration from peak outward, stopping at baseline crossing
- All scratch on stack, output into a fixed-size `WaveResult`

## Data Structures

`EventData` (~1.4 MB, allocate once and reuse):
```
EventData
  └─ rocs[MAX_ROCS=10]
       └─ slots[MAX_SLOTS=22]
            └─ channels[MAX_CHANNELS=16]
                 └─ uint16_t samples[MAX_SAMPLES=200]
```

`WaveResult`:
```
WaveResult
  ├─ ped: {mean, rms}
  ├─ npeaks
  └─ peaks[MAX_PEAKS=8]: {pos, time, height, integral, left, right, overflow}
```

## Usage

```cpp
#include "EvChannel.h"
#include "Fadc250Data.h"
#include "WaveAnalyzer.h"

evc::EvChannel ch;
ch.Open("data.evio");

fdec::EventData event;
fdec::WaveAnalyzer ana;
fdec::WaveResult wres;

while (ch.Read() == evc::status::success) {
    if (!ch.Scan()) continue;
    for (int i = 0; i < ch.GetNEvents(); ++i) {
        ch.DecodeEvent(i, event);
        for (int r = 0; r < event.nrocs; ++r) {
            auto &roc = event.rocs[r];
            for (int s = 0; s < fdec::MAX_SLOTS; ++s) {
                if (!roc.slots[s].present) continue;
                for (int c = 0; c < fdec::MAX_CHANNELS; ++c) {
                    if (!(roc.slots[s].channel_mask & (1u << c))) continue;
                    auto &cd = roc.slots[s].channels[c];
                    ana.Analyze(cd.samples, cd.nsamples, wres);
                    // wres.ped.mean, wres.peaks[0].integral, ...
                }
            }
        }
    }
}
```

## Evio Data Format

Composite bank tag `0xe101`, format string `c,i,l,N(c,Ns)`, packed native little-endian:

| Field | Type | Description |
|-------|------|-------------|
| `c` | uint8 | Slot number (3–20) |
| `i` | int32 | Event/trigger number |
| `l` | int64 | 48-bit timestamp |
| `N` | uint32 | Channels fired in this slot |
| `c` | uint8 | Channel number (0–15) |
| `N` | uint32 | Number of samples |
| `s` | int16[] | ADC values |

Multiple slots are packed back-to-back in one ROC's composite payload. ROC banks use tags `0x80`–`0x8c` (7 HyCal crates).

## Dependencies

- [evio](https://github.com/JeffersonLab/evio) (branch `evio-6.0`)
- [et](https://github.com/JeffersonLab/et) (for EtChannel)

Both can be fetched by CMake or linked from a prebuilt CODA installation.

## Files

```
include/
    EvStruct.h            Bank header parsers, EvNode, data type enums
    EvChannel.h           File reader + bank tree scanner + event decoder
    EtChannel.h           ET system reader (inherits EvChannel)
    EtConfigWrapper.h     C++ wrapper for ET config structs
    Fadc250Data.h         EventData / SlotData / ChannelData / Peak / Pedestal
    Fadc250Decoder.h      Composite payload decoder
    WaveAnalyzer.h        WaveConfig / WaveResult / WaveAnalyzer
src/
    EvChannel.cpp
    EtChannel.cpp
    Fadc250Decoder.cpp
    WaveAnalyzer.cpp
```
