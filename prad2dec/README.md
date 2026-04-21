# prad2dec

Static library for reading CODA EVIO data and decoding detector electronics.

## Components

- **EvChannel** — EVIO file reader with bank tree scanning, lazy per-product
  decoding, sequential + random-access modes
- **EtChannel** — Subclass for live ET system reading (`-DWITH_ET=ON`)
- **Fadc250Decoder** — FADC250 composite waveform decoding
- **Fadc250RawDecoder** — FADC250 raw hardware format decoding
- **Adc1881mDecoder** — ADC1881M decoding (original PRad)
- **SspDecoder** — SSP/MPD fiber data for GEM readout
- **WaveAnalyzer** — Pedestal subtraction, peak search, integration

## Usage

### Sequential read

```cpp
#include "EvChannel.h"

evc::EvChannel ch;
ch.SetConfig(cfg);                    // evc::DaqConfig
ch.Open("data.evio");

while (ch.Read() == evc::status::success) {
    if (!ch.Scan()) continue;
    if (ch.GetEventType() != evc::EventType::Physics) continue;

    for (int i = 0; i < ch.GetNEvents(); ++i) {
        ch.SelectEvent(i);                   // picks sub-event, clears cache
        const auto &info = ch.Info();        // cheapest: TI/trigger metadata
        const auto &fadc = ch.Fadc();        // decoded on first call, then cached
        // also available: ch.Gem(), ch.Tdc(), ch.Vtp()
    }
}
ch.Close();
```

`Info()` / `Fadc()` / `Gem()` / `Tdc()` / `Vtp()` each decode on the first
call after `SelectEvent()` and return a cached reference on subsequent
calls for the same sub-event, so asking for `Fadc()` then `Gem()` on one
event costs nothing extra if you only need one of them.

### Random access (evio "ra" mode)

`OpenRandomAccess` mmaps the file and has evio build an internal
event-pointer table during open — subsequent `ReadEventByIndex(i)`
calls jump to any event in O(1). No close/reopen needed to go
backwards.

```cpp
evc::EvChannel ch;
ch.SetConfig(cfg);
ch.OpenRandomAccess("data.evio");

int n = ch.GetRandomAccessEventCount();      // total evio events in the file
for (int i : {0, n / 2, n - 1}) {
    if (ch.ReadEventByIndex(i) != evc::status::success) continue;
    if (!ch.Scan()) continue;
    ch.SelectEvent(0);
    auto &fadc = ch.Fadc();
    // …
}
```

Random-access events are **evio events** (blocks); a CODA built-trigger
block can contain multiple physics sub-events, so after
`ReadEventByIndex` + `Scan` use `GetNEvents()` + `SelectEvent(i)` as
usual. See `src/evio_data_source.cpp` for the two-pass pattern the
server uses (one Scan-only indexing pass to record `{evio_event,
sub_event}` pairs; subsequent random access via the index).

### Legacy API (compatibility)

`DecodeEvent(i, event, ssp=nullptr, vtp=nullptr, tdc=nullptr)` still
exists as a thin compatibility wrapper that writes directly into
caller-owned structs without touching the lazy cache. New code should
prefer `SelectEvent()` + `Fadc()` / `Gem()` / `Tdc()` / `Vtp()`.

## Dependencies

- [evio](https://github.com/JeffersonLab/evio) (evio-6.0) — required
- [et](https://github.com/JeffersonLab/et) — optional, for EtChannel

Both are resolved from the Hall-B CODA installation by default; if not
found, CMake fetches from GitHub. Override with `-DEVIO_SOURCE=fetch` /
`-DET_SOURCE=fetch`.
