# prad2dec

Static library for reading CODA EVIO data and decoding FADC250/SSP waveforms.

## Components

- **EvChannel** — Reads EVIO files, scans bank trees, decodes FADC250 composite data and SSP/MPD banks
- **EtChannel** — Subclass that reads from a live ET system (requires `-DWITH_ET=ON`)
- **Fadc250Decoder** — Decodes composite bank payload into `SlotData`/`ChannelData` structs
- **SspDecoder** — Decodes SSP fiber data into `MpdData`/`ApvData` (GEM readout)
- **WaveAnalyzer** — Waveform analysis: pedestal, peak search, integration

## Usage

```cpp
evc::EvChannel ch;
ch.Open("data.evio");
fdec::EventData event;

while (ch.Read() == evc::status::success) {
    if (!ch.Scan()) continue;
    for (int i = 0; i < ch.GetNEvents(); ++i) {
        ch.DecodeEvent(i, event);
        // access event.rocs[r].slots[s].channels[c].samples[]
    }
}
```

## Dependencies

- [evio](https://github.com/JeffersonLab/evio) (evio-6.0) — always required
- [et](https://github.com/JeffersonLab/et) — optional, for EtChannel (`-DWITH_ET=ON`)

Both fetched by CMake or linked from prebuilt CODA installation.
