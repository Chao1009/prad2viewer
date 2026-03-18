# PRad2 Decoder

FADC250 waveform decoder and event viewer for PRad-II at Jefferson Lab.

Reads CODA evio files, decodes composite data banks (`c,i,l,N(c,Ns)` format), performs waveform analysis (pedestal, peak finding, integration), and provides a web-based event display.

## Building

```bash
mkdir build && cd build
cmake ..
make -j$(nproc)
```

Requires CMake ≥ 3.14 and a C++17 compiler. Dependencies (`evio`, `et`, `nlohmann/json`, `websocketpp`, `asio`) are fetched automatically.

To use prebuilt CODA libraries instead of fetching:
```bash
cmake .. -DEVIO_SOURCE=prebuilt -DET_SOURCE=prebuilt
```

## Usage

### Command-line test
```bash
# Channel hit counts (JSON output)
./bin/evc_test data.evio

# Verbose: print all waveforms
./bin/evc_test data.evio -v

# Print bank tree
./bin/evc_test data.evio -t
```

### Event viewer
```bash
./bin/evc_viewer data.evio [port]
# Open http://localhost:5050
```

The viewer auto-discovers `database/hycal_modules.json` and `database/daq_map.json` at the compile-time `DATABASE_DIR` path.

#### Histogram mode
Use `--hist` to do a one-time pass over all events and accumulate peak integral histograms per channel. The largest peak within a configurable sample range is selected for each channel per event.

```bash
# Use default config from database/hist_config.json
./bin/evc_viewer data.evio --hist

# Use a custom config file
./bin/evc_viewer data.evio --hist my_hist_config.json
```

Example `hist_config.json`:
```json
{
    "hist": {
        "time_min": 170,
        "time_max": 190,
        "bin_min": 0,
        "bin_max": 20000,
        "bin_step": 100,
        "threshold": 3.0
    }
}
```

| Field | Description |
|-------|-------------|
| `time_min`, `time_max` | Peak time range in nanoseconds |
| `bin_min`, `bin_max`, `bin_step` | Histogram axis range and bin width |
| `threshold` | Minimum peak height (ADC above pedestal) to count |

Without `--hist`, the histogram panel in the GUI shows empty.

### Online monitor (ET)
```bash
./bin/evc_monitor [port] [--config online_config.json]
# Open http://localhost:5051
```

Connects to an ET system and monitors events in real time. The viewer auto-discovers `database/online_config.json` at the compile-time `DATABASE_DIR` path.

Features:
- Background thread reads events from ET continuously
- Histograms accumulate until the user clicks **Clear Hist** in the GUI
- Latest N events stored in a ring buffer (default 20, configurable)
- WebSocket push notifies the viewer on each new event
- Viewer auto-follows the latest event; select older events from the dropdown
- Arrow keys navigate the ring buffer; press **F** to resume auto-follow

Example `online_config.json`:
```json
{
    "et": {
        "host": "localhost",
        "port": 11111,
        "et_file": "/tmp/et_sys_prad2",
        "station": "prad2_monitor"
    },
    "ring_buffer_size": 20,
    "hist": { ... same fields as hist_config.json ... }
}
```

## Project Structure

```
CMakeLists.txt
database/
    daq_map.json              DAQ channel map (crate/slot/channel → module)
    hycal_modules.json        Module geometry (name, type, position, size)
    hist_config.json          Histogram settings for file viewer (optional)
    online_config.json        ET + histogram settings for online monitor
prad2dec/                     Static library: libprad2dec.a
    include/
        EvStruct.h            Evio header parsers and EvNode tree structure
        EvChannel.h           Evio file reader and bank tree scanner
        EtChannel.h           ET system live reader
        EtConfigWrapper.h     ET configuration C++ wrapper
        Fadc250Data.h         Pre-allocated event data (zero-alloc flat arrays)
        Fadc250Decoder.h      Composite bank decoder
        WaveAnalyzer.h        Waveform analysis (pedestal, peaks, integration)
    src/
        EvChannel.cpp
        EtChannel.cpp
        Fadc250Decoder.cpp
        WaveAnalyzer.cpp
resources/
    viewer.html               Web-based event display (file + online modes)
src/
    evc_viewer.cpp            File viewer HTTP server (websocketpp/asio)
    evc_monitor.cpp           Online monitor: ET reader + WebSocket push
test/
    test_main.cpp             CLI test and channel counting tool
```

## Data Format

Each evio event contains ROC banks (tags `0x80`–`0x8c`), each holding a composite bank (tag `0xe101`) with format string `c,i,l,N(c,Ns)`:

| Field | Type | Description |
|-------|------|-------------|
| `c` | uint8 | Slot number (3–20) |
| `i` | int32 | Event number |
| `l` | int64 | 48-bit timestamp |
| `N` | uint32 | Number of channels fired |
| `c` | uint8 | Channel number (0–15) |
| `N` | uint32 | Number of samples |
| `s` | int16[] | ADC sample values |

Slots repeat back-to-back within a single ROC's composite payload.

## Library API

```cpp
#include "EvChannel.h"
#include "Fadc250Data.h"
#include "WaveAnalyzer.h"

evc::EvChannel ch;
ch.Open("data.evio");

fdec::EventData event;          // allocate once, reuse
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
                    // wres.ped.mean, wres.peaks[0].height, ...
                }
            }
        }
    }
}
```
