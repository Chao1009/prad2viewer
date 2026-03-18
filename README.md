# PRad2 Decoder

FADC250 waveform decoder and event viewer for PRad-II at Jefferson Lab.

## Building

```bash
mkdir build && cd build
cmake ..
make -j$(nproc)
```

CMake ≥ 3.14, C++17. Dependencies (`evio`, `et`, `nlohmann/json`, `websocketpp`, `asio`) fetched automatically. For prebuilt CODA libraries:

```bash
cmake .. -DEVIO_SOURCE=prebuilt -DET_SOURCE=prebuilt
```

## Event Viewer

```bash
evc_viewer <evio_file> [-p port] [-H] [-c hist_config.json] [-d data_dir]
```

| Option | Description |
|--------|-------------|
| `-p` | Server port (default 5050) |
| `-H, --hist` | Build per-channel histograms and occupancy on startup |
| `-c, --hist-config` | Histogram config file (implies `-H`, default: `database/hist_config.json`) |
| `-d, --data-dir` | Enable in-browser file picker sandboxed to this directory |

Open `http://localhost:5050` in a browser.

**GUI layout:**

Left panel is a HyCal geometry view. Color metric is selectable (peak integral, height, time, pedestal, occupancy) with editable range. Right panel has two histogram plots on top (integral with time cut, peak position), and a peaks table + waveform plot on the bottom. All dividers are draggable.

When `--data-dir` is set, the 📂 Open button lets you switch files from the browser. A "Process histograms" checkbox controls whether the new file gets a full histogram pass. Loading happens in the background with a progress bar.

## Online Monitor

```bash
evc_monitor [-p port] [-c online_config.json]
```

Connects to a running ET system. Same GUI as the file viewer, plus:

- ET connection status indicator
- Ring buffer of recent events (default 20) with dropdown selector
- **Clear Hist** button resets histograms and occupancy
- Auto-follows latest event; press **F** to resume after browsing

Events arrive via a background reader thread. The browser gets WebSocket push notifications, throttled to avoid flooding.

## Test Tools

```bash
evc_test <evio_file> [-v] [-t]          # decode and count channels
et_feeder <evio_file> [-h host] [-p port] [-f et_file] [-i interval_ms]
evchan_test <evio_file> [-h host] [-p port] [-f et_file] [-i interval_ms]
```

`et_feeder` replays an evio file into ET event-by-event at a configurable rate. `evchan_test` reads the same file from both an ET station and disk, printing a word-by-word comparison.

## Configuration

`database/hist_config.json`:
```json
{
    "hist": {
        "time_min": 170, "time_max": 190,
        "bin_min": 0, "bin_max": 20000, "bin_step": 100,
        "threshold": 3.0,
        "pos_min": 0, "pos_max": 400, "pos_step": 4
    }
}
```

`database/online_config.json`:
```json
{
    "et": { "host": "localhost", "port": 11111,
            "et_file": "/tmp/et_sys_prad2", "station": "prad2_monitor" },
    "ring_buffer_size": 20,
    "hist": { "...same as above..." }
}
```

| Field | Description |
|-------|-------------|
| `time_min/max` | Time window (ns) for integral histogram and time-cut occupancy |
| `bin_min/max/step` | Integral histogram binning |
| `pos_min/max/step` | Peak position histogram binning (ns) |
| `threshold` | Min peak height (ADC above pedestal) for histograms/occupancy |

## Project Structure

```
CMakeLists.txt
database/
    daq_map.json  hycal_modules.json  hist_config.json  online_config.json
prad2dec/                     → libprad2dec.a (see prad2dec/README.md)
resources/
    viewer.html  viewer.css  viewer.js
src/
    evc_viewer.cpp            File viewer (HTTP + evio decoder)
    evc_monitor.cpp           Online monitor (ET + WebSocket)
test/
    test_main.cpp             CLI decode/count tool
    et_feeder.cpp             Replay evio → ET
    evchan_test.cpp           EvChannel vs EtChannel comparison
```
