# PRad2 Event Viewer

FADC250 waveform decoder and event viewer for PRad-II at Jefferson Lab. Also supports original PRad (ADC1881M) via DAQ configuration.

## Building

```bash
mkdir build && cd build
cmake ..
make -j$(nproc)
```

CMake >= 3.14, C++17. Dependencies (`evio`, `et`, `nlohmann/json`, `websocketpp`, `asio`) fetched automatically. To speed up rebuilds, set a persistent fetch cache:

```bash
export FETCHCONTENT_BASE_DIR=$HOME/.cmake/fetchcontent
```

For prebuilt CODA libraries:
```bash
cmake .. -DEVIO_SOURCE=prebuilt -DET_SOURCE=prebuilt
```

To build the Qt5 WebEngine monitor client:
```bash
cmake .. -DBUILD_GUI=ON
```

## Event Viewer

```bash
evc_viewer [evio_file] [-p port] [-H] [-c hist_config.json] [-d data_dir] [-D daq_config.json]
```

| Option | Description |
|--------|-------------|
| `-p` | Server port (default 5050) |
| `-H` | Build per-channel histograms and occupancy on startup |
| `-c` | Histogram config file (implies `-H`) |
| `-d` | Enable in-browser file picker sandboxed to this directory |
| `-D` | DAQ configuration file (for PRad: `database/prad1/prad_daq_config.json`) |

Examples:
```bash
evc_viewer data.evio -H                            # PRad-II with histograms
evc_viewer prad.evio -D database/prad1/prad_daq_config.json -H  # PRad
evc_viewer -d /data/stage6 -H                      # browse and pick from GUI
```

Open `http://localhost:5050`. Three tabs:

- **Data Quality** — HyCal geometry colored by peak metrics, waveform display, per-channel histograms
- **Clustering** — Island clustering with energy histogram, cluster table, per-event reconstruction
- **LMS Monitoring** — Per-module LMS signal vs time, reference channel normalization, drift warnings

Geo view supports scroll-wheel zoom, left-drag pan, double-click or Reset button to restore.

## Online Monitor

```bash
evc_monitor [-p port] [-c online_config.json] [-H hist_config.json] [-D daq_config.json]
```

Connects to a running ET system. Same GUI as the viewer, plus:

- ET connection status indicator
- Ring buffer of recent events with dropdown selector
- **Clear Hist** button resets histograms and occupancy
- Auto-follows latest event; press **F** to resume after browsing
- LMS monitoring runs continuously regardless of client connections

## Monitor Client (Qt)

A lightweight Qt5 WebEngine wrapper.

```bash
prad2qtmon                        # http://localhost:5051
prad2qtmon -H clonpc19 -p 8080   # http://clonpc19:8080
```

Requires `-DBUILD_GUI=ON`.

## Test Tools

See [test/README.md](test/README.md) for detailed usage of `evio_dump`, `evc_test`, `et_feeder`, `evchan_test`, and `ped_calc`.

## Configuration

### `database/hist_config.json`

Histogram and waveform analysis settings (shared by viewer and monitor):

| Field | Description |
|-------|-------------|
| `time_min/max` | Time window (ns) for integral histogram and time-cut occupancy |
| `bin_min/max/step` | Integral histogram binning |
| `pos_min/max/step` | Peak position histogram binning (ns) |
| `threshold` | Min peak height (ADC above pedestal) |
| `min_peak_ratio` | Secondary peak must rise >= this fraction of main peak |

### `database/reconstruction.json`

Clustering, LMS monitoring, and calibration:

| Section | Key fields |
|---------|------------|
| `clustering` | `min_module_energy`, `min_cluster_energy`, `skip_trigger_bits`, `energy_hist` |
| `lms_monitor` | `trigger_bit`, `warn_threshold`, `max_history`, `reference_channels` |
| `calibration` | `adc_to_mev`, `calibration_file` |

### `database/online_config.json`

ET connection settings (monitor only):

```json
{
    "et": { "host": "localhost", "port": 11111,
            "et_file": "/tmp/et_sys_prad2", "station": "prad2qtmon" },
    "ring_buffer_size": 20
}
```

### PRad Support

PRad-specific data files are in `database/prad1/`:

| File | Description |
|------|-------------|
| `prad_daq_config.json` | DAQ config: event tags, ROC IDs, TI format, bank tags, pedestal/DAQ map refs |
| `adc1881m_pedestals.json` | Per-channel pedestals (1733 channels) |
| `prad_daq_map.json` | HyCal DAQ mapping for Fastbus ADC1881M |
| `prad_calibration.json` | Per-module calibration constants |

Use `-D database/prad1/prad_daq_config.json` with both viewer and monitor.

## Project Structure

```
CMakeLists.txt
database/
    hist_config.json  reconstruction.json  online_config.json
    daq_map.json  hycal_modules.json  daq_config.json
    prad1/                  PRad-specific config and calibration
prad2dec/                   libprad2dec.a (EVIO decoder library)
prad2ana/                   libprad2ana.a (HyCal clustering + analysis)
resources/
    viewer.html  viewer.css  viewer.js
src/
    evc_viewer.cpp          File viewer (HTTP server)
    evc_monitor.cpp         Online monitor (ET + WebSocket)
    app_state.h/cpp         Shared state between viewer and monitor
    viewer_utils.h          Common types and helpers
test/                       Diagnostic tools (see test/README.md)
prad2qtmon/              Qt5 WebEngine client (optional)
```
