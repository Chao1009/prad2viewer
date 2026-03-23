# PRad2 Event Viewer & Monitor

FADC250 waveform decoder, event viewer, and online monitor for PRad-II at Jefferson Lab. Also supports original PRad (ADC1881M) via DAQ configuration.

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

For prebuilt CODA libraries: `cmake .. -DEVIO_SOURCE=prebuilt -DET_SOURCE=prebuilt`

For offline analysis tools (requires ROOT 6.0+): `cmake .. -DBUILD_ANALYSIS=ON`

For the Qt5 monitor client: `cmake .. -DBUILD_GUI=ON`

## Event Viewer

```bash
evc_viewer [evio_file] [-p port] [-H] [-c config.json] [-d data_dir] [-D daq_config.json]
```

| Option | Description |
|--------|-------------|
| `-p` | Server port (default 5050) |
| `-H` | Build histograms on startup |
| `-c` | Main config file override |
| `-d` | Enable in-browser file picker sandboxed to this directory |
| `-D` | DAQ config (for PRad: `database/prad1/prad_daq_config.json`) |

Examples:
```bash
evc_viewer data.evio -H                                                    # PRad-II
evc_viewer prad.evio -D database/prad1/prad_daq_config.json -H             # PRad
evc_viewer -d /data/stage6 -H                                              # file browser
```

Open `http://localhost:5050`. Three tabs:

- **Waveform Data** — HyCal geometry colored by peak metrics, waveform display, per-channel histograms
- **Clustering** — Island clustering with energy/ncluster/nblocks histograms, cluster table
- **Gain Monitoring** — Per-module LMS signal vs time, reference channel correction, drift warnings

Geo view: scroll-wheel zoom, left-drag pan, double-click or Reset to restore.

## Online Monitor

```bash
evc_monitor [-p port] [-c config.json] [-D daq_config.json]
```

Connects to a running ET system. Same GUI as the viewer, plus:

- ET connection status indicator
- Ring buffer of recent events with dropdown selector
- Per-tab **Clear** buttons reset tab-specific data
- Auto-follows latest event; press **F** to resume after browsing
- Gain monitoring and clustering accumulate continuously

### Testing with et_feeder

Replay an EVIO file into a local ET system to test the monitor:

```bash
# Terminal 1: start ET system
et_start -f /tmp/test_et -s 100000 -n 500

# Terminal 2: start monitor (PRad example)
./bin/evc_monitor -D ../database/prad1/prad_daq_config.json \
  -c ../database/prad1/prad_config.json

# Terminal 3: feed events
./bin/et_feeder prad.evio -f /tmp/test_et -i 50 -s 10000 -n 5000
```

The `-c` flag points to a config file with ET connection settings and experiment-specific parameters (trigger bits, histogram ranges, etc.).

## Monitor Client (Qt)

```bash
prad2qtmon                        # http://localhost:5051
prad2qtmon -H clonpc19 -p 8080   # http://clonpc19:8080
```

Requires `-DBUILD_GUI=ON`.

## Test Tools

See [test/README.md](test/README.md) for `evio_dump`, `evc_test`, `et_feeder`, `evchan_test`, and `ped_calc`.

## Configuration

### `database/config.json`

Main configuration file for both viewer and monitor:

| Section | Key fields |
|---------|------------|
| `online` | `et_host`, `et_port`, `et_file`, `et_station`, `ring_buffer_size`, `refresh_ms` |
| `waveform` | `time_cut`, `integral_hist`, `time_hist`, `thresholds` |
| `clustering` | `min_module_energy`, `min_cluster_energy`, `skip_trigger_bits`, `energy_hist`, `nclusters_hist`, `nblocks_hist` |
| `lms_monitor` | `trigger_bit`, `warn_threshold`, `warn_min_mean`, `max_history`, `reference_channels` |
| `color_ranges` | Per-tab:metric color range defaults (e.g. `"dq:integral": [0, 10000]`) |
| `calibration` | `adc_to_mev`, `calibration_file` |

### PRad Support

PRad-specific files in `database/prad1/`:

| File | Description |
|------|-------------|
| `prad_daq_config.json` | Event tags, ROC IDs, TI format, bank tags, pedestal/DAQ map refs |
| `adc1881m_pedestals.json` | Per-channel pedestals (1733 channels) |
| `prad_daq_map.json` | HyCal DAQ mapping for Fastbus ADC1881M |
| `prad_calibration.json` | Per-module calibration constants |

Use `-D database/prad1/prad_daq_config.json` with both viewer and monitor.

## Project Structure

```
database/
    config.json                 Main config (waveform, clustering, LMS, calibration, online)
    daq_config.json             PRad-II DAQ defaults
    daq_map.json                PRad-II HyCal DAQ mapping
    hycal_modules.json          HyCal module geometry
    prad1/                      PRad-specific config and calibration
prad2dec/                       libprad2dec.a (EVIO decoder library)
prad2ana/                       libprad2ana.a (HyCal clustering + analysis)
resources/
    viewer.html/css/js          Web GUI (shared by viewer and monitor)
src/
    evc_viewer.cpp              File viewer (HTTP server)
    evc_monitor.cpp             Online monitor (ET + WebSocket)
    app_state.h/cpp             Shared application state
    viewer_utils.h              Common types and helpers
analysis/                       Replay + physics analysis (requires ROOT, see analysis/README.md)
test/                           Diagnostic tools (see test/README.md)
prad2qtmon/                     Qt5 WebEngine client (optional)
```

## Contributors
Chao Peng — Argonne National Laboratory\
Yuan Li — Shandong University\
Mingyu Li — Shandong University

