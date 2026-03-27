# PRad2 Event Viewer & Monitor

FADC250 waveform decoder, event viewer, and online monitor for PRad-II at Jefferson Lab. Also supports original PRad (ADC1881M) via DAQ configuration.

## Building

### Linux

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
```

CMake >= 3.14, C++17. Dependencies (`evio`, `et`, `nlohmann/json`, `websocketpp`, `asio`) fetched automatically.

Optional: `cmake -B build -DBUILD_ANALYSIS=ON` (ROOT 6.0+), `cmake -B build -DBUILD_GUI=ON` (Qt5).

For prebuilt CODA libraries: `cmake -B build -DEVIO_SOURCE=prebuilt -DET_SOURCE=prebuilt`

### Windows

File-based tools (`evc_viewer`, `gem_dump`, `evio_dump`, `ped_calc`) build on Windows with `-DWITH_ET=OFF`, which skips the ET library and live monitor.

**MSYS2 setup:**
```bash
pacman -S mingw-w64-x86_64-toolchain mingw-w64-x86_64-cmake mingw-w64-x86_64-ninja mingw-w64-x86_64-expat
```

**Build (PowerShell with MinGW in PATH):**
```powershell
$env:PATH = "C:\msys64\mingw64\bin;" + $env:PATH
cmake -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DWITH_ET=OFF
cmake --build build
```

**Alternative:** build natively in WSL2.

## Event Viewer

```bash
evc_viewer [evio_file] [-p port] [-H] [-c config.json] [-d data_dir] [-D daq_config.json]
```

Opens a web GUI at `http://localhost:5050` with tabs: Waveform Data, Clustering, Gain Monitoring (LMS), EPICS.

```bash
evc_viewer data.evio -H                                          # PRad-II
evc_viewer prad.evio -D database/prad1/prad_daq_config.json -H   # PRad
evc_viewer -d /data/stage6 -H                                    # file browser
```

## Online Monitor

```bash
evc_monitor [-p port] [-c config.json] [-D daq_config.json]
```

Same GUI as the viewer, connected to a live ET system. Includes event ring buffer, auto-follow, histogram accumulation, and elog report generation.

Test with `et_feeder`:
```bash
et_start -f /tmp/test_et -s 100000 -n 500
./bin/evc_monitor -D ../database/prad1/prad_daq_config.json -c ../database/prad1/prad_config.json
./bin/et_feeder prad.evio -f /tmp/test_et -i 50 -n 5000
```

## Test & Analysis Tools

See [test/README.md](test/README.md) (`evio_dump`, `ped_calc`, `gem_dump`, `evc_test`, `et_feeder`) and [analysis/README.md](analysis/README.md) (ROOT replay tools).

GEM visualization scripts: [scripts/README.md](scripts/README.md).

## Configuration

`database/config.json` — main config for viewer/monitor (waveform, clustering, LMS, EPICS, elog, color ranges).

PRad support: use `-D database/prad1/prad_daq_config.json` with viewer/monitor.

## Project Structure

```
database/           Config, DAQ maps, calibration, gem_map.json
prad2dec/           libprad2dec (EVIO/ET decoder, FADC250, SSP, waveform analysis)
prad2ana/           libprad2ana (HyCal/GEM clustering, reconstruction)
resources/          Web GUI (HTML/CSS/JS), report generation
src/                evc_viewer, evc_monitor, app_state
analysis/           Replay + physics analysis (optional, requires ROOT)
test/               Diagnostic tools (evio_dump, gem_dump, ped_calc, etc.)
scripts/            Python visualization (gem_layout, gem_cluster_view)
prad2qtmon/         Qt5 WebEngine client (optional)
```

## Contributors
Chao Peng — Argonne National Laboratory\
Yuan Li — Shandong University\
Mingyu Li — Shandong University
