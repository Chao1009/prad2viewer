# PRad2 Event Viewer & Monitor

FADC250 waveform decoder, event viewer, and online monitor for PRad-II at Jefferson Lab. Also supports original PRad (ADC1881M) via DAQ configuration.

## Building

### Linux

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
```

CMake >= 3.14, C++17. By default, `evio` and `et` are resolved from the Hall-B CODA installation (`EVIO_SOURCE=hallb`, `ET_SOURCE=hallb`). If the Hall-B libraries are not found, CMake automatically falls back to fetching them from GitHub. Other dependencies (`nlohmann/json`, `websocketpp`, `asio`) are always fetched automatically.

Optional flags:
- `-DBUILD_ANALYSIS=ON` — ROOT-based replay tools (requires ROOT 6.0+)
- `-DBUILD_GUI=ON` — Qt standalone viewer and remote client (requires Qt6 or Qt5 WebEngine)

To force fetching from source: `cmake -B build -DEVIO_SOURCE=fetch -DET_SOURCE=fetch`

### Windows

File-based tools (`prad2_server`, `gem_dump`, `evio_dump`, `ped_calc`) build on Windows with `-DWITH_ET=OFF`, which skips the ET library and live monitor.

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

## Server

`prad2_server` is a unified HTTP server that supports both file-based viewing and online ET monitoring. It replaces the former `prad2_viewer` and `prad2_monitor` executables.

```bash
prad2_server [evio_file] [-p port] [-H] [-i] [-c config.json] [-d data_dir] [-D daq_config.json] [--et]
```

Opens a web GUI at `http://localhost:5051` with tabs: Waveform Data, Clustering, Gain Monitoring (LMS), EPICS, Physics, GEM. See [docs/API.md](docs/API.md) for the full HTTP/WebSocket API reference and interactive CLI commands (`-i`).

### File mode (default)

```bash
prad2_server data.evio -H                                          # PRad-II
prad2_server prad.evio -D database/prad1/prad_daq_config.json -H   # PRad
prad2_server -d /data/stage6 -H                                    # file browser
```

### Online mode

```bash
prad2_server --et                       # connect to ET system (config from database/config.json)
prad2_server data.evio --et -H          # start with file, "Go Online" button available in UI
```

ET connection settings are read from the `"online"` section of `config.json`.

### Mode switching

The server maintains separate data accumulators for file and online modes -- switching modes does not discard data from the other mode. Use "Clear All" in the UI to reset the active mode's data.

- Loading a file switches to file mode
- The "Go Online" button in the web UI (or `--et` flag) switches to online mode
- The toggle button appears when both capabilities are available (compiled with `WITH_ET`)

Test with `et_feeder`:
```bash
et_start -f /tmp/test_et -s 100000 -n 500
./bin/prad2_server --et -D ../database/prad1/prad_daq_config.json -c ../database/prad1/prad_config.json
./bin/et_feeder prad.evio -f /tmp/test_et -i 50 -n 5000
```

### Remote access (JLab ifarm)

To run the server on a JLab ifarm node and view in your local browser:

**1. Start the server on ifarm** (pick a port, e.g. 5678):
```bash
ssh username@ifarm2402
cd /path/to/prad2evviewer/build
./bin/prad2_server data.evio -H -p 5678
```

**2. SSH tunnel from your local machine:**
```bash
ssh -L 5678:ifarm2402:5678 -J username@scilogin.jlab.org username@ifarm2402
```

**3. Open in your local browser:**
```
http://localhost:5678
```

This also works for online monitoring -- start `prad2_server --et -p 5678` on the ifarm node where the ET system is accessible.

## Qt Standalone Viewer (optional)

Build with `-DBUILD_GUI=ON` (uses Qt6 WebEngine by default, falls back to Qt5 if Qt6 is not found).

`prad2evviewer` embeds the server and web frontend in a native Qt window. No separate server process needed -- updates to the web frontend or server API are automatically reflected.

```bash
prad2evviewer                           # empty, use File > Open
prad2evviewer data.evio -H              # open file with histograms
prad2evviewer -d /data/stage6           # enable file browser
```

Features: native file dialogs (File > Open), drag-and-drop `.evio` files, View > Go Online (if compiled with ET), status bar with loading progress.

`prad2_client` is a thin Qt WebEngine wrapper that connects to a remote `prad2_server` instance:

```bash
prad2_client                            # connect to localhost:5051
prad2_client -H clonpc19 -p 8080       # connect to remote server
```

## Test & Analysis Tools

See [test/README.md](test/README.md) (`evio_dump`, `ped_calc`, `gem_dump`, `evc_test`, `et_feeder`) and [analysis/README.md](analysis/README.md) (ROOT replay tools).

GEM visualization scripts: [scripts/README.md](scripts/README.md).

## Installation

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/opt/prad2
cmake --build build -j$(nproc)
cmake --install build
```

Set up the environment (sets `PATH`, `LD_LIBRARY_PATH`, `PRAD2_DATABASE_DIR`, `PRAD2_RESOURCE_DIR`):

```bash
source /opt/prad2/bin/setup.sh    # bash/zsh
source /opt/prad2/bin/setup.csh   # csh/tcsh
```

## Configuration

`database/config.json` -- main config for the server (waveform, clustering, LMS, EPICS, elog, online/ET, color ranges).

PRad support: use `-D database/prad1/prad_daq_config.json` with the server.

## Project Structure

```
database/           Config, DAQ maps, calibration, gem_map.json
prad2dec/           libprad2dec (EVIO/ET decoder, FADC250, SSP, waveform analysis)
prad2det/           libprad2det (HyCal/GEM clustering, reconstruction)
resources/          Web GUI (HTML/CSS/JS), report generation
src/                prad2_server, prad2evviewer, prad2_client, viewer_server, app_state
analysis/           Replay + physics analysis (optional, requires ROOT)
test/               Diagnostic tools (evio_dump, gem_dump, ped_calc, etc.)
scripts/            Python visualization (gem_layout, gem_cluster_view)
```

## Contributors
Chao Peng -- Argonne National Laboratory\
Yuan Li -- Shandong University\
Mingyu Li -- Shandong University
