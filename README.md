# PRad-II Event Viewer & Monitor

EVIO decoder, event viewer, and online monitor for the PRad-II experiment at Jefferson Lab.
Decodes FADC250 waveforms (HyCal) and SSP/MPD data (GEM), with a web-based GUI for
waveform inspection, clustering, gain monitoring, and physics replay.
Also supports original PRad (ADC1881M) via DAQ configuration.

## Building

### Linux

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
```

CMake >= 3.14, C++17. Dependencies (`evio`, `et`) are resolved from the Hall-B CODA installation by default; if not found, CMake falls back to fetching from GitHub. Other dependencies (`nlohmann/json`, `websocketpp`, `asio`) are always fetched automatically.

Optional flags:
- `-DBUILD_ANALYSIS=ON` — ROOT-based replay and analysis tools (requires ROOT 6.0+)
- `-DBUILD_GUI=ON` — Qt standalone viewer and remote client (requires Qt6 or Qt5 WebEngine)
- `-DBUILD_PYTHON=ON` — pybind11 bindings `prad2py` (enables the Python GEM tools)
- `-DPython_EXECUTABLE=<path>` — force a specific Python interpreter (needed on systems with multiple `python3` versions in PATH — CMake otherwise may pick the wrong one)
- `-DWITH_ET=OFF` — disable ET (live monitoring); required on Windows
- `-DEVIO_SOURCE=fetch` — force fetching evio from GitHub (no CODA needed)

### Windows (MSYS2)

```powershell
# MSYS2 packages
pacman -S mingw-w64-x86_64-toolchain mingw-w64-x86_64-cmake mingw-w64-x86_64-ninja mingw-w64-x86_64-expat

# Build (PowerShell with MinGW in PATH)
$env:PATH = "C:\msys64\mingw64\bin;" + $env:PATH
cmake -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DWITH_ET=OFF
cmake --build build
```

Or build natively in WSL2.

## Server

`prad2_server` is a unified HTTP server that supports both file-based viewing and online ET monitoring. It replaces the former `prad2_viewer` and `prad2_monitor` executables.

```bash
prad2_server [evio_file] [-p port] [-H] [-i] [-f filter.json] [-c config.json] [-d data_dir] [-D daq_config.json] [--et]
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

## Tools

- [test/README.md](test/README.md) — generic EVIO diagnostics (`evio_dump`, `ped_calc`, plus dev tools in `test/dev/`).
- [gem/README.md](gem/README.md) — GEM tracker tools (`gem_dump`, `gem_event_viewer`, `gem_cluster_view`, `gem_layout`, …) plus detector reference notes.
- [analysis/README.md](analysis/README.md) — ROOT-based replay and physics analysis (binaries are installed with a `prad2ana_` prefix).
- [calibration/README.md](calibration/README.md) — HyCal gain-equalizer and operator calibration scan procedures.
- [scripts/README.md](scripts/README.md) — HyCal scaler / pedestal / gain / map-builder GUIs, waveform viewer, tagger viewer, plus `daq_tool/` (dev-only DAQ editors) and `shell/` (operator scripts).

## Installation

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/opt/prad2
cmake --build build -j$(nproc)
cmake --install build
```

Set up the environment (sets `PATH`, `LD_LIBRARY_PATH`, `PYTHONPATH`, `PRAD2_DATABASE_DIR`, `PRAD2_RESOURCE_DIR`):

```bash
source /opt/prad2/bin/prad2_setup.sh    # bash/zsh
source /opt/prad2/bin/prad2_setup.csh   # csh/tcsh
```

`prad2_setup.csh` bakes the install prefix in at configure time, so sourcing it from inside another script (e.g. a JLab farm wrapper that also `module load`s ROOT) works correctly. If you later move the install tree, `setenv PRAD2_DIR <new-prefix>` before sourcing to override.

### Python tools (optional)

The GEM Python tools (`gem_event_viewer`, `gem_cluster_view`, …) need `matplotlib`, `numpy`, and `PyQt6`. On shared systems where you can't `pip install` into site-packages, create a venv next to the install and point the wrapper at it:

```bash
python3 -m venv /opt/prad2/venv
source /opt/prad2/venv/bin/activate     # or activate.csh for tcsh
pip install -r /path/to/prad2evviewer/scripts/requirements.txt
```

Then source the venv's `activate` *before* sourcing `prad2_setup.csh` in your personal wrapper script.

**Watch out on RHEL 9** — CMake's `find_package(Python)` may pick up `/usr/bin/python3.12` over the venv's 3.9. Pass `-DPython_EXECUTABLE=$(which python3)` at configure time to pin the interpreter explicitly, otherwise `prad2py.cpython-312-*.so` won't be importable from your 3.9 venv.

## Configuration

`database/config.json` -- main config for the server (waveform, clustering, LMS, EPICS, elog, online/ET, color ranges).

PRad support: use `-D database/prad1/prad_daq_config.json` with the server.

## Project Structure

```
prad2dec/           libprad2dec — EVIO/ET reader, FADC250/SSP/ADC1881M decoders
prad2det/           libprad2det — HyCal/GEM clustering and reconstruction
python/             pybind11 bindings (prad2py)
src/                Server, Qt GUI, data source layer
resources/          Web frontend (HTML/CSS/JS)
database/           DAQ config, channel maps, calibration constants
test/               Generic EVIO diagnostic CLI tools (test/dev/ = not installed)
gem/                GEM tracker: gem_dump binary + Python tools + README
analysis/           ROOT-based replay and physics analysis (optional)
calibration/        HyCal calibration scan tools
scripts/            HyCal / tagger Python utilities (installed)
scripts/daq_tool/   Dev-only DAQ config editors (not installed)
scripts/shell/      Operator shell scripts (not installed)
cmake/              Build helpers (PradHelpers, WebDeps, bin-wrapper template)
docs/               ROL references, API documentation
```

## Contributors
Chao Peng -- Argonne National Laboratory\
Weizhi Xiong, Yuan Li, Mingyu Li -- Shandong University
