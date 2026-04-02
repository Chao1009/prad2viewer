# PRad-II Event Viewer — Server API Reference

Base URL: `http://localhost:<port>` (default port 5051).

All responses are JSON unless noted otherwise. The server also pushes real-time updates over a WebSocket connection on the same port.

---

## Server

| Method | Endpoint | Description |
|--------|----------|-------------|
| GET | `/api/config` | Server configuration and capabilities |
| GET | `/api/progress` | File-loading progress (`loading`, `phase`, `current`, `total`, `file`) |

## Mode Switching

| Method | Endpoint | Description |
|--------|----------|-------------|
| POST | `/api/mode/online` | Switch to ET/online mode. Optional JSON body: `{"host","port","et_file","station"}` |
| POST | `/api/mode/file` | Switch to file/offline mode |

## Events

| Method | Endpoint | Description |
|--------|----------|-------------|
| GET | `/api/event/<n>` | Decoded event `n` (1-based in file mode; seq number in online mode) |
| GET | `/api/event/latest` | Latest event from the ring buffer (online mode) |
| GET | `/api/waveform/<n>/<roc_slot_ch>` | On-demand waveform samples for a single channel (file mode) |
| GET | `/api/clusters/<n>` | Cluster reconstruction for event `n` |
| GET | `/api/ring` | Ring buffer summary: list of seq numbers and latest seq |

## File Browser

| Method | Endpoint | Description |
|--------|----------|-------------|
| GET | `/api/files` | List available files under `data_dir` |
| GET | `/api/load?file=<path>&hist=0\|1` | Load an evio file (relative to `data_dir`). `hist=1` enables histograms |

## Histograms & Occupancy

| Method | Endpoint | Description |
|--------|----------|-------------|
| GET | `/api/occupancy` | Per-module occupancy (hit counts and integrals) |
| GET | `/api/cluster_hist` | Cluster-level histograms |
| GET | `/api/hist/<module_key>` | Amplitude histogram for a module |
| GET | `/api/poshist/<module_key>` | Position histogram for a module |
| POST | `/api/hist/clear` | Clear all amplitude/position histograms |

## LMS (Laser Monitoring)

| Method | Endpoint | Description |
|--------|----------|-------------|
| GET | `/api/lms/summary[?ref=<ch>]` | LMS summary across all modules. Optional `ref` channel for normalization |
| GET | `/api/lms/<module>[?ref=<ch>]` | LMS time series for a single module |
| GET | `/api/lms/refs` | List available LMS reference channels |
| POST | `/api/lms/clear` | Clear LMS accumulator |

## EPICS

| Method | Endpoint | Description |
|--------|----------|-------------|
| GET | `/api/epics/channels` | List of known EPICS channel names |
| GET | `/api/epics/latest` | Latest values for all EPICS channels |
| GET | `/api/epics/channel/<name>` | Time series for a single EPICS channel |
| GET | `/api/epics/batch?ch=<n1>&ch=<n2>` | Batch fetch multiple channels |
| POST | `/api/epics/clear` | Clear EPICS history |

## GEM

| Method | Endpoint | Description |
|--------|----------|-------------|
| GET | `/api/gem/config` | GEM detector geometry and strip mapping |
| GET | `/api/gem/hits` | GEM hits for the current event |
| GET | `/api/gem/occupancy` | GEM strip occupancy histograms |
| GET | `/api/gem/hist` | GEM amplitude histograms |

## Physics

| Method | Endpoint | Description |
|--------|----------|-------------|
| GET | `/api/physics/energy_angle` | Energy vs. angle distribution |
| GET | `/api/physics/moller` | Moller scattering analysis |

## Elog

| Method | Endpoint | Description |
|--------|----------|-------------|
| POST | `/api/elog/post` | Post to the electronic logbook. JSON body: `{"xml": "<elog XML>"}` |

---

## WebSocket Messages

Connect to `ws://localhost:<port>` for real-time push notifications.

### Server → Client

| `type` | Fields | Description |
|--------|--------|-------------|
| `status` | `connected`, `waiting`, `retries` | ET connection status change |
| `new_event` | `seq` | New event available in ring buffer |
| `mode` | `mode` | Mode changed (`file`, `online`, `idle`) |
| `file_loaded` | (full config) | File finished loading |
| `load_progress` | `phase`, `current`, `total` | File load progress update |
| `hist_cleared` | — | Histograms were cleared |
| `lms_cleared` | — | LMS data was cleared |
| `lms_event` | `count` | New LMS trigger event |
| `epics_cleared` | — | EPICS data was cleared |
| `epics_event` | `count` | New EPICS event received |

---

## CLI Interactive Mode

Start the server with `-i` (or `--interactive`) to enable the stdin command interface:

```
prad2_server data.evio -H -i
```

| Command | Description |
|---------|-------------|
| `status` | Show current mode, file info, ET connection |
| `load <path> [1]` | Load an evio file (append `1` for histograms) |
| `online` | Switch to ET/online mode |
| `offline` | Switch to file/offline mode |
| `clear hist\|lms\|epics` | Clear accumulators |
| `quit` / `exit` | Stop the server |
| `help` | Show command list |
