# Test & Diagnostic Tools

## evio_dump

EVIO file structure diagnostic tool.

```bash
evio_dump <file> [options]
```

| Option | Description |
|--------|-------------|
| `-m <mode>` | Analysis mode (default: summary) |
| `-n <N>` | Event count (tree mode) or event number (event mode) |
| `-D <file>` | Load DAQ configuration (for PRad etc.) |

Modes:

| Mode | Description |
|------|-------------|
| (default) | Summary: count events by tag |
| `tree` | Print bank tree for first N events (default N=5) |
| `tags` | List all unique bank tags with stats |
| `epics` | Dump all EPICS event text |
| `event` | Detailed dump of record N (1-based) |
| `triggers` | List trigger info for all events |

Examples:
```bash
evio_dump data.evio                                # event tag summary
evio_dump data.evio -m tree -n 10                  # bank tree for 10 events
evio_dump data.evio -m tags                        # all unique tags
evio_dump data.evio -m triggers -D prad_daq.json   # trigger bits (PRad config)
evio_dump data.evio -m event -n 100                # detailed dump of event 100
```

## evc_test

Basic smoke-test for the evc library. Reads EVIO buffers, decodes events, or connects to ET.

```bash
evc_test <file> [options]
```

| Option | Description |
|--------|-------------|
| `-m <mode>` | Mode: `scan` or `et` (default: read buffers) |
| `-s <N>` | Start event for scan mode (default: 1) |
| `-n <N>` | Number of events for scan mode (default: 50) |
| `-H <host>` | ET host (default: localhost) |
| `-P <port>` | ET port (default: 11111) |
| `-f <file>` | ET system file |
| `-S <station>` | ET station name |

Examples:
```bash
evc_test data.evio                                 # read and count buffers
evc_test data.evio -m scan -s 100 -n 20            # scan events 100-119
evc_test -m et -H localhost -P 11111 -f /tmp/et -S mon  # read from ET
```

## et_feeder

Replays an EVIO file into an ET system event-by-event at a configurable rate.

```bash
et_feeder <evio_file> [options]
```

| Option | Description |
|--------|-------------|
| `-h` | ET host (default: localhost) |
| `-p` | ET port (default: 11111) |
| `-f` | ET system file (default: /tmp/et_feeder) |
| `-i` | Interval between events in ms (default: 100) |
| `-s` | Start event number, 1-based (default: 1) |
| `-n` | Number of events to feed (default: all) |

Examples:
```bash
et_feeder data.evio -f /tmp/et_sys_prad -i 50       # feed all at 20 Hz
et_feeder data.evio -f /tmp/et_sys -s 1000 -n 500   # feed events 1000-1499
```

## evchan_test

Reads an EVIO file from both an ET station and disk simultaneously, printing word-by-word comparison. Used to verify ET transport fidelity.

```bash
evchan_test <evio_file> [-h host] [-p port] [-f et_file] [-i interval_ms]
```

Options are the same as `et_feeder`.

## ped_calc

Compute per-channel pedestals from EVIO data by selecting events with a specific trigger bit.

```bash
ped_calc <evio_file> -D <daq_config.json> [options]
```

| Option | Description |
|--------|-------------|
| `-D <file>` | DAQ configuration (required for PRad) |
| `-t <bit>` | Trigger bit to select (default: 3 = LMS_Alpha for PRad) |
| `-o <file>` | Output JSON file (default: pedestals_out.json) |
| `-n <N>` | Max events to process (default: all) |

Trigger bits (PRad):
- 0 = PHYS_LeadGlassSum (0x01)
- 1 = PHYS_TotalSum (0x02)
- 2 = LMS_Led (0x04)
- 3 = LMS_Alpha / pedestal (0x08)

Example:
```bash
ped_calc prad.evio -D prad_daq_config.json -t 3 -o pedestals.json
```

## gem_dump

GEM data diagnostic tool. Decodes SSP/MPD banks, runs GEM reconstruction, prints diagnostics or dumps events to JSON for visualization.

```bash
gem_dump <evio_file> [options]
```

| Option | Description |
|--------|-------------|
| `-m <mode>` | `summary` (default), `raw`, `hits`, `clusters`, `ped`, `evdump` |
| `-D <file>` | DAQ config (auto-searches `database/daq_config.json`) |
| `-G <file>` | GEM map (default: `gem_map.json` next to DAQ config) |
| `-P <file>` | Pedestal file (required for `hits`/`clusters`/`evdump`) |
| `-n <N>` | Max events (default: 10, 0=all) |
| `-e <N>` | Single event N (1-based) |
| `-t <bit>` | Trigger bit filter (-1=all) |
| `-o <file>` | Output file (`ped`/`evdump` modes) |
| `-z <sigma>` | Override zero-suppression threshold |
| `-f <filter>` | APV filter for `evdump` (see below) |

The `evdump` mode outputs JSON with three pipeline layers: `raw_apvs` (raw ADC), `zs_apvs` (zero-suppressed channels by APV address), and per-detector `clusters` + `hits_2d`. Use with `scripts/gem_cluster_view.py`.

By default evdump selects the first event with 2D hits. Use `-n K` for first K matching events, `-e N` for a specific event, or `-f` for custom APV-based filtering:

```
-f field=val[,val]:min_dets
```
Fields: `pos`, `plane` (X/Y), `match` (+Y/-Y), `orient`, `det`.

Examples:
```bash
gem_dump data.evio -m ped -o gem_ped.json                       # pedestals
gem_dump data.evio -P gem_ped.json -m summary -n 50             # overview
gem_dump data.evio -P gem_ped.json -m evdump -e 42              # dump event 42
gem_dump data.evio -P gem_ped.json -m evdump -f pos=10,11:3 -n 5  # beam-hole APVs in >=3 GEMs
python scripts/gem_cluster_view.py gem_event.json database/gem_map.json
```
