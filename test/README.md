# Test & Diagnostic Tools

## evio_dump

EVIO file structure inspector.

```bash
evio_dump <file> [-m mode] [-n N] [-D daq_config.json]
```

Modes: (default) summary by tag, `tree`, `tags`, `epics`, `event` (detail for event N), `triggers`.

## ped_calc

Compute HyCal per-channel pedestals from trigger-selected events.

```bash
ped_calc <evio_file> -D <daq_config.json> [-t bit] [-o file.json] [-n N]
```

Default trigger bit 3 (LMS_Alpha/pedestal for PRad).

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

## ET tools (requires `-DWITH_ET=ON`)

**evc_test** — Smoke-test: read EVIO buffers, scan events, or connect to ET.

**et_feeder** — Replay EVIO into ET at configurable rate: `et_feeder data.evio -f /tmp/et -i 50`

**evchan_test** — Compare ET-streamed vs disk-read events word-by-word.
