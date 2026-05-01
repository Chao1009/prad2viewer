# prad2det

Static library for PRad-II detector reconstruction — takes decoded data
from `prad2dec` and produces physics-level outputs (clusters, 2-D hits).

## Components

- **HyCalSystem** — HyCal geometry, DAQ map, calibration constants, sector-grid neighbor lookup
- **HyCalCluster** — Island clustering + log-weighted center-of-gravity
  on top of a `HyCalSystem`
  (see [`docs/technical_notes/hycal_clustering/hycal_clustering.md`](../docs/technical_notes/hycal_clustering/hycal_clustering.md)
  for a worked example with plots)
- **GemSystem** — GEM detector hierarchy, pedestal / common-mode / zero-suppression, strip mapping (full-readout vs online-ZS auto-detected per APV)
- **GemCluster** — Per-plane strip clustering + X/Y matching (Cartesian
  with cuts or ADC-sorted)
  (see [`docs/technical_notes/gem_clustering/gem_clustering.md`](../docs/technical_notes/gem_clustering/gem_clustering.md)
  for a worked example with plots)
- **DetectorTransform** — 2×3 rotation + translation matrix for detector → lab-frame geometry
- **EpicsStore** — EPICS slow-control snapshot accumulator with look-up by event number
- **EventData** — Shared HyCal module-hit / cluster-hit data structures

## Usage

```cpp
// GEM — per-event pipeline
gem::GemSystem gsys;
gsys.Init("database/gem_daq_map.json");
gsys.LoadPedestals("gem_ped.json");     // required for full-readout data
gem::GemCluster gcl;

for each event {
    gsys.Clear();
    gsys.ProcessEvent(ssp_event_data);  // from prad2dec
    gsys.Reconstruct(gcl);
    for (auto &h : gsys.GetAllHits()) { /* h.x, h.y, h.det_id, ... */ }
}

// HyCal — feed per-module energies, cluster, read back
fdec::HyCalSystem hsys;
hsys.Init("database/hycal_modules.json", "database/hycal_daq_map.json");
hsys.LoadCalibration("database/hycal_calib.json");
fdec::HyCalCluster hcl(hsys);

hcl.Clear();
for (auto &[idx, E] : hits) hcl.AddHit(idx, E);
hcl.FormClusters();
std::vector<fdec::ClusterHit> out;
hcl.ReconstructHits(out);
```

## Dependencies

- [prad2dec](../prad2dec/README.md) — for the SSP/FADC event-data types
- [nlohmann/json](https://github.com/nlohmann/json) — fetched automatically for config loading
