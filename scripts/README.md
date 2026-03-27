# Scripts

Python utilities for GEM visualization. Requires matplotlib and numpy.

## gem_layout.py

Visualize GEM strip layout from `gem_map.json`. Shows X/Y strips, APV boundaries, and beam hole.

```bash
python scripts/gem_layout.py [gem_map.json]
```

## gem_cluster_view.py

Visualize GEM clustering from `gem_dump -m evdump` JSON. Strip geometry is derived from APV addresses via `gem_strip_map.py`, so beam-hole half-strips are drawn correctly.

```bash
python scripts/gem_cluster_view.py <event.json> [gem_map.json] [--det N] [-o file.png]
```

Shows: fired X strips (blue) and Y strips (red) color-coded by charge, cluster center markers, 2D hit positions, and beam hole. Prints a cluster summary table to the terminal.

```bash
gem_dump data.evio -m ped -o gem_ped.json
gem_dump data.evio -P gem_ped.json -m evdump -e 42
python scripts/gem_cluster_view.py gem_event.json database/gem_map.json
```
