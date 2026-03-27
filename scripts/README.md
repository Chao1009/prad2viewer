# Scripts

Utility scripts for PRad-II detector visualization and analysis. Requires Python 3 with matplotlib and numpy.

```bash
pip install -r scripts/requirements.txt
```

## gem_layout.py

Visualize GEM strip layout from `gem_map.json`. Shows X/Y strips and APV boundaries in a 2x2 detector grid.

```bash
python scripts/gem_layout.py [path/to/gem_map.json]
```

Defaults to `database/gem_map.json`. Saves `gem_layout.png`.
