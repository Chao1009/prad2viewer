# Offline Analysis Tools

Replay and physics analysis tools for PRad2. **Requires ROOT 6.0+**.

Not built by default. Enable with:

```bash
cmake .. -DBUILD_ANALYSIS=ON
make -j$(nproc)
```

ROOT must be discoverable by CMake. On JLab machines:
```bash
source /apps/root/6.28.06/setroot_CUE.bash
cmake .. -DBUILD_ANALYSIS=ON && make -j$(nproc)
```

## Tools

### replay_rawdata

Convert EVIO file to ROOT tree with per-channel waveform data.

```bash
replay_rawdata <input.evio> [-o output.root] [-n max_events] [-p]
```

| Option | Description |
|--------|-------------|
| `-o` | Output ROOT file (default: input with .root extension) |
| `-n` | Max events to process (default: all) |
| `-p` | Include peak analysis branches (height, time, integral) |

### replay_hycalRecon

HyCal reconstruction replay — runs clustering and writes reconstructed hits to ROOT tree with per-module energy histograms.

```bash
replay_hycalRecon <input.evio> [-o output.root] [-c config.json] [-D daq_config.json] [-n N]
```

## Shared Sources

Compiled into each tool (no separate library):

- **`Replay`** (`include/Replay.h`, `src/Replay.cpp`) — EVIO to ROOT tree conversion
- **`PhysicsTools`** (`include/PhysicsTools.h`, `src/PhysicsTools.cpp`) — per-module energy histograms, Gaussian peak/resolution fitting, elastic kinematics (e-p, e-e), energy loss corrections

## Structure

```
analysis/
    include/            Shared headers
    src/                Shared implementations
    tools/              Executable entry points
    CMakeLists.txt      Isolated ROOT dependency
```

## Adding a New Tool

1. Create `tools/my_tool.cpp` with `main()`
2. Add one line to `CMakeLists.txt`:
   ```cmake
   add_analysis_tool(my_tool tools/my_tool.cpp)
   ```

All shared sources (`Replay.cpp`, `PhysicsTools.cpp`) and dependencies (`prad2dec`, `prad2ana`, ROOT) are linked automatically.

## Contributors
Yuan Li — Shandong University
