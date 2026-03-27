# Offline Analysis Tools

Replay and physics analysis for PRad2. **Requires ROOT 6.0+.**

```bash
cmake -B build -DBUILD_ANALYSIS=ON
cmake --build build -j$(nproc)
```

## Tools

**replay_rawdata** — EVIO to ROOT tree with per-channel waveform data.
```bash
replay_rawdata <input.evio> [-o output.root] [-n max_events] [-p]
```

**replay_recon** — HyCal reconstruction replay with clustering and per-module energy histograms.
```bash
replay_recon <input.evio> [-o output.root] [-c config.json] [-D daq_config.json] [-n N]
```

## Adding a Tool

Create `tools/my_tool.cpp`, then add to `CMakeLists.txt`:
```cmake
add_analysis_tool(my_tool tools/my_tool.cpp)
```

Shared sources (`Replay.cpp`, `PhysicsTools.cpp`) and dependencies linked automatically.

## Contributors
Yuan Li — Shandong University
