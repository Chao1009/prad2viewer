//=============================================================================
// replay_rawdata — convert EVIO file to ROOT tree
//
// Usage: replay_rawdata <input.evio> -o output_dir [-n max_events] [-p]
//   -o  output directory (REQUIRED)
//   -n  max events to process (default: all)
//   -p  include peak analysis branches
//=============================================================================

#include "Replay.h"
#include "InstallPaths.h"

#include <iostream>
#include <string>
#include <filesystem>
#include <cstdlib>
#include <getopt.h>

#ifndef DATABASE_DIR
#define DATABASE_DIR "."
#endif

static std::string makeOutputFile(const std::string &evio_path)
{
    std::string out = std::filesystem::path(evio_path).filename().string();
    auto pos = out.find(".evio");
    if (pos != std::string::npos)
        out = out.substr(0, pos) + out.substr(pos + 5);
    out += "_raw.root";
    return out;
}

int main(int argc, char *argv[])
{
    std::string input, output_dir, daq_config;
    int max_events = -1;
    bool peaks = false;

    std::string db_dir = prad2::resolve_data_dir(
        "PRAD2_DATABASE_DIR",
        {"../share/prad2evviewer/database"},
        DATABASE_DIR);
    daq_config = db_dir + "/daq_config.json"; // default DAQ config for PRad2

    int opt;
    while ((opt = getopt(argc, argv, "o:n:D:p")) != -1) {
        switch (opt) {
            case 'o': output_dir = optarg; break;
            case 'n': max_events = std::atoi(optarg); break;
            case 'D': daq_config = optarg; break;
            case 'p': peaks = true; break;
        }
    }
    if (optind < argc) input = argv[optind];

    if (input.empty() || output_dir.empty()) {
        std::cerr << "Usage: replay_rawdata <input.evio> -o output_dir [-D daq_config.json] [-n N] [-p]\n";
        std::cerr << "  -o  output directory (REQUIRED)\n";
        std::cerr << "  -D  DAQ config JSON (default: <db>/daq_config.json)\n";
        std::cerr << "  -n  max events to process (default: all)\n";
        std::cerr << "  -p  include peak analysis branches(required at most time)\n";
        return 1;
    }

    output_dir += "/" + makeOutputFile(input);

    analysis::Replay replay;
    if (!daq_config.empty()) replay.LoadDaqConfig(daq_config);
    replay.LoadDaqMap(db_dir + "/daq_map.json");
    std::cerr << "Using DAQ map: " << db_dir + "/daq_map.json" << "\n";

    if (!replay.Process(input, output_dir, max_events, peaks, daq_config))
        return 1;

    return 0;
}
