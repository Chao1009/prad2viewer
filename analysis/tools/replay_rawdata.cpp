//=============================================================================
// replay_rawdata — convert EVIO file to ROOT tree
//
// Usage: replay_rawdata <input.evio> [-o output.root] [-n max_events] [-p]
//   -o  output ROOT file (default: input with _raw.root extension)
//   -n  max events to process (default: all)
//   -p  include peak analysis branches
//=============================================================================

#include "Replay.h"
#include <iostream>
#include <string>
#include <cstdlib>
#include <getopt.h>

#ifndef DATABASE_DIR
#define DATABASE_DIR "."
#endif

int main(int argc, char *argv[])
{
    std::string input, output, daq_config;
    int max_events = -1;
    bool peaks = false;

    std::string db_dir = DATABASE_DIR;
    if (const char *env = std::getenv("PRAD2_DATABASE_DIR"))  db_dir = env;
    daq_config = db_dir + "/daq_config.json"; // default DAQ config for PRad2

    int opt;
    while ((opt = getopt(argc, argv, "o:n:D:p")) != -1) {
        switch (opt) {
            case 'o': output = optarg; break;
            case 'n': max_events = std::atoi(optarg); break;
            case 'D': daq_config = optarg; break;
            case 'p': peaks = true; break;
        }
    }
    if (optind < argc) input = argv[optind];

    if (input.empty()) {
        std::cerr << "Usage: replay_rawdata <input.evio> [-o output.root] [-D daq_config.json] [-n N] [-p]\n";
        return 1;
    }

    if (output.empty()) {
        output = input;
        auto pos = output.find(".evio");
        if (pos != std::string::npos) output = 
            output.substr(0, pos) + output.substr(pos + 5);
        output += "_raw.root";
    }

    analysis::Replay replay;
    if (!daq_config.empty()) replay.LoadDaqConfig(daq_config);
    replay.LoadDaqMap(db_dir + "/daq_map.json");
    std::cerr << "Using DAQ map: " << db_dir + "/daq_map.json" << "\n";

    if (!replay.Process(input, output, max_events, peaks, daq_config))
        return 1;

    return 0;
}
