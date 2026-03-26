//=============================================================================
// replay_rawdata — convert EVIO file to ROOT tree
//
// Usage: replay_rawdata <input.evio> [-o output.root] [-n max_events] [-p]
//   -o  output ROOT file (default: input with .root extension)
//   -n  max events to process (default: all)
//   -p  include peak analysis branches
//=============================================================================

#include "Replay.h"
#include <iostream>
#include <string>
#include <getopt.h>

int main(int argc, char *argv[])
{
    std::string input, output, daq_config;
    int max_events = -1;
    bool peaks = false;

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
        output += ".root";
    }

    analysis::Replay replay;
    if (!daq_config.empty()) replay.LoadDaqConfig(daq_config);

    if (!replay.Process(input, output, max_events, peaks))
        return 1;

    return 0;
}
