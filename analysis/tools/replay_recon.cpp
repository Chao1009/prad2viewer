//=============================================================================
// replay_recon — decode EVIO file and do detector reconstruction write to ROOT tree
//
// Usage: replay_recon <input.evio> [-o output.root] [-D daq_config.json] [-p]
//                                  [-g gem_pedestal.json] [-z zerosup_threshold]
//   -o  output ROOT file (default: input with _recon.root extension)
//   -D  DAQ config file (default: PRAD2_DATABASE_DIR/daq_config.json)
//   -p  read prad1 data and do not include GEM
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
    std::string input, output, daq_config, gem_ped_file;
    float zerosup_override = 0.f;
    bool prad1 = false;

    std::string db_dir = DATABASE_DIR;
    if (const char *env = std::getenv("PRAD2_DATABASE_DIR"))  db_dir = env;
    daq_config = db_dir + "/daq_config.json"; // default DAQ config for PRad2

    int opt;
    while ((opt = getopt(argc, argv, "o:D:g:z:p")) != -1) {
        switch (opt) {
            case 'o': output = optarg; break;
            case 'D': daq_config = optarg; break;
            case 'p': prad1 = true; break;
            case 'g': gem_ped_file = optarg; break;
            case 'z': zerosup_override = std::atof(optarg); break;
        }
    }
    if (optind < argc) input = argv[optind];

    if (input.empty()) {
        std::cerr << "Usage: replay_recon <input.evio> [-o output.root] [-D daq_config.json] [-p]\n";
        return 1;
    }

    if (output.empty()) {
        output = input;
        auto pos = output.find(".evio");
        if (pos != std::string::npos) output = output.substr(0, pos);
        output += "_recon.root";
    }

    analysis::Replay replay;
    if (!daq_config.empty()) replay.LoadDaqConfig(daq_config);

    if (!replay.ProcessWithRecon(input, output, daq_config, gem_ped_file, zerosup_override, prad1))
        return 1;

    return 0;
}