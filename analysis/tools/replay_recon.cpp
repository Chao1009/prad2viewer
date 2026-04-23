//=============================================================================
// replay_recon — decode EVIO file and do detector reconstruction write to ROOT tree
//
// Usage: replay_recon <input.evio> -o output_dir [-D daq_config.json] [-p]
//                                  [-g gem_pedestal.json] [-z zerosup_threshold]
//   -o  output directory (REQUIRED)
//   -D  DAQ config file (default: PRAD2_DATABASE_DIR/daq_config.json)
//   -p  read prad1 data and do not include GEM
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
    out += "_recon.root";
    return out;
}

int main(int argc, char *argv[])
{
    std::string input, output_dir, daq_config, gem_ped_file;
    float zerosup_override = 0.f;
    bool prad1 = false;

    std::string db_dir = prad2::resolve_data_dir(
        "PRAD2_DATABASE_DIR",
        {"../share/prad2evviewer/database"},
        DATABASE_DIR);
    daq_config = db_dir + "/daq_config.json"; // default DAQ config for PRad2

    int opt;
    while ((opt = getopt(argc, argv, "o:D:g:z:p")) != -1) {
        switch (opt) {
            case 'o': output_dir = optarg; break;
            case 'D': daq_config = optarg; break;
            case 'p': prad1 = true; break;
            case 'g': gem_ped_file = optarg; break;
            case 'z': zerosup_override = std::atof(optarg); break;
        }
    }
    if (optind < argc) input = argv[optind];

    if (input.empty() || output_dir.empty()) {
        std::cerr << "Usage: replay_recon <input.evio> -o output_dir [-D daq_config.json] [-p]\n";
        std::cerr << "  -o  output directory (REQUIRED)\n";
        std::cerr << "  -D  DAQ config JSON (default: <db>/daq_config.json)\n";
        std::cerr << "  -g  GEM pedestal JSON\n";
        std::cerr << "  -z  zero-suppression threshold override\n";
        std::cerr << "  -p  PRad1 mode (no GEM)\n";
        return 1;
    }

    std::string output_file = output_dir + "/" + makeOutputFile(input);

    analysis::Replay replay;
    if (!daq_config.empty()) replay.LoadDaqConfig(daq_config);
    replay.LoadDaqMap(db_dir + "/daq_map.json");
    std::cerr << "Using DAQ map: " << db_dir + "/daq_map.json" << "\n";

    if (!replay.ProcessWithRecon(input, output_file, daq_config, gem_ped_file, zerosup_override, prad1))
        return 1;

    return 0;
}