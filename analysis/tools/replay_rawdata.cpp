//=============================================================================
// replay_rawdata — convert multiple EVIO files to ROOT trees (multi-threaded)
//
// Usage: replay_rawdata <evio_dir> -o output_dir [-f max_files] [-n max_events] [-p] [-j num_threads]
//   -o  output directory (REQUIRED)
//   -f  max files to process (default: all)
//   -n  max events per file (default: all)
//   -p  include peak analysis branches
//   -j  number of threads (default: 4)
//   -D  DAQ configuration file
//=============================================================================

#include "Replay.h"
#include "InstallPaths.h"

#include <iostream>
#include <string>
#include <cstdlib>
#include <getopt.h>
#include <filesystem>
#include <algorithm>
#include <vector>
#include <thread>
#include <atomic>
#include <mutex>

#include <TFileMerger.h>
#include <TClass.h>
#include <TROOT.h>

#ifndef DATABASE_DIR
#define DATABASE_DIR "."
#endif

static std::vector<std::string> getFilesInDir(const std::string &dir_path)
{
    std::vector<std::string> files;
    for (auto &entry : std::filesystem::directory_iterator(dir_path)) {
        if (entry.is_regular_file()) {
            if (entry.path().filename().string().find(".evio") != std::string::npos)
                files.push_back(entry.path().string());
        }
    }
    std::sort(files.begin(), files.end());
    return files;
}

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
    // Initialize ROOT for multi-threading
    ROOT::EnableThreadSafety();
    
    // Force ROOT dictionary initialization in main thread
    // This prevents concurrent TClass::Init calls
    TClass::GetClass("TTree");
    TClass::GetClass("TFile");
    TClass::GetClass("TBranch");

    std::string input, daq_config, merged_output, output_dir;
    int max_events = -1;
    int max_files = -1;
    bool peaks = false;
    int num_threads = 4;

    std::string db_dir = prad2::resolve_data_dir(
        "PRAD2_DATABASE_DIR",
        {"../share/prad2evviewer/database"},
        DATABASE_DIR);
    daq_config = db_dir + "/daq_config.json"; // default DAQ config for PRad2

    int opt;
    while ((opt = getopt(argc, argv, "o:f:n:D:j:p")) != -1) {
        switch (opt) {
            case 'o': output_dir = optarg; break;
            case 'f': max_files = std::atoi(optarg); break;
            case 'n': max_events = std::atoi(optarg); break;
            case 'D': daq_config = optarg; break;
            case 'j': num_threads = std::atoi(optarg); break;
            case 'p': peaks = true; break;
        }
    }
    if (optind < argc) input = argv[optind];

    if (input.empty() || output_dir.empty()) {
        std::cerr << "Usage: replay_rawdata <evio_dir> -o output_dir [-f max_files] [-j threads]"
                  << " [-D daq_config.json] [-n N] [-p]\n";
        std::cerr << "  -o  output directory (REQUIRED)\n";
        std::cerr << "  -f  max files to process (default: all)\n";
        std::cerr << "  -j  number of threads (default: 4)\n";
        std::cerr << "  -D  DAQ config JSON (default: <db>/daq_config.json)\n";
        std::cerr << "  -n  max events per file (default: all)\n";
        std::cerr << "  -p  include peak analysis branches\n";
        return 1;
    }

    std::vector<std::string> evio_files = getFilesInDir(input);
    if (evio_files.empty()) {
        std::cerr << "No EVIO files found in: " << input << "\n";
        return 1;
    }
    int num_files = static_cast<int>(evio_files.size());
    if (max_files > 0) num_files = std::min(num_files, max_files);
    num_threads = std::max(1, std::min(num_threads, num_files));

    std::cout << "Processing " << num_files << " files with "
              << num_threads << " threads\n";

    // shared work queue: atomic index into file list
    std::atomic<int> next_file{0};
    std::mutex io_mtx;
    std::atomic<int> errors{0};
    std::vector<std::string> merged_files;

    auto worker = [&]() {
        // each thread gets its own Replay instance (own EvChannel, own buffers)
        analysis::Replay replay;
        if (!daq_config.empty()) replay.LoadDaqConfig(daq_config);
        replay.LoadDaqMap(db_dir + "/daq_map.json");
        std::cerr << "Using DAQ map: " << db_dir + "/daq_map.json" << "\n";

        while (true) {
            int idx = next_file.fetch_add(1);
            if (idx >= num_files) break;

            std::string out = output_dir + "/" + makeOutputFile(evio_files[idx]);
            bool ok = replay.Process(evio_files[idx], out, max_events, peaks, daq_config);

            std::lock_guard<std::mutex> lk(io_mtx);
            if (ok) {
                std::cout << "  [" << (idx + 1) << "/" << num_files << "] "
                          << evio_files[idx] << " -> " << out << "\n";
            } else {
                std::cerr << "  [" << (idx + 1) << "/" << num_files << "] FAILED: "
                          << evio_files[idx] << "\n";
                errors++;
            }
        }
    };

    std::vector<std::thread> threads;
    threads.reserve(num_threads);
    for (int i = 0; i < num_threads; ++i)
        threads.emplace_back(worker);
    for (auto &t : threads)
        t.join();

    std::cout << "Done: " << num_files << " files"
              << (errors > 0 ? ", " + std::to_string(errors.load()) + " errors" : "")
              << "\n";

    return errors > 0 ? 1 : 0;
}
