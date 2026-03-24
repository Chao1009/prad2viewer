//=============================================================================
// replay_rawdata — convert EVIO file to ROOT tree
//Multithreaded version of replay_rawdata
// Usage: replay_rawdata <input.evio> [-f files_number] [-n max_events] [-p] [-j num_threads]
//   -n  max events to process (default: all)
//   -p  include peak analysis branches
//   -j  number of threads to use (default: 1)
//=============================================================================

#include "Replay.h"

#include <iostream>
#include <string>
#include <getopt.h>
#include <filesystem>
#include <algorithm>
#include <vector>

#include <unistd.h>
#include <sys/wait.h>

std::vector<std::string> getFilesInDir(const std::string &dir_path)
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

int main(int argc, char *argv[])
{
    std::string input, daq_config;
    int max_events = -1;
    int max_files = -1;
    bool peaks = false;

    int num_threads = 1;

    int opt;
    while ((opt = getopt(argc, argv, "f:n:D:j:p")) != -1) {
        switch (opt) {
            case 'f': max_files = std::atoi(optarg); break;
            case 'n': max_events = std::atoi(optarg); break;
            case 'D': daq_config = optarg; break;
            case 'j': num_threads = std::atoi(optarg); break;
            case 'p': peaks = true; break;
        }
    }
    if (optind < argc) input = argv[optind];

    if (input.empty()) {
        std::cerr << "Usage: replay_rawdata <evio_dir> [-f <files_number>] [-j <num_threads>] [-D daq_config.json] [-n N] [-p]\n";
        return 1;
    }

    // Get list of EVIO files in the directory
    std::vector<std::string> evio_files = getFilesInDir(input);
    int num_files = evio_files.size();
    if (max_files > 0) {
        num_files = std::min(num_files, max_files);
    }
    if(num_files == 0){
        std::cerr << "No files found in directory: " << input << "\n";
        return 1;
    }
    // Limit number of threads to number of files
    num_threads = std::min(num_threads, num_files);
    std::vector<pid_t> pids(num_threads);
    int files_per_process = (num_files + num_threads - 1) / num_threads; // round up

    std::cout << "Processing " << num_files << " files with " << num_threads << " threads (" 
              << files_per_process << " files/thread)\n";

    for (int i = 0; i < num_threads; ++i) {
        int start = i * files_per_process;
        int end   = std::min(start + files_per_process, num_files);

        pid_t pid = fork();
        if (pid < 0) {
            std::cerr << "process " << i << " fork failed\n";
            continue;
        }
        else if (pid == 0) {
            // child process: handle files from start to end
            analysis::Replay replay;
            if (!daq_config.empty()) replay.LoadDaqConfig(daq_config);

            for (int f = start; f < end; ++f) {
                // each child process writes to its own output file to avoid conflicts
                std::string out = evio_files[f];
                auto pos = out.find(".evio");
                if (pos != std::string::npos) out = 
                    out.substr(0, pos) + out.substr(pos + 5);
                out += ".root";
                if(replay.Process(evio_files[f], out, max_events, peaks))
                    std::cout << "Processed " << evio_files[f] << " -> " << out << "\n";
                else{
                    std::cerr << "Failed to process " << evio_files[f] << "\n";
                    std::exit(1);
                }
            }
            std::exit(0);
        }
        else pids[i] = pid;
    }

    // parent process waits for all child processes to finish
    for (int i = 0; i < num_threads; ++i)
        if(pids[i] > 0) waitpid(pids[i], nullptr, 0);

    return 0;
}
