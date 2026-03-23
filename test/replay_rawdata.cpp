
#include "replay.h"
static void usage(const char *prog)
{
    std::cerr << "Usage:\n"
              << "  " << prog << " <evio_file> <output_root_file> [--num N] [--peaks]   -- decode events (default: all, mode: raw samples)\n"
              << "  Options:\n"
              << "    --num N    max number of events to decode\n"
              << "    --peaks    save peak info instead of raw samples\n";
}

int main(int argc, char *argv[]){
    if (argc < 3) { usage(argv[0]); return 1; }
    std::string evio_file = argv[1];
    std::string output_root_file = argv[2];
    int max_events = -1;  // default: decode all events
    bool mode_peaks = false;  // default: save raw samples
    for (int i = 3; i < argc; i++) {
        if (std::string(argv[i]) == "--num" && i + 1 < argc) {
            max_events = std::stoi(argv[++i]);
        } 
        else if (std::string(argv[i]) == "--peaks") {
            mode_peaks = true;
        } else {
            std::cerr << "Unknown option: " << argv[i] << "\n";
            usage(argv[0]);
            return 1;
        }
    }
    replay rep;
    if(!rep.rawdata2root(evio_file, output_root_file, max_events, mode_peaks)) {
        std::cerr << "Failed to convert raw data to ROOT file.\n";
        return 1;
    }
    return 0;
}
