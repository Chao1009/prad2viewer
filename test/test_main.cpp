// test/test_main.cpp
// Usage: evc_test <evio_file> [max_events] [-v]

#include "EvChannel.h"
#include "Fadc250Data.h"
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstring>

using namespace evc;

int main(int argc, char *argv[])
{
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <evio_file> [max_events] [-v]\n";
        return 1;
    }

    const char *filename = argv[1];
    int max_ev = 0;
    bool verbose = false;
    for (int a = 2; a < argc; ++a) {
        if (std::strcmp(argv[a], "-v") == 0) verbose = true;
        else max_ev = std::atoi(argv[a]);
    }

    EvChannel ch;
    if (ch.Open(filename) != status::success) {
        std::cerr << "Failed to open " << filename << "\n";
        return 1;
    }

    fdec::EventData event;
    int total = 0, nread = 0;

    while (ch.Read() == status::success) {
        ++nread;
        if (!ch.Scan()) {
            std::cerr << "Read " << nread << ": Scan failed\n";
            continue;
        }

        auto hdr = ch.GetEvHeader();
        std::cout << "Read " << nread
                  << "  tag=0x" << std::hex << hdr.tag << std::dec
                  << "  type=0x" << std::hex << hdr.type << std::dec
                  << "  num=" << hdr.num
                  << "  length=" << hdr.length
                  << "  nevents=" << ch.GetNEvents() << "\n";

        if (ch.GetNEvents() == 0) continue;

        int nevt = ch.GetNEvents();
        for (int i = 0; i < nevt; ++i) {
            if (!ch.DecodeEvent(i, event)) continue;
            ++total;

            std::cout << "=== Event " << total << ": " << event.nrocs << " ROCs ===\n";
            for (int r = 0; r < event.nrocs; ++r) {
                auto &roc = event.rocs[r];
                std::cout << "  ROC 0x" << std::hex << roc.tag << std::dec
                          << " (" << roc.nslots << " slots)\n";

                for (int s = 0; s < fdec::MAX_SLOTS; ++s) {
                    if (!roc.slots[s].present) continue;
                    auto &slot = roc.slots[s];

                    std::cout << "    slot " << std::setw(2) << s
                              << "  ev=" << slot.trigger
                              << "  ts=0x" << std::hex << (uint64_t)slot.timestamp << std::dec
                              << "  nch=" << slot.nchannels << "\n";

                    for (int c = 0; c < fdec::MAX_CHANNELS; ++c) {
                        if (!(slot.channel_mask & (1u << c))) continue;
                        auto &cd = slot.channels[c];

                        if (verbose) {
                            std::cout << "      ch " << std::setw(2) << c
                                      << " [" << cd.nsamples << "]:";
                            for (int j = 0; j < cd.nsamples; ++j)
                                std::cout << " " << cd.samples[j];
                            std::cout << "\n";
                        } else {
                            std::cout << "      ch " << std::setw(2) << c
                                      << " [" << cd.nsamples << "]\n";
                        }
                    }
                }
            }

            if (max_ev > 0 && total >= max_ev) goto done;
        }
    }
done:
    std::cout << "Read " << nread << " buffers, decoded " << total << " event(s).\n";
    ch.Close();
    return 0;
}
