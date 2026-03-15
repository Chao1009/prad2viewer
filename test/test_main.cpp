// test/test_main.cpp
// Usage: evc_test <evio_file> [max_events]

#include "EvChannel.h"
#include "Fadc250Decoder.h"
#include <iostream>
#include <cstdlib>

using namespace evc;

int main(int argc, char *argv[])
{
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <evio_file> [max_events]\n";
        return 1;
    }
    int max_ev = (argc >= 3) ? std::atoi(argv[2]) : 0;

    EvChannel ch;
    if (ch.Open(argv[1]) != status::success) {
        std::cerr << "Failed to open " << argv[1] << "\n";
        return 1;
    }

    int n = 0;
    while (ch.Read() == status::success) {
        ++n;
        if (!ch.Scan()) continue;

        auto hdr = ch.GetEvHeader();
        std::cout << "=== Event " << n
                  << "  tag=0x" << std::hex << hdr.tag << std::dec
                  << "  length=" << hdr.length << " ===\n";

        // print full bank tree
        ch.PrintTree(std::cout);

        // decode all FADC250 composite banks (tag 0xe101)
        for (auto *node : ch.FindByTag(0xe101)) {
            size_t nbytes;
            auto *payload = ch.GetCompositePayload(*node, nbytes);
            if (!payload) continue;

            auto slots = fdec::Fadc250Decoder::Decode(payload, nbytes);
            for (auto &s : slots) {
                std::cout << "    FADC slot=" << (int)s.slot
                          << " trig=" << s.trigger
                          << " ts=0x" << std::hex << (uint64_t)s.timestamp << std::dec
                          << " |";
                for (auto &c : s.channels)
                    std::cout << " ch" << (int)c.channel << ":" << c.samples.size();
                std::cout << "\n";
            }
        }
        std::cout << "\n";
        if (max_ev > 0 && n >= max_ev) break;
    }

    std::cout << "Done. " << n << " event(s).\n";
    ch.Close();
    return 0;
}
