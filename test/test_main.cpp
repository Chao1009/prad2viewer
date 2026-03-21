// test/test_main.cpp
// Basic smoke-test for the evc library.
//
// Usage:
//   evc_test <evio_file>                          read buffers
//   evc_test <evio_file> -m scan [-s N] [-n N]    per-event details
//   evc_test -m et -H <host> -P <port> -f <et_file> -S <station>

#include "EvChannel.h"
#include "EtChannel.h"
#include "EvStruct.h"
#include "Fadc250Data.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <cstdlib>
#include <cstdint>
#include <getopt.h>

static void usage(const char *prog)
{
    std::cerr << "Usage:\n"
              << "  " << prog << " <evio_file>                                  Read buffers\n"
              << "  " << prog << " <evio_file> -m scan [-s start] [-n num]      Per-event details\n"
              << "  " << prog << " -m et -H <host> -P <port> -f <file> -S <station>  ET system\n";
}

// ---- scan mode: per-event details + duplicate detection ------------------
static int scanFile(const std::string &path, int start_ev, int num_ev)
{
    evc::EvChannel ch;
    if (ch.Open(path) != evc::status::success) {
        std::cerr << "Failed to open: " << path << "\n"; return 1;
    }

    fdec::EventData event;
    int ev_seq = 0;      // decoded event counter (matches monitor seq)
    int printed = 0;
    uint32_t prev_hash = 0;
    int buf_num = 0;

    while (ch.Read() == evc::status::success) {
        ++buf_num;
        if (!ch.Scan()) continue;

        int nevt = ch.GetNEvents();
        for (int i = 0; i < nevt; ++i) {
            ++ev_seq;
            if (ev_seq < start_ev) continue;
            if (printed >= num_ev) break;

            uint32_t *raw = ch.GetRawBuffer();
            auto hdr = evc::BankHeader(raw);

            // simple hash of data words (detect identical content)
            uint32_t hash = 0;
            int words = std::min((int)hdr.length + 1, 256);
            for (int j = 0; j < words; ++j) hash ^= raw[j] * 2654435761u;

            bool dup = (printed > 0 && hash == prev_hash);

            // decode this sub-event
            int nrocs = 0, nchannels = 0;
            if (ch.DecodeEvent(i, event)) {
                for (int r = 0; r < event.nrocs; ++r) {
                    if (!event.rocs[r].present) continue;
                    nrocs++;
                    for (int s = 0; s < fdec::MAX_SLOTS; ++s) {
                        auto &slot = event.rocs[r].slots[s];
                        if (!slot.present) continue;
                        for (int c = 0; c < fdec::MAX_CHANNELS; ++c)
                            if (slot.channel_mask & (1ull << c)) nchannels++;
                    }
                }
            }

            std::cout << "ev " << std::setw(5) << ev_seq
                      << "  buf=" << std::setw(4) << buf_num
                      << "  sub=" << i << "/" << nevt
                      << "  tag=0x" << std::hex << std::setw(4) << std::setfill('0') << hdr.tag
                      << std::dec << std::setfill(' ')
                      << "  num=" << std::setw(4) << hdr.num
                      << "  len=" << std::setw(8) << hdr.length
                      << "  rocs=" << nrocs
                      << "  ch=" << std::setw(4) << nchannels
                      << "  hash=" << std::hex << hash << std::dec
                      << (dup ? "  ** DUP" : "")
                      << "\n";

            prev_hash = hash;
            ++printed;
        }
        if (printed >= num_ev) break;
    }
    ch.Close();
    std::cout << "Printed " << printed << " events (seq " << start_ev << "-" << (start_ev + printed - 1) << ")\n";
    return 0;
}

// ---- file mode -----------------------------------------------------------
static int testFile(const std::string &path)
{
    evc::EvChannel ch;

    if (ch.Open(path) != evc::status::success) {
        std::cerr << "Failed to open: " << path << "\n";
        return 1;
    }

    int nevents = 0;
    evc::status st;
    while ((st = ch.Read()) == evc::status::success) {
        auto hdr = ch.GetEvHeader();
        std::cout << "Event " << ++nevents
                  << "  tag=" << hdr.tag
                  << "  type=" << hdr.type
                  << "  length=" << hdr.length << "\n";
    }

    std::cout << "Done. Read " << nevents << " event(s). Final status: "
              << static_cast<int>(st) << "\n";
    ch.Close();
    return 0;
}

// ---- ET mode -------------------------------------------------------------
static int testET(const std::string &ip, int port,
                  const std::string &et_file, const std::string &station)
{
    evc::EtChannel ch;

    if (ch.Connect(ip, port, et_file) != evc::status::success) {
        std::cerr << "Failed to connect to ET at " << ip << ":" << port << "\n";
        return 1;
    }
    if (ch.Open(station) != evc::status::success) {
        std::cerr << "Failed to open station: " << station << "\n";
        ch.Disconnect();
        return 1;
    }

    int nevents = 0, max_events = 20;
    evc::status st;
    while (nevents < max_events && (st = ch.Read()) != evc::status::failure) {
        if (st == evc::status::empty) continue;
        auto hdr = ch.GetEvHeader();
        std::cout << "ET Event " << ++nevents
                  << "  tag=" << hdr.tag
                  << "  length=" << hdr.length << "\n";
    }

    std::cout << "Done. Read " << nevents << " event(s).\n";
    ch.Disconnect();
    return 0;
}

// --------------------------------------------------------------------------
int main(int argc, char *argv[])
{
    std::string mode;
    std::string et_host = "localhost", et_file, et_station;
    int et_port = 11111;
    int start = 1, num = 50;

    int opt;
    while ((opt = getopt(argc, argv, "m:s:n:H:P:f:S:h")) != -1) {
        switch (opt) {
        case 'm': mode = optarg; break;
        case 's': start = std::atoi(optarg); break;
        case 'n': num = std::atoi(optarg); break;
        case 'H': et_host = optarg; break;
        case 'P': et_port = std::atoi(optarg); break;
        case 'f': et_file = optarg; break;
        case 'S': et_station = optarg; break;
        default:  usage(argv[0]); return 1;
        }
    }

    if (mode == "et") {
        if (et_file.empty() || et_station.empty()) { usage(argv[0]); return 1; }
        return testET(et_host, et_port, et_file, et_station);
    }

    if (optind >= argc) { usage(argv[0]); return 1; }
    std::string path = argv[optind];

    if (mode == "scan") return scanFile(path, start, num);
    return testFile(path);
}
