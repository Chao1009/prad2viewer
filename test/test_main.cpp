// src/test_main.cpp
// Basic smoke-test for the evc library.
// Usage:
//   evc_test <evio_file>            -- read from file
//   evc_test --et <ip> <port> <file> <station>  -- read from ET system

#include "EvChannel.h"
#include "EtChannel.h"
#include <iostream>
#include <string>
#include <cstdlib>

static void usage(const char *prog)
{
    std::cerr << "Usage:\n"
              << "  " << prog << " <evio_file>\n"
              << "  " << prog << " --et <ip> <port> <et_file> <station>\n";
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
    if (argc < 2) { usage(argv[0]); return 1; }

    std::string first = argv[1];

    if (first == "--et") {
        if (argc < 6) { usage(argv[0]); return 1; }
        return testET(argv[2], std::atoi(argv[3]), argv[4], argv[5]);
    }

    return testFile(first);
}
