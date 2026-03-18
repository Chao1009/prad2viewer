// test/evchan_test.cpp — Test consistency between EvChannel and EtChannel
//
// Start this program first (connects to ET), then use et_feeder to feed
// the same evio file to the same ET system. Set interval large enough
// to keep the two channels synchronized.
//
// Usage: evchan_test <evio_file> [-h host] [-p port] [-f et_file] [-i interval_ms]

#include "EtConfigWrapper.h"
#include "EvChannel.h"
#include "EtChannel.h"
#include <csignal>
#include <thread>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <unistd.h>

#define PROGRESS_COUNT 100

using namespace std::chrono;

volatile std::sig_atomic_t gSignalStatus;

void signal_handler(int signal) { gSignalStatus = signal; }

static void usage(const char *prog) {
    std::cerr << "Usage: " << prog << " <evio_file> [-h host] [-p port] [-f et_file] [-i interval_ms]\n";
}

int main(int argc, char* argv[])
{
    std::string host = "localhost";
    int port = 11111;
    std::string et_file = "/tmp/et_feeder";
    int interval = 100;

    int opt;
    while ((opt = getopt(argc, argv, "h:p:f:i:")) != -1) {
        switch (opt) {
        case 'h': host = optarg; break;
        case 'p': port = std::atoi(optarg); break;
        case 'f': et_file = optarg; break;
        case 'i': interval = std::atoi(optarg); break;
        default:  usage(argv[0]); return 1;
        }
    }
    if (optind >= argc) { usage(argv[0]); return 1; }
    std::string evio_file = argv[optind];

    // ET channel reader
    evc::EtChannel et_chan;
    if (et_chan.Connect(host, port, et_file) != evc::status::success ||
        et_chan.Open("MONITOR") != evc::status::success) {
        std::cerr << "Failed to open ET channel\n";
        return -1;
    }

    // evio file reader
    evc::EvChannel ev_chan;
    if (ev_chan.Open(evio_file) != evc::status::success) {
        std::cerr << "Failed to open coda file \"" << evio_file << "\"\n";
        return -1;
    }

    // install signal handler
    std::signal(SIGINT, signal_handler);
    int count = 0;
    bool loop = true;
    std::cout << "Listening to ET system, now you can feed the data to it.\n";
    while (loop) {
        if (gSignalStatus == SIGINT) {
            std::cout << "Received control-C, exiting...\n";
            break;
        }
        switch (et_chan.Read()) {
        case evc::status::success:
            break;
        case evc::status::empty:
            std::this_thread::sleep_for(std::chrono::seconds(1));
            continue;
        default:
            loop = false;
            continue;
        }
        ev_chan.Read();
        system_clock::time_point start(system_clock::now());
        system_clock::time_point next(start + std::chrono::milliseconds(interval));

        if (++count % PROGRESS_COUNT == 0) {
            std::cout << "Received " << count << " events.\r" << std::flush;
        }
        auto ev_buf = ev_chan.GetRawBuffer();
        auto et_buf = et_chan.GetRawBuffer();
        size_t et_len = et_buf[0] + 1;
        size_t ev_len = ev_buf[0] + 1;
        std::cout << "Event " << count << ": ET len=" << et_len
                  << ", EV len=" << ev_len << "\n";

        std::cout << std::hex << std::setfill('0');
        for (size_t i = 0; i < std::max(ev_len, et_len); ++i) {
            std::cout << "0x" << std::setw(8);
            if (i < et_len) { std::cout << et_buf[i]; }
            if (i < ev_len) { std::cout << ", 0x" << std::setw(8) << ev_buf[i]; }
            std::cout << "\n";
        }
        std::cout << std::dec;
        std::this_thread::sleep_until(next);
    }
    std::cout << "Received " << count << " events\n";

    et_chan.Close();
    et_chan.Disconnect();
    ev_chan.Close();
    return 0;
}
