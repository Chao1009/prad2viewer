// test/et_feeder.cpp — Feed an evio file to an ET system event-by-event
//
// Usage: et_feeder <evio_file> [-h host] [-p port] [-f et_file] [-i interval_ms] [-s start] [-n num]

#include "EtConfigWrapper.h"
#include "EvChannel.h"
#include <csignal>
#include <thread>
#include <chrono>
#include <iostream>
#include <cstring>
#include <cstdlib>
#include <unistd.h>

#define PROGRESS_COUNT 100

using namespace std::chrono;

volatile std::sig_atomic_t gSignalStatus;

void signal_handler(int signal) { gSignalStatus = signal; }

static void usage(const char *prog) {
    std::cerr << "Usage: " << prog << " <evio_file> [options]\n\n"
              << "Options:\n"
              << "  -h <host>     ET host (default: localhost)\n"
              << "  -p <port>     ET port (default: 11111)\n"
              << "  -f <file>     ET system file (default: /tmp/et_feeder)\n"
              << "  -i <ms>       Interval between events in ms (default: 100)\n"
              << "  -s <N>        Start event number, 1-based (default: 1)\n"
              << "  -n <N>        Number of events to feed (default: all)\n";
}

int main(int argc, char* argv[])
{
    std::string host = "localhost";
    int port = 11111;
    std::string et_file = "/tmp/et_feeder";
    int interval = 100;
    int start_ev = 1;
    int max_events = 0;  // 0 = all

    int opt;
    while ((opt = getopt(argc, argv, "h:p:f:i:s:n:")) != -1) {
        switch (opt) {
        case 'h': host = optarg; break;
        case 'p': port = std::atoi(optarg); break;
        case 'f': et_file = optarg; break;
        case 'i': interval = std::atoi(optarg); break;
        case 's': start_ev = std::atoi(optarg); break;
        case 'n': max_events = std::atoi(optarg); break;
        default:  usage(argv[0]); return 1;
        }
    }
    if (optind >= argc) { usage(argv[0]); return 1; }
    std::string evio_file = argv[optind];

    et_sys_id et_id;
    et_att_id att_id;

    // open ET system
    et_wrap::OpenConfig conf;
    conf.set_cast(ET_DIRECT);
    conf.set_host(host.c_str());
    conf.set_serverport(port);

    std::vector<char> fname(et_file.begin(), et_file.end());
    fname.push_back('\0');
    auto status = et_open(&et_id, fname.data(), conf.configure().get());

    if (status != ET_OK) {
        std::cerr << "Cannot open ET at " << host << ":" << port << " with " << et_file << "\n";
        return -1;
    }

    // attach to GRAND CENTRAL
    status = et_station_attach(et_id, ET_GRANDCENTRAL, &att_id);
    if (status != ET_OK) {
        std::cerr << "Failed to attach to the ET Grand Central Station.\n";
        return -1;
    }

    // evio file reader
    evc::EvChannel chan;
    if (chan.Open(evio_file) != evc::status::success) {
        std::cerr << "Failed to open coda file \"" << evio_file << "\"\n";
        return -1;
    }

    // install signal handler
    std::signal(SIGINT, signal_handler);
    int total = 0, fed = 0;
    et_event *ev;
    while ((chan.Read() == evc::status::success) && et_alive(et_id)) {
        if (gSignalStatus == SIGINT) {
            std::cout << "Received control-C, exiting...\n";
            break;
        }
        ++total;

        // skip to start event
        if (total < start_ev) continue;

        // check max events
        if (max_events > 0 && fed >= max_events) break;

        system_clock::time_point t0(system_clock::now());
        system_clock::time_point next(t0 + std::chrono::milliseconds(interval));

        if (++fed % PROGRESS_COUNT == 0) {
            std::cout << "Fed " << fed << " events (from #" << start_ev << ") to ET.\r" << std::flush;
        }

        uint32_t *buf = chan.GetRawBuffer();
        size_t nbytes = (buf[0] + 1) * sizeof(uint32_t);

        status = et_event_new(et_id, att_id, &ev, ET_SLEEP, nullptr, nbytes);
        if (status != ET_OK) {
            std::cerr << "Failed to add new event to the ET system.\n";
            return -1;
        }
        void *data;
        et_event_getdata(ev, &data);
        memcpy(data, buf, nbytes);
        et_event_setlength(ev, nbytes);

        status = et_event_put(et_id, att_id, ev);
        if (status != ET_OK) {
            std::cerr << "Failed to put event back to the ET system.\n";
            return -1;
        }

        std::this_thread::sleep_until(next);
    }
    std::cout << "Fed " << fed << " events (starting from #" << start_ev << ") to ET\n";

    chan.Close();
    return 0;
}
