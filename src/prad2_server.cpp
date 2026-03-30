// =========================================================================
// prad2_server — CLI entry point for PRad-II event viewer/monitor server
//
// Usage:
//   prad2_server [evio_file] [-p port] [-H] [-c config.json]
//                [-d data_dir] [-D daq_config.json] [--et]
//
// Examples:
//   prad2_server data.evio -H              # view file with histograms
//   prad2_server -d /data/stage6 -H        # browse and pick files
//   prad2_server --et                      # online ET monitoring
//   prad2_server data.evio --et -H         # file mode, ET available via toggle
// =========================================================================

#include "viewer_server.h"

#include <iostream>
#include <string>
#include <csignal>
#include <cstdlib>
#include <getopt.h>

#ifndef DATABASE_DIR
#define DATABASE_DIR "."
#endif
#ifndef RESOURCE_DIR
#define RESOURCE_DIR "."
#endif

static ViewerServer *g_server = nullptr;

int main(int argc, char *argv[])
{
    ViewerServer::Config cfg;
    cfg.database_dir = DATABASE_DIR;
    cfg.resource_dir = RESOURCE_DIR;
    if (const char *env = std::getenv("PRAD2_DATABASE_DIR"))  cfg.database_dir = env;
    if (const char *env = std::getenv("PRAD2_RESOURCE_DIR"))  cfg.resource_dir = env;

    static struct option long_opts[] = {
        {"port",       required_argument, nullptr, 'p'},
        {"hist",       no_argument,       nullptr, 'H'},
        {"config",     required_argument, nullptr, 'c'},
        {"data-dir",   required_argument, nullptr, 'd'},
        {"daq-config", required_argument, nullptr, 'D'},
        {"et",         no_argument,       nullptr, 'E'},
        {"help",       no_argument,       nullptr, '?'},
        {nullptr, 0, nullptr, 0},
    };

    int opt;
    while ((opt = getopt_long(argc, argv, "p:Hc:d:D:", long_opts, nullptr)) != -1) {
        switch (opt) {
        case 'p': cfg.port = std::atoi(optarg); break;
        case 'H': cfg.hist_enabled = true; break;
        case 'c': cfg.config_file = optarg; break;
        case 'd': cfg.data_dir = optarg; break;
        case 'D': cfg.daq_config_file = optarg; break;
        case 'E': cfg.start_online = true; break;
        default:
            std::cerr << "Usage: " << argv[0]
                      << " [evio_file] [-p port] [-H] [-c config.json]"
                      << " [-d data_dir] [-D daq_config.json] [--et]\n";
            return 1;
        }
    }
    if (optind < argc) cfg.initial_file = argv[optind];

    ViewerServer server;
    g_server = &server;

    std::signal(SIGINT, [](int) {
        if (g_server) g_server->stop();
    });

    server.init(cfg);
    server.run();

    return 0;
}
