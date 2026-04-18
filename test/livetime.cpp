// livetime — calculate DAQ live time from EVIO files
//
// Two independent methods:
//   1. DSC2 scalers: gated (busy) vs ungated (total) from the 0xe115 bank.
//      "Gated" counts only during DAQ busy (dead time), "ungated" counts
//      all the time.  Live time = 1 - gated/ungated.
//      (Same convention as PRadAnalyzer/calcLiveTime.)
//   2. Pulser counting: accepted 100 Hz pulser events vs expected from
//      elapsed time.
//
// Accepts a single file, a base name (auto-discovers .00000, .00001, ...),
// or a directory (processes all .evio* files sorted by name).
//
// Usage: livetime <input> [-D daq_config.json] [-f pulser_freq] [-t interval]

#include "EvChannel.h"
#include "DaqConfig.h"
#include "load_daq_config.h"
#include "Fadc250Data.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <map>
#include <algorithm>
#include <getopt.h>

// portable directory listing
#ifdef _WIN32
#include <windows.h>
#else
#include <dirent.h>
#include <sys/stat.h>
#endif

using namespace evc;

#ifndef DATABASE_DIR
#define DATABASE_DIR "."
#endif

// ---- DSC2 scaler data layout (dsc2ReadScalers with rflag=0xFF) ----
//
// Grp1 = externally gated (gate source = BUSY signal -> counts dead time)
// Grp2 = ungated (gate source = constant -> counts all time)
//
// Per slot: 67 words
//   [0]      header  0xDCA00000 | (slot<<8) | rflag
//   [1..16]  TRG  Grp1 (gated/busy)  — 16 channels
//   [17..32] TDC  Grp1 (gated/busy)  — 16 channels
//   [33..48] TRG  Grp2 (ungated)     — 16 channels
//   [49..64] TDC  Grp2 (ungated)     — 16 channels
//   [65]     Ref  Grp1 (gated/busy)  — 125 MHz clock
//   [66]     Ref  Grp2 (ungated)     — 125 MHz clock

static constexpr uint32_t DSC2_HEADER_MASK = 0xFFFF0000;
static constexpr uint32_t DSC2_HEADER_ID   = 0xDCA00000;
static constexpr int DSC2_WORDS_PER_SLOT   = 67;
static constexpr int DSC2_NCH              = 16;

struct Dsc2Slot {
    uint32_t slot;
    uint32_t trg_gated[DSC2_NCH];
    uint32_t tdc_gated[DSC2_NCH];
    uint32_t trg_ungated[DSC2_NCH];
    uint32_t tdc_ungated[DSC2_NCH];
    uint32_t ref_gated;
    uint32_t ref_ungated;
};

static std::vector<Dsc2Slot> parse_dsc2_bank(const uint32_t *data, size_t nwords)
{
    std::vector<Dsc2Slot> slots;
    size_t pos = 0;
    while (pos + DSC2_WORDS_PER_SLOT <= nwords) {
        uint32_t hdr = data[pos];
        if ((hdr & DSC2_HEADER_MASK) != DSC2_HEADER_ID)
            break;
        Dsc2Slot s{};
        s.slot = (hdr >> 8) & 0xFF;
        const uint32_t *p = &data[pos + 1];
        std::memcpy(s.trg_gated,   p,      DSC2_NCH * 4);
        std::memcpy(s.tdc_gated,   p + 16, DSC2_NCH * 4);
        std::memcpy(s.trg_ungated, p + 32, DSC2_NCH * 4);
        std::memcpy(s.tdc_ungated, p + 48, DSC2_NCH * 4);
        s.ref_gated   = p[64];
        s.ref_ungated = p[65];
        slots.push_back(s);
        pos += DSC2_WORDS_PER_SLOT;
    }
    return slots;
}

// ---- file discovery ----------------------------------------------------------

static bool is_regular_file(const std::string &path)
{
#ifdef _WIN32
    DWORD attr = GetFileAttributesA(path.c_str());
    return attr != INVALID_FILE_ATTRIBUTES && !(attr & FILE_ATTRIBUTE_DIRECTORY);
#else
    struct stat st;
    return stat(path.c_str(), &st) == 0 && S_ISREG(st.st_mode);
#endif
}

static bool is_directory(const std::string &path)
{
#ifdef _WIN32
    DWORD attr = GetFileAttributesA(path.c_str());
    return attr != INVALID_FILE_ATTRIBUTES && (attr & FILE_ATTRIBUTE_DIRECTORY);
#else
    struct stat st;
    return stat(path.c_str(), &st) == 0 && S_ISDIR(st.st_mode);
#endif
}

static std::vector<std::string> list_dir(const std::string &dir)
{
    std::vector<std::string> entries;
#ifdef _WIN32
    WIN32_FIND_DATAA fd;
    HANDLE h = FindFirstFileA((dir + "\\*").c_str(), &fd);
    if (h == INVALID_HANDLE_VALUE) return entries;
    do {
        if (!(fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY))
            entries.push_back(dir + "/" + fd.cFileName);
    } while (FindNextFileA(h, &fd));
    FindClose(h);
#else
    DIR *d = opendir(dir.c_str());
    if (!d) return entries;
    while (auto *e = readdir(d)) {
        if (e->d_name[0] == '.') continue;
        std::string full = dir + "/" + e->d_name;
        if (is_regular_file(full)) entries.push_back(full);
    }
    closedir(d);
#endif
    return entries;
}

// Discover EVIO files to process from a user-supplied path.
// Handles: single file, base name with splits (.00000, .00001, ...),
// or a directory (all .evio files sorted).
static std::vector<std::string> discover_files(const std::string &path)
{
    std::vector<std::string> files;

    if (is_regular_file(path)) {
        // single file — also look for sibling splits
        // e.g. "prad_023527.evio.00000" -> find all "prad_023527.evio.*"
        auto dot = path.rfind('.');
        if (dot != std::string::npos) {
            std::string suffix = path.substr(dot + 1);
            // check if suffix is all digits (split number)
            bool is_split = !suffix.empty() &&
                suffix.find_first_not_of("0123456789") == std::string::npos;
            if (is_split) {
                std::string base = path.substr(0, dot); // "path/prad_023527.evio"
                auto slash = base.rfind('/');
                if (slash == std::string::npos) slash = base.rfind('\\');
                std::string dir = (slash != std::string::npos) ? base.substr(0, slash) : ".";
                std::string base_name = (slash != std::string::npos) ? base.substr(slash + 1) : base;

                for (auto &f : list_dir(dir)) {
                    auto fslash = f.rfind('/');
                    if (fslash == std::string::npos) fslash = f.rfind('\\');
                    std::string fname = (fslash != std::string::npos) ? f.substr(fslash + 1) : f;
                    if (fname.size() > base_name.size() && fname.substr(0, base_name.size() + 1) == base_name + ".")
                        files.push_back(f);
                }
                std::sort(files.begin(), files.end());
                if (!files.empty()) return files;
            }
        }
        files.push_back(path);
        return files;
    }

    // path might be a base name without split suffix: "prad_023527.evio"
    if (!is_directory(path)) {
        auto slash = path.rfind('/');
        if (slash == std::string::npos) slash = path.rfind('\\');
        std::string dir = (slash != std::string::npos) ? path.substr(0, slash) : ".";
        std::string base_name = (slash != std::string::npos) ? path.substr(slash + 1) : path;

        for (auto &f : list_dir(dir)) {
            auto fslash = f.rfind('/');
            if (fslash == std::string::npos) fslash = f.rfind('\\');
            std::string fname = (fslash != std::string::npos) ? f.substr(fslash + 1) : f;
            if (fname.size() > base_name.size() && fname.substr(0, base_name.size() + 1) == base_name + ".")
                files.push_back(f);
        }
        if (files.empty()) {
            // try as exact file
            if (is_regular_file(path)) files.push_back(path);
        }
        std::sort(files.begin(), files.end());
        return files;
    }

    // directory — collect all .evio files
    for (auto &f : list_dir(path)) {
        if (f.find(".evio") != std::string::npos)
            files.push_back(f);
    }
    std::sort(files.begin(), files.end());
    return files;
}

// ---- output helpers ----------------------------------------------------------

static const char *tag_name(uint32_t tag)
{
    switch (tag) {
    case 0x00A9: return "SSP RawSum";
    case 0x00B0: return "Pulser 100Hz";
    case 0x00B9: return "LMS";
    case 0x00BA: return "Alpha";
    case 0x00BC: return "Master OR";
    case 0x00FA: return "Cluster";
    default:     return nullptr;
    }
}

// ---- main --------------------------------------------------------------------

int main(int argc, char *argv[])
{
    std::string input, daq_config_file;
    double pulser_freq = 100.0;   // Hz
    double report_interval = 10.0; // seconds

    std::string db_dir = DATABASE_DIR;
    if (const char *env = std::getenv("PRAD2_DATABASE_DIR")) db_dir = env;
    daq_config_file = db_dir + "/daq_config.json";

    int opt;
    while ((opt = getopt(argc, argv, "D:f:t:h")) != -1) {
        switch (opt) {
        case 'D': daq_config_file = optarg; break;
        case 'f': pulser_freq = std::atof(optarg); break;
        case 't': report_interval = std::atof(optarg); break;
        default:
            std::cerr << "Usage: " << argv[0]
                      << " <input> [-D daq_config.json] [-f freq_hz] [-t interval_sec]\n"
                      << "  <input>: file, base name (finds .00000 .00001 ...), or directory\n";
            return opt == 'h' ? 0 : 1;
        }
    }
    if (optind < argc) input = argv[optind];
    if (input.empty()) {
        std::cerr << "Usage: " << argv[0]
                  << " <input> [-D daq_config.json] [-f freq_hz] [-t interval_sec]\n";
        return 1;
    }

    DaqConfig cfg;
    if (!load_daq_config(daq_config_file, cfg)) {
        std::cerr << "Failed to load DAQ config: " << daq_config_file << "\n";
        return 1;
    }

    auto files = discover_files(input);
    if (files.empty()) {
        std::cerr << "No EVIO files found for: " << input << "\n";
        return 1;
    }
    std::cerr << "Processing " << files.size() << " file(s):\n";
    for (auto &f : files) std::cerr << "  " << f << "\n";

    // ---- state ----
    static constexpr double TI_TICK_SEC = 4e-9;
    static constexpr uint32_t PULSER_TAG = 0x00B0;
    static constexpr uint32_t DSC2_BANK_TAG = 0xE115;

    auto event = std::make_unique<fdec::EventData>();
    std::map<uint32_t, uint64_t> tag_counts;
    uint64_t first_ts = 0, last_ts = 0;
    uint64_t total_physics = 0;
    uint32_t run_number = 0;
    uint32_t unix_start = 0, unix_end = 0;
    int sync_count = 0;

    // DSC2 cumulative (latest SYNC)
    std::vector<Dsc2Slot> dsc2_data;

    // periodic reporting state
    double next_report = report_interval;

    // print periodic report header
    auto print_header = []() {
        std::cout << std::left
                  << std::setw(10) << "Time(s)"
                  << std::setw(12) << "LT_DSC2(%)"
                  << std::setw(14) << "LT_Pulser(%)"
                  << std::setw(10) << "Physics"
                  << std::setw(10) << "Pulser"
                  << std::setw(8)  << "Syncs"
                  << std::setw(0)  << "File"
                  << "\n" << std::string(78, '-') << "\n";
    };

    // emit one periodic row
    std::string current_file;
    auto print_row = [&](double elapsed) {
        // cumulative DSC2 live time (from latest SYNC)
        double lt_dsc2 = -1;
        if (!dsc2_data.empty() && dsc2_data[0].ref_ungated > 0) {
            lt_dsc2 = (1.0 - static_cast<double>(dsc2_data[0].ref_gated)
                           / dsc2_data[0].ref_ungated) * 100.0;
        }
        // cumulative pulser live time
        double lt_pulser = -1;
        uint64_t pulser_total = 0;
        auto it = tag_counts.find(PULSER_TAG);
        if (it != tag_counts.end()) pulser_total = it->second;
        double expected = pulser_freq * elapsed;
        if (expected > 0)
            lt_pulser = 100.0 * static_cast<double>(pulser_total) / expected;

        // extract just the filename from path
        auto slash = current_file.rfind('/');
        if (slash == std::string::npos) slash = current_file.rfind('\\');
        std::string fname = (slash != std::string::npos) ? current_file.substr(slash + 1) : current_file;

        std::cout << std::fixed
                  << std::setw(10) << std::setprecision(1) << elapsed;
        if (lt_dsc2 >= 0)
            std::cout << std::setw(12) << std::setprecision(2) << lt_dsc2;
        else
            std::cout << std::setw(12) << "--";
        if (lt_pulser >= 0)
            std::cout << std::setw(14) << std::setprecision(2) << lt_pulser;
        else
            std::cout << std::setw(14) << "--";
        std::cout << std::setw(10) << total_physics
                  << std::setw(10) << pulser_total
                  << std::setw(8)  << sync_count
                  << fname
                  << "\n";
    };

    print_header();

    // ---- process files ----
    EvChannel ch;
    ch.SetConfig(cfg);

    for (auto &file : files) {
        current_file = file;
        if (ch.Open(file) != status::success) {
            std::cerr << "Warning: cannot open " << file << ", skipping\n";
            continue;
        }

        while (ch.Read() == status::success) {
            if (!ch.Scan()) continue;
            auto evtype = ch.GetEventType();

            if (evtype == EventType::Prestart || evtype == EventType::Go) {
                uint32_t ct = ch.Sync().unix_time;
                if (ct != 0 && unix_start == 0) unix_start = ct;
            }
            if (evtype == EventType::End) {
                uint32_t ct = ch.Sync().unix_time;
                if (ct != 0) unix_end = ct;
            }

            // DSC2 scalers
            if (evtype == EventType::Sync || evtype == EventType::Physics) {
                const EvNode *dsc2_node = ch.FindFirstByTag(DSC2_BANK_TAG);
                if (dsc2_node && dsc2_node->data_words > 0) {
                    auto parsed = parse_dsc2_bank(ch.GetData(*dsc2_node), dsc2_node->data_words);
                    if (!parsed.empty()) {
                        dsc2_data = std::move(parsed);
                        if (evtype == EventType::Sync) sync_count++;
                    }
                }
            }

            if (evtype != EventType::Physics) continue;

            for (int ie = 0; ie < ch.GetNEvents(); ++ie) {
                event->clear();
                if (!ch.DecodeEvent(ie, *event)) continue;

                tag_counts[event->info.event_tag]++;
                total_physics++;

                if (run_number == 0 && event->info.run_number != 0)
                    run_number = event->info.run_number;

                uint64_t ts = event->info.timestamp;
                if (ts != 0) {
                    if (first_ts == 0) first_ts = ts;
                    last_ts = ts;

                    // periodic report
                    double elapsed = static_cast<double>(ts - first_ts) * TI_TICK_SEC;
                    while (elapsed >= next_report) {
                        print_row(next_report);
                        next_report += report_interval;
                    }
                }
            }
        }
        ch.Close();
    }

    // ---- final report ----
    double elapsed_ti = (first_ts != 0 && last_ts > first_ts)
        ? static_cast<double>(last_ts - first_ts) * TI_TICK_SEC : 0.0;
    double elapsed_unix = (unix_start != 0 && unix_end > unix_start)
        ? static_cast<double>(unix_end - unix_start) : 0.0;
    double elapsed = (elapsed_ti > 0) ? elapsed_ti : elapsed_unix;

    // final row (partial interval)
    if (elapsed > 0)
        print_row(elapsed);

    // ---- summary ----
    std::cout << "\n=== Summary ===\n";
    if (run_number != 0)
        std::cout << "Run number     : " << run_number << "\n";
    std::cout << "Files          : " << files.size() << "\n";
    std::cout << "Total physics  : " << total_physics << "\n";
    std::cout << "SYNC events    : " << sync_count << "\n";
    if (elapsed_ti > 0)
        std::cout << std::fixed << std::setprecision(2)
                  << "Elapsed (TI)   : " << elapsed_ti << " sec\n";
    if (elapsed_unix > 0)
        std::cout << std::fixed << std::setprecision(0)
                  << "Elapsed (unix) : " << elapsed_unix << " sec\n";

    std::cout << "\n--- Trigger counts ---\n";
    std::cout << std::left << std::setw(10) << "Tag"
              << std::setw(20) << "Name"
              << std::right << std::setw(10) << "Count" << "\n";
    std::cout << std::string(40, '-') << "\n";
    for (auto &[tag, count] : tag_counts) {
        const char *name = tag_name(tag);
        char hex[16];
        snprintf(hex, sizeof(hex), "0x%04X", tag);
        std::cout << std::left << std::setw(10) << hex
                  << std::setw(20) << (name ? name : "unknown")
                  << std::right << std::setw(10) << count << "\n";
    }

    // DSC2 summary
    std::cout << "\n--- DSC2 scaler live time (cumulative) ---\n";
    if (dsc2_data.empty()) {
        std::cout << "(no DSC2 scaler bank 0xE115 found)\n";
    } else {
        static constexpr double DSC2_REF_FREQ = 125.0e6;
        for (auto &s : dsc2_data) {
            double lt = (s.ref_ungated > 0)
                ? (1.0 - static_cast<double>(s.ref_gated) / s.ref_ungated) * 100.0
                : 0.0;
            std::cout << std::fixed << std::setprecision(2);
            std::cout << "  DSC2 slot " << s.slot
                      << ": ref_gated=" << s.ref_gated
                      << "  ref_ungated=" << s.ref_ungated
                      << "  live=" << std::setprecision(3) << lt << "%\n";

            bool any = false;
            for (int c = 0; c < DSC2_NCH; ++c) {
                if (s.trg_ungated[c] == 0) continue;
                if (!any) {
                    std::cout << "    TRG ch  gated(busy)  ungated(total)  live%\n";
                    any = true;
                }
                double cl = (1.0 - static_cast<double>(s.trg_gated[c])
                                 / s.trg_ungated[c]) * 100.0;
                std::cout << "      " << std::setw(2) << c
                          << std::setw(13) << s.trg_gated[c]
                          << std::setw(16) << s.trg_ungated[c]
                          << std::setw(9) << std::setprecision(2) << cl << "\n";
            }
        }
    }

    // pulser summary
    uint64_t pulser_count = 0;
    auto pit = tag_counts.find(PULSER_TAG);
    if (pit != tag_counts.end()) pulser_count = pit->second;
    double expected_pulser = pulser_freq * elapsed;
    double lt_pulser = (expected_pulser > 0)
        ? 100.0 * static_cast<double>(pulser_count) / expected_pulser : 0.0;

    std::cout << "\n--- Pulser counting live time ---\n";
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "  Pulser freq    : " << pulser_freq << " Hz\n";
    std::cout << "  Pulser accepted: " << pulser_count << "\n";
    if (expected_pulser > 0)
        std::cout << "  Pulser expected: " << std::setprecision(0) << expected_pulser << "\n";
    std::cout << std::setprecision(2);
    if (lt_pulser > 0)
        std::cout << "  Live time      : " << lt_pulser << " %\n";
    else
        std::cout << "  Live time      : N/A\n";

    return 0;
}
