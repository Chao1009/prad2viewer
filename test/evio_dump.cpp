// test/evio_dump.cpp
// Diagnostic tool to inspect EVIO file structure.
//
// Modes:
//   evio_dump <file>                          -- summary: count events by type/tag
//   evio_dump <file> --tree [--num N]         -- print bank tree for first N events
//   evio_dump <file> --tags                   -- list all unique bank tags with counts
//   evio_dump <file> --epics                  -- dump EPICS text from all EPICS events
//   evio_dump <file> --event N                -- detailed dump of event N (1-based)

#include "EvChannel.h"
#include "EvStruct.h"
#include "Fadc250Data.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <map>
#include <set>
#include <vector>
#include <bitset>
#include <cstdlib>

using namespace evc;

// --- helpers ----------------------------------------------------------------
static std::string hex(uint32_t v)
{
    char buf[16];
    snprintf(buf, sizeof(buf), "0x%04X", v);
    return buf;
}

static std::string tag_label(uint32_t tag)
{
    // CODA3 physics events (page 27/29)
    if (tag == 0xFF50) return "PHYSICS(PEB)";
    if (tag == 0xFF58) return "PHYSICS(PEB+sync)";
    if (tag == 0xFF70) return "PHYSICS(SEB)";
    if (tag == 0xFF78) return "PHYSICS(SEB+sync)";
    if (tag >= 0xFF50 && tag <= 0xFF8F) return "PHYSICS";

    // CODA3 trigger banks (page 26)
    if (tag >= 0xFF10 && tag <= 0xFF1F) return "TRIGGER(raw)";
    if (tag >= 0xFF20 && tag <= 0xFF2F) return "TRIGGER(built)";
    if (tag == 0xFF4F) return "TRIGGER(bad)";

    // CODA3 control events (page 20/29)
    if (tag == 0xFFD0) return "SYNC";
    if (tag == 0xFFD1) return "PRESTART";
    if (tag == 0xFFD2) return "GO";
    if (tag == 0xFFD3) return "PAUSE";
    if (tag == 0xFFD4) return "END";
    if (tag >= 0xFFD0 && tag <= 0xFFDF) return "CONTROL";

    // Legacy CODA2 tags (may appear in older data)
    if (tag == 0x11) return "PRESTART(legacy)";
    if (tag == 0x12) return "GO(legacy)";
    if (tag == 0x14) return "END(legacy)";
    if (tag == 0xC1) return "SYNC(legacy)";

    // JLab single-event physics
    if (tag == 0xB1) return "PHYSICS(single)";
    if (tag == 0xFE) return "PHYSICS(single)";

    // JLab-specific banks
    if (tag == 0xC000) return "TRIGGER_BANK";

    // EPICS
    if (tag == 0x1F) return "EPICS";

    return "";
}

// --- mode: summary ----------------------------------------------------------
static int doSummary(EvChannel &ch)
{
    std::map<uint32_t, int> tag_counts;
    int total = 0;

    while (ch.Read() == status::success) {
        auto hdr = ch.GetEvHeader();
        tag_counts[hdr.tag]++;
        total++;
    }

    std::cout << "=== Event Summary ===\n";
    std::cout << "Total EVIO records: " << total << "\n\n";
    std::cout << std::setw(12) << "Tag" << std::setw(12) << "Count"
              << "  " << "Label" << "\n";
    std::cout << std::string(40, '-') << "\n";

    for (auto &[tag, cnt] : tag_counts) {
        std::cout << std::setw(12) << hex(tag) << std::setw(12) << cnt
                  << "  " << tag_label(tag) << "\n";
    }
    return 0;
}

// --- mode: tree -------------------------------------------------------------
static int doTree(EvChannel &ch, int num)
{
    int count = 0;
    while (ch.Read() == status::success) {
        if (!ch.Scan()) continue;
        if (++count > num) break;

        auto hdr = ch.GetEvHeader();
        std::cout << "========== Record " << count
                  << "  tag=" << hex(hdr.tag)
                  << " (" << tag_label(hdr.tag) << ")"
                  << "  num=" << hdr.num
                  << "  len=" << hdr.length << "w"
                  << " ==========\n";
        ch.PrintTree(std::cout);
        std::cout << "\n";
    }
    std::cout << "Printed " << std::min(count, num) << " record(s).\n";
    return 0;
}

// --- mode: tags (deep scan of all bank tags) --------------------------------
struct TagInfo {
    uint32_t type;
    int      count;
    int      min_depth, max_depth;
    size_t   min_words, max_words;
    std::set<uint32_t> parent_tags;
};

static int doTags(EvChannel &ch)
{
    std::map<uint32_t, TagInfo> all_tags;
    int nrecords = 0;

    while (ch.Read() == status::success) {
        if (!ch.Scan()) continue;
        nrecords++;

        for (auto &n : ch.GetNodes()) {
            auto it = all_tags.find(n.tag);
            if (it == all_tags.end()) {
                TagInfo ti;
                ti.type = n.type;
                ti.count = 1;
                ti.min_depth = ti.max_depth = n.depth;
                ti.min_words = ti.max_words = n.data_words;
                if (n.parent >= 0)
                    ti.parent_tags.insert(ch.GetNodes()[n.parent].tag);
                all_tags[n.tag] = ti;
            } else {
                it->second.count++;
                it->second.min_depth = std::min(it->second.min_depth, n.depth);
                it->second.max_depth = std::max(it->second.max_depth, n.depth);
                it->second.min_words = std::min(it->second.min_words, n.data_words);
                it->second.max_words = std::max(it->second.max_words, n.data_words);
                if (n.parent >= 0)
                    it->second.parent_tags.insert(ch.GetNodes()[n.parent].tag);
            }
        }
    }

    std::cout << "=== All Bank Tags (across " << nrecords << " records) ===\n\n";
    std::cout << std::setw(12) << "Tag"
              << std::setw(10) << "Type"
              << std::setw(8) << "Count"
              << std::setw(8) << "Depth"
              << std::setw(16) << "Data words"
              << "  Parents\n";
    std::cout << std::string(80, '-') << "\n";

    for (auto &[tag, ti] : all_tags) {
        std::cout << std::setw(12) << hex(tag)
                  << std::setw(10) << TypeName(ti.type)
                  << std::setw(8) << ti.count;

        if (ti.min_depth == ti.max_depth)
            std::cout << std::setw(8) << ti.min_depth;
        else
            std::cout << std::setw(3) << ti.min_depth << "-" << std::setw(3) << ti.max_depth << " ";

        if (ti.min_words == ti.max_words)
            std::cout << std::setw(16) << ti.min_words;
        else
            std::cout << std::setw(7) << ti.min_words << "-" << std::setw(7) << ti.max_words << "  ";

        std::cout << "  ";
        for (auto pt : ti.parent_tags)
            std::cout << hex(pt) << " ";
        std::cout << "\n";
    }
    return 0;
}

// --- mode: epics ------------------------------------------------------------
static int doEpics(EvChannel &ch)
{
    int count = 0, record = 0;

    while (ch.Read() == status::success) {
        record++;
        if (!ch.Scan()) continue;

        auto hdr = ch.GetEvHeader();
        // check for EPICS: common tags are 0x1F, but also scan for string banks
        bool is_epics = (hdr.tag == 0x1F || hdr.tag == 0x1f);

        if (!is_epics) {
            // also check if any child bank has string data
            for (auto &n : ch.GetNodes()) {
                if (n.depth == 1 &&
                    (n.type == DATA_CHARSTAR8 || n.type == DATA_CHAR8) &&
                    n.data_words > 4)
                {
                    is_epics = true;
                    break;
                }
            }
        }

        if (!is_epics) continue;

        count++;
        std::cout << "--- EPICS record " << count
                  << " (file record " << record
                  << ", tag=" << hex(hdr.tag)
                  << ", num=" << hdr.num << ") ---\n";

        // print tree structure
        ch.PrintTree(std::cout);

        // extract and print all string data
        for (auto &n : ch.GetNodes()) {
            if ((n.type == DATA_CHARSTAR8 || n.type == DATA_CHAR8) && n.data_words > 0) {
                const char *raw = reinterpret_cast<const char*>(ch.GetData(n));
                size_t max_len = n.data_words * 4;
                size_t len = 0;
                while (len < max_len && raw[len] != '\0') ++len;

                std::cout << "\n  [String data from tag=" << hex(n.tag)
                          << ", " << len << " bytes]:\n";

                // print first 2000 chars
                size_t show = std::min(len, size_t(2000));
                std::string text(raw, show);
                std::cout << text;
                if (show < len) std::cout << "\n  ... (" << len - show << " more bytes)";
                std::cout << "\n";
            }
        }
        std::cout << "\n";
    }

    std::cout << "Found " << count << " EPICS record(s) in " << record << " total records.\n";
    return 0;
}

// --- mode: single event detail ----------------------------------------------
static int doEvent(EvChannel &ch, int target)
{
    int record = 0;

    while (ch.Read() == status::success) {
        record++;
        if (!ch.Scan()) continue;

        if (record != target) continue;

        auto hdr = ch.GetEvHeader();
        std::cout << "=== Record " << record
                  << "  tag=" << hex(hdr.tag)
                  << " (" << tag_label(hdr.tag) << ")"
                  << "  type=0x" << std::hex << hdr.type << std::dec
                  << "  num=" << hdr.num
                  << "  length=" << hdr.length << "w"
                  << " ===\n\n";

        // full tree
        std::cout << "--- Bank Tree ---\n";
        ch.PrintTree(std::cout);

        // if physics, try decoding
        if (ch.GetNEvents() > 0) {
            std::cout << "\n--- Physics Decode ---\n";
            std::cout << "Sub-events in block: " << ch.GetNEvents() << "\n";

            fdec::EventData evt;
            for (int i = 0; i < ch.GetNEvents(); ++i) {
                if (!ch.DecodeEvent(i, evt)) {
                    std::cout << "  sub-event " << i << ": decode failed\n";
                    continue;
                }

                std::cout << "  sub-event " << i
                          << ": event#=" << evt.info.event_number
                          << " trigger#=" << evt.info.trigger_number
                          << " trigger_bits=0x" << std::hex
                          << evt.info.trigger_bits << std::dec
                          << " timestamp=" << evt.info.timestamp
                          << " run=" << evt.info.run_number
                          << " unix_time=" << evt.info.unix_time
                          << " rocs=" << evt.nrocs << "\n";

                for (int r = 0; r < evt.nrocs; ++r) {
                    auto &roc = evt.rocs[r];
                    if (!roc.present) continue;
                    std::cout << "    ROC tag=" << hex(roc.tag)
                              << " slots=" << roc.nslots << "\n";

                    for (int s = 0; s < fdec::MAX_SLOTS; ++s) {
                        auto &slot = roc.slots[s];
                        if (!slot.present) continue;

                        int nch = 0;
                        for (int c = 0; c < fdec::MAX_CHANNELS; ++c)
                            if (slot.channel_mask & (1u << c)) nch++;

                        std::cout << "      slot=" << std::setw(2) << s
                                  << " trigger=" << slot.trigger
                                  << " timestamp=" << slot.timestamp
                                  << " channels=" << nch << " [";

                        for (int c = 0; c < fdec::MAX_CHANNELS; ++c) {
                            if (!(slot.channel_mask & (1u << c))) continue;
                            std::cout << c << "(" << slot.channels[c].nsamples << ")";
                            if (c < fdec::MAX_CHANNELS - 1) std::cout << " ";
                        }
                        std::cout << "]\n";
                    }
                }
            }
        }

        // dump raw words of any non-container leaf banks
        std::cout << "\n--- Leaf Bank Data (first 16 words) ---\n";
        for (auto &n : ch.GetNodes()) {
            if (n.child_count > 0 || n.data_words == 0) continue;
            if (IsContainer(n.type) || n.type == DATA_COMPOSITE) continue;

            std::cout << "  tag=" << hex(n.tag)
                      << " type=" << TypeName(n.type)
                      << " depth=" << n.depth
                      << " words=" << n.data_words << " |";

            const uint32_t *d = ch.GetData(n);
            int show = std::min<int>(n.data_words, 16);
            for (int i = 0; i < show; ++i)
                std::cout << " " << std::hex << std::setw(8) << std::setfill('0')
                          << d[i] << std::setfill(' ') << std::dec;
            if (n.data_words > 16) std::cout << " ...";
            std::cout << "\n";
        }

        return 0;
    }

    std::cerr << "Record " << target << " not found (file has " << record << " records).\n";
    return 1;
}

// --- mode: triggers ---------------------------------------------------------
static int doTriggers(EvChannel &ch)
{
    fdec::EventData evt;
    int record = 0, decoded = 0;
    std::map<uint32_t, int> trig_counts;

    std::cout << std::setw(8) << "event#"
              << std::setw(10) << "trigger#"
              << std::setw(14) << "trigger_bits"
              << std::setw(18) << "timestamp"
              << std::setw(8) << "rocs"
              << "\n";
    std::cout << std::string(58, '-') << "\n";

    while (ch.Read() == status::success) {
        record++;
        if (!ch.Scan()) continue;
        if (ch.GetNEvents() == 0) continue;

        for (int i = 0; i < ch.GetNEvents(); ++i) {
            if (!ch.DecodeEvent(i, evt)) continue;
            decoded++;
            trig_counts[evt.info.trigger_bits]++;

            std::cout << std::setw(8) << evt.info.event_number
                      << std::setw(10) << evt.info.trigger_number
                      << "    0x" << std::hex << std::setw(8)
                      << std::setfill('0') << evt.info.trigger_bits
                      << std::dec << std::setfill(' ')
                      << std::setw(18) << evt.info.timestamp
                      << std::setw(8) << evt.nrocs
                      << "\n";
        }
    }

    std::cout << "\n=== Trigger Bits Summary (" << decoded << " events) ===\n";
    for (auto &[bits, cnt] : trig_counts) {
        std::cout << "  0x" << std::hex << std::setw(8) << std::setfill('0')
                  << bits << std::dec << std::setfill(' ')
                  << "  count=" << cnt << "\n";
    }

    return 0;
}

// --- main -------------------------------------------------------------------
static void usage(const char *prog)
{
    std::cerr
        << "EVIO file structure diagnostic tool\n\n"
        << "Usage:\n"
        << "  " << prog << " <file>                      Summary: count events by tag\n"
        << "  " << prog << " <file> --tree [--num N]      Print bank tree (default N=5)\n"
        << "  " << prog << " <file> --tags                List all unique bank tags with stats\n"
        << "  " << prog << " <file> --epics               Dump all EPICS event text\n"
        << "  " << prog << " <file> --event N             Detailed dump of record N (1-based)\n"
        << "  " << prog << " <file> --triggers            List trigger info for all events\n";
}

int main(int argc, char *argv[])
{
    if (argc < 2) { usage(argv[0]); return 1; }

    std::string path = argv[1];
    std::string mode = (argc > 2) ? argv[2] : "";

    EvChannel ch;
    if (ch.Open(path) != status::success) {
        std::cerr << "Failed to open: " << path << "\n";
        return 1;
    }

    int rc;
    if (mode == "--tree") {
        int num = 5;
        for (int i = 3; i < argc; ++i)
            if (std::string(argv[i]) == "--num" && i + 1 < argc) num = std::atoi(argv[++i]);
        rc = doTree(ch, num);
    }
    else if (mode == "--tags") {
        rc = doTags(ch);
    }
    else if (mode == "--epics") {
        rc = doEpics(ch);
    }
    else if (mode == "--triggers") {
        rc = doTriggers(ch);
    }
    else if (mode == "--event") {
        if (argc < 4) { usage(argv[0]); return 1; }
        rc = doEvent(ch, std::atoi(argv[3]));
    }
    else {
        rc = doSummary(ch);
    }

    ch.Close();
    return rc;
}
