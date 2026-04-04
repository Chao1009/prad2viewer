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
#include "load_daq_config.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <map>
#include <set>
#include <vector>
#include <getopt.h>
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
    if ((tag >= 0x00A0 && tag <= 0x00BF) || (tag >= 0xFF50 && tag <= 0xFF8F)) return "PHYSICS(built)";

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
                            if (slot.channel_mask & (1ull << c)) nch++;

                        std::cout << "      slot=" << std::setw(2) << s
                                  << " trigger=" << slot.trigger
                                  << " timestamp=" << slot.timestamp
                                  << " channels=" << nch << " [";

                        for (int c = 0; c < fdec::MAX_CHANNELS; ++c) {
                            if (!(slot.channel_mask & (1ull << c))) continue;
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

// --- mode: trig-debug -------------------------------------------------------
// Cross-correlates three trigger layers:
//   1. Event tag (top-level) = 0x80 + TI_event_type
//   2. TI event_type (d[0] bits 31:24) = TS trigger decision
//   3. FP trigger bits (TI master d[5]) = raw 32-bit front panel input snapshot
//
// FP bit assignments (from prad_v0.trg):
//   Bits 8-15:  SSP PRAD TRGBIT 0-7  (P2 outputs)
//   Bit 23:     v1495 OR from SD/FADC
//   Bit 24:     v1495 LMS
//   Bit 25:     v1495 alpha
//   Bit 26:     v1495 Faraday
//   Bit 27:     v1495 Master OR
//
// TS monitors mask: 0x0F00FF00 (bits 8-15 and 24-27)

struct TrigDebugEntry {
    uint32_t event_tag;         // top-level bank tag
    uint32_t ti_event_type;     // d[0] >> 24 (from any TI bank)
    uint32_t ti_nwords;         // nwords from TI event header
    uint32_t fp_trigger_bits;   // d[5] from TI master (0 if unavailable)
    int32_t  event_number;      // from 0xC000 or TI d[1]
    int      nrocs;             // number of FADC composite banks found
    bool     has_ti_master;     // TI master bank found
    bool     tag_matches;       // event_tag == 0x80 + ti_event_type
};

static const uint32_t TS_FP_MASK = 0x0F00FF00;

struct FpBitInfo {
    int bit;
    uint32_t mask;
    const char *name;
};

static const FpBitInfo fp_bits[] = {
    {  8, 0x00000100, "SSP TRGBIT0 (RawSum>1000)" },
    {  9, 0x00000200, "SSP TRGBIT1 (1clus>1GeV)"  },
    { 10, 0x00000400, "SSP TRGBIT2 (2clus>1GeV)"  },
    { 11, 0x00000800, "SSP TRGBIT3 (3clus>1GeV)"  },
    { 12, 0x00001000, "SSP TRGBIT4 (disabled)"     },
    { 13, 0x00002000, "SSP TRGBIT5 (disabled)"     },
    { 14, 0x00004000, "SSP TRGBIT6 (disabled)"     },
    { 15, 0x00008000, "SSP TRGBIT7 (100Hz pulser)" },
    { 23, 0x00800000, "v1495 OR from SD/FADC"      },
    { 24, 0x01000000, "v1495 LMS"                   },
    { 25, 0x02000000, "v1495 alpha"                 },
    { 26, 0x04000000, "v1495 Faraday"               },
    { 27, 0x08000000, "v1495 Master OR"             },
};
static const int N_FP_BITS = sizeof(fp_bits) / sizeof(fp_bits[0]);

static int doTrigDebug(EvChannel &ch, bool verbose)
{
    const auto &cfg = ch.GetConfig();
    std::vector<TrigDebugEntry> entries;

    // per-tag accumulators
    struct TagStats {
        int count = 0;
        int has_fadc = 0;
        uint32_t fp_or = 0;         // OR of all FP bits seen with this tag
        uint32_t fp_and = 0xFFFFFFFF; // AND of all FP bits
        int fp_bit_counts[32] = {};
    };
    std::map<uint32_t, TagStats> tag_stats;

    int record = 0, physics = 0;

    while (ch.Read() == status::success) {
        record++;
        if (!ch.Scan()) continue;
        if (ch.GetNEvents() == 0) continue;

        auto hdr = ch.GetEvHeader();
        auto &nodes = ch.GetNodes();

        for (int iev = 0; iev < ch.GetNEvents(); ++iev) {
            TrigDebugEntry e = {};
            e.event_tag = hdr.tag;

            // --- extract from 0xC000 trigger bank ---
            for (auto &n : nodes) {
                if (n.tag == cfg.trigger_bank_tag && n.type == DATA_UINT32 && n.data_words >= 1) {
                    e.event_number = static_cast<int32_t>(ch.GetData(n)[0]);
                    break;
                }
            }

            // --- extract TI event_type from first 0xE10A bank ---
            for (auto &n : nodes) {
                if (n.tag == cfg.ti_bank_tag && n.type == DATA_UINT32 && n.data_words >= 1) {
                    const uint32_t *d = ch.GetData(n);
                    e.ti_event_type = (d[0] >> 24) & 0xFF;
                    e.ti_nwords = d[0] & 0xFFFF;
                    if (n.data_words >= 2 && e.event_number == 0)
                        e.event_number = static_cast<int32_t>(d[1]);
                    break;
                }
            }

            // --- extract FP trigger bits from TI master's 0xE10A ---
            for (auto &n : nodes) {
                if (n.depth == 1 && n.tag == cfg.ti_master_tag) {
                    e.has_ti_master = true;
                    for (size_t ci = 0; ci < n.child_count; ++ci) {
                        auto &child = nodes[n.child_first + ci];
                        if (child.tag == cfg.ti_bank_tag && child.type == DATA_UINT32) {
                            const uint32_t *d = ch.GetData(child);
                            size_t nw = child.data_words;
                            // d[4] = 8-bit trigger type byte
                            // d[5] = 32-bit FP trigger inputs (if FP readout enabled)
                            if (nw > 5)
                                e.fp_trigger_bits = d[5];
                            break;
                        }
                    }
                    break;
                }
            }

            // --- count FADC composite banks ---
            for (auto &n : nodes) {
                if (n.tag == cfg.fadc_composite_tag && n.type == DATA_COMPOSITE)
                    e.nrocs++;
            }

            // --- verify tag = 0x80 + event_type ---
            e.tag_matches = (e.event_tag == (0x80u + e.ti_event_type));

            // --- accumulate ---
            physics++;
            entries.push_back(e);

            auto &ts = tag_stats[e.event_tag];
            ts.count++;
            if (e.nrocs > 0) ts.has_fadc++;
            ts.fp_or |= e.fp_trigger_bits;
            ts.fp_and &= e.fp_trigger_bits;
            for (int b = 0; b < 32; ++b)
                if (e.fp_trigger_bits & (1u << b)) ts.fp_bit_counts[b]++;
        }
    }

    // === Per-event detail (-v) ===
    if (verbose) {
        std::cout << std::setw(8) << "event#"
                  << std::setw(8) << "tag"
                  << std::setw(8) << "TItype"
                  << std::setw(6) << "chk"
                  << std::setw(12) << "FP_bits"
                  << std::setw(12) << "FP_masked"
                  << std::setw(6) << "FADC"
                  << "  active FP signals"
                  << "\n";
        std::cout << std::string(90, '-') << "\n";

        for (auto &e : entries) {
            uint32_t masked = e.fp_trigger_bits & TS_FP_MASK;
            std::cout << std::setw(8) << e.event_number
                      << "  0x" << std::hex << std::setw(4) << std::setfill('0') << e.event_tag
                      << "    0x" << std::setw(2) << e.ti_event_type
                      << std::setfill(' ') << std::dec
                      << std::setw(6) << (e.tag_matches ? "OK" : "FAIL")
                      << "  0x" << std::hex << std::setw(8) << std::setfill('0') << e.fp_trigger_bits
                      << "  0x" << std::setw(8) << masked
                      << std::setfill(' ') << std::dec
                      << std::setw(6) << e.nrocs;

            // list active FP bit names
            std::cout << "  ";
            bool first = true;
            for (int i = 0; i < N_FP_BITS; ++i) {
                if (e.fp_trigger_bits & fp_bits[i].mask) {
                    if (!first) std::cout << ", ";
                    std::cout << "b" << fp_bits[i].bit;
                    first = false;
                }
            }
            std::cout << "\n";
        }
        std::cout << "\n";
    }

    // === Per-tag summary ===
    std::cout << "=== Trigger Debug Summary (" << physics << " physics events, "
              << record << " records) ===\n\n";

    std::cout << "--- Event Tag Summary ---\n";
    std::cout << std::setw(8) << "tag"
              << std::setw(8) << "TItype"
              << std::setw(8) << "count"
              << std::setw(8) << "w/FADC"
              << std::setw(12) << "FP_OR"
              << std::setw(12) << "FP_AND"
              << "\n";
    std::cout << std::string(56, '-') << "\n";

    int mismatch_total = 0;
    for (auto &[tag, ts] : tag_stats) {
        uint32_t expected_type = tag - 0x80;
        std::cout << "  0x" << std::hex << std::setw(4) << std::setfill('0') << tag
                  << "    0x" << std::setw(2) << expected_type
                  << std::setfill(' ') << std::dec
                  << std::setw(8) << ts.count
                  << std::setw(8) << ts.has_fadc
                  << "  0x" << std::hex << std::setw(8) << std::setfill('0') << ts.fp_or
                  << "  0x" << std::setw(8) << ts.fp_and
                  << std::setfill(' ') << std::dec
                  << "\n";
    }

    // === Tag verification ===
    for (auto &e : entries)
        if (!e.tag_matches) mismatch_total++;

    std::cout << "\n--- Tag Verification ---\n";
    std::cout << "  tag == 0x80 + TI_event_type: "
              << (physics - mismatch_total) << " OK, "
              << mismatch_total << " MISMATCH\n";
    if (mismatch_total > 0) {
        std::cout << "  First mismatches:\n";
        int shown = 0;
        for (auto &e : entries) {
            if (!e.tag_matches && shown < 10) {
                std::cout << "    event#=" << e.event_number
                          << " tag=0x" << std::hex << e.event_tag
                          << " TI_type=0x" << e.ti_event_type
                          << " expected_tag=0x" << (0x80 + e.ti_event_type)
                          << std::dec << "\n";
                shown++;
            }
        }
    }

    // === FP bit activity per tag ===
    std::cout << "\n--- FP Bit Activity (TS mask = 0x"
              << std::hex << std::setw(8) << std::setfill('0') << TS_FP_MASK
              << std::setfill(' ') << std::dec << ") ---\n";

    // header
    std::cout << std::setw(5) << "bit" << std::setw(12) << "hex"
              << "  " << std::setw(30) << std::left << "name" << std::right
              << std::setw(8) << "total";
    for (auto &[tag, ts] : tag_stats)
        std::cout << std::setw(8) << ("0x" + hex(tag).substr(2));
    std::cout << "\n";
    std::cout << std::string(55 + 8 * tag_stats.size(), '-') << "\n";

    for (int i = 0; i < N_FP_BITS; ++i) {
        int bit = fp_bits[i].bit;
        bool in_mask = (TS_FP_MASK & fp_bits[i].mask) != 0;

        // count total across all tags
        int total = 0;
        for (auto &[tag, ts] : tag_stats)
            total += ts.fp_bit_counts[bit];

        if (total == 0 && !in_mask) continue; // skip unused bits not in mask

        std::cout << std::setw(5) << bit
                  << "  0x" << std::hex << std::setw(8) << std::setfill('0')
                  << fp_bits[i].mask << std::setfill(' ') << std::dec
                  << (in_mask ? "* " : "  ")
                  << std::setw(30) << std::left << fp_bits[i].name << std::right
                  << std::setw(8) << total;

        for (auto &[tag, ts] : tag_stats)
            std::cout << std::setw(8) << ts.fp_bit_counts[bit];
        std::cout << "\n";
    }
    std::cout << "\n  * = in TS_FP_INPUT_MASK\n";

    // === Check d[4] vs d[0] consistency ===
    std::cout << "\n--- TI Event Header vs d[4] ---\n";
    int d4_match = 0, d4_mismatch = 0, d4_unavail = 0;
    for (auto &e : entries) {
        if (e.ti_nwords <= 3) { d4_unavail++; continue; }
        // We don't have d[4] stored; but we can verify TI type from d[0]
        // matches event_tag. Already done above.
    }
    std::cout << "  (d[4] check requires direct bank access — use -m event -n N for individual inspection)\n";
    std::cout << "  TI event header nwords distribution:\n";
    std::map<uint32_t, int> nwords_dist;
    for (auto &e : entries) nwords_dist[e.ti_nwords]++;
    for (auto &[nw, cnt] : nwords_dist)
        std::cout << "    nwords=" << nw << ": " << cnt << " events\n";

    return 0;
}

// --- mode: triggers ---------------------------------------------------------
static int doTriggers(EvChannel &ch, bool verbose)
{
    fdec::EventData evt;
    int record = 0, decoded = 0;
    std::map<uint32_t, int> trig_counts;

    if (verbose) {
        std::cout << std::setw(8) << "event#"
                  << std::setw(10) << "trigger#"
                  << std::setw(14) << "trigger_bits"
                  << std::setw(18) << "timestamp"
                  << std::setw(8) << "rocs"
                  << "\n";
        std::cout << std::string(58, '-') << "\n";
    }

    while (ch.Read() == status::success) {
        record++;
        if (!ch.Scan()) continue;
        if (ch.GetNEvents() == 0) continue;

        for (int i = 0; i < ch.GetNEvents(); ++i) {
            ch.DecodeEvent(i, evt);
            decoded++;
            trig_counts[evt.info.trigger_bits]++;

            if (verbose) {
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
    }

    std::cout << "=== Trigger Bits Summary (" << decoded << " events) ===\n";
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
        << "  " << prog << " <file> [options]\n\n"
        << "Modes (default: summary):\n"
        << "  -m tree       Print bank tree\n"
        << "  -m tags       List all unique bank tags with stats\n"
        << "  -m epics      Dump all EPICS event text\n"
        << "  -m event      Detailed dump of a single record\n"
        << "  -m triggers   List trigger bit counts (add -v for per-event detail)\n"
        << "  -m trig-debug Cross-correlate event tags, TI event type, and FP trigger bits\n\n"
        << "Options:\n"
        << "  -D <file>     DAQ configuration (auto-searches daq_config.json if omitted)\n"
        << "  -n <N>        Number of events (tree mode, default 5) or event number (event mode)\n"
        << "  -v            Verbose output (triggers mode: print every event)\n";
}

int main(int argc, char *argv[])
{
    std::string daq_config_file;
    std::string mode;
    int num = 5;
    bool verbose = false;

    int opt;
    while ((opt = getopt(argc, argv, "D:m:n:vh")) != -1) {
        switch (opt) {
        case 'D': daq_config_file = optarg; break;
        case 'm': mode = optarg; break;
        case 'n': num = std::atoi(optarg); break;
        case 'v': verbose = true; break;
        default:  usage(argv[0]); return 1;
        }
    }
    if (optind >= argc) { usage(argv[0]); return 1; }
    std::string path = argv[optind];

    // auto-search for daq_config.json if not specified
    if (daq_config_file.empty()) {
        for (auto p : {"daq_config.json", "database/daq_config.json", "../database/daq_config.json"}) {
            std::ifstream f(p);
            if (f.good()) { daq_config_file = p; break; }
        }
    }

    evc::DaqConfig daq_cfg;
    if (daq_config_file.empty() || !evc::load_daq_config(daq_config_file, daq_cfg)) {
        std::cerr << "Error: failed to load DAQ config"
                  << (daq_config_file.empty() ? " (not found)" : ": " + daq_config_file)
                  << "\n";
        return 1;
    }
    std::cerr << "DAQ config: " << daq_config_file
              << " (adc_format=" << daq_cfg.adc_format << ")\n";

    EvChannel ch;
    ch.SetConfig(daq_cfg);
    if (ch.Open(path) != status::success) {
        std::cerr << "Failed to open: " << path << "\n";
        return 1;
    }

    int rc;
    if      (mode == "tree")       rc = doTree(ch, num);
    else if (mode == "tags")       rc = doTags(ch);
    else if (mode == "epics")      rc = doEpics(ch);
    else if (mode == "triggers")   rc = doTriggers(ch, verbose);
    else if (mode == "trig-debug") rc = doTrigDebug(ch, verbose);
    else if (mode == "event")      rc = doEvent(ch, num);
    else                           rc = doSummary(ch);

    ch.Close();
    return rc;
}
