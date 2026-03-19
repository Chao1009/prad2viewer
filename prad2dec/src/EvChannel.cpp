#include "EvChannel.h"
#include "Fadc250Decoder.h"
#include "evio.h"
#include <cstring>
#include <iostream>
#include <iomanip>

using namespace evc;

// --- evio C library status --------------------------------------------------
static inline status evio_status(int code)
{
    if (static_cast<unsigned>(code) == S_EVFILE_UNXPTDEOF) return status::incomplete;
    switch (code) {
    case S_SUCCESS:      return status::success;
    case EOF:            return status::eof;
    case S_EVFILE_TRUNC: return status::incomplete;
    default:             return status::failure;
    }
}

// --- open / close / read ----------------------------------------------------
EvChannel::EvChannel(size_t buflen) : fHandle(-1) { buffer.resize(buflen); }

status EvChannel::Open(const std::string &path)
{
    if (fHandle > 0) Close();
    char *cp = strdup(path.c_str()), *cm = strdup("r");
    int st = evOpen(cp, cm, &fHandle);
    free(cp); free(cm);
    return evio_status(st);
}

void EvChannel::Close() { evClose(fHandle); fHandle = -1; }
status EvChannel::Read() { return evio_status(evRead(fHandle, buffer.data(), buffer.size())); }

// === Scan ===================================================================
bool EvChannel::Scan()
{
    nodes.clear();
    BankHeader evh(&buffer[0]);
    if (evh.length + 1 > buffer.size()) return false;

    scanBank(0, 0, -1);

    // Classify event type using DaqConfig
    evtype = classify_event(evh.tag, config);

    // Determine number of events in this buffer.
    if (config.is_control(evh.tag)) {
        nevents = 0;
    } else if (config.is_physics(evh.tag)) {
        // CODA built-trigger (0xFF50-0xFF8F): num = event count in block.
        // Single-event mode (0xFE etc.): num = session ID, always 1 event.
        if (evh.tag >= 0xFF50 && evh.tag <= 0xFF8F)
            nevents = std::max<int>(evh.num, 1);
        else
            nevents = 1;
    } else {
        // EPICS, sync, and other non-physics events: single "event"
        nevents = (evtype == EventType::Epics || evtype == EventType::Sync) ? 1 : 0;
    }

    return true;
}

// --- scan a BANK (2-word header) --------------------------------------------
size_t EvChannel::scanBank(size_t off, int depth, int parent)
{
    BankHeader h(&buffer[off]);
    size_t total = h.length + 1;

    int idx = static_cast<int>(nodes.size());
    nodes.push_back({h.tag, h.type, h.num, depth, parent,
                     off + BankHeader::size(), h.data_words(), 0, 0});

    if (IsContainer(h.type)) {
        scanChildren(off + BankHeader::size(), h.data_words(), h.type, depth + 1, idx);
    } else if (h.type == DATA_COMPOSITE) {
        size_t doff = off + BankHeader::size();
        size_t dwords = h.data_words();
        size_t first_child = nodes.size();

        if (dwords >= 1) {
            size_t consumed = scanTagSegment(doff, depth + 1, idx);
            if (consumed < dwords)
                scanBank(doff + consumed, depth + 1, idx);
        }

        nodes[idx].child_first = first_child;
        nodes[idx].child_count = nodes.size() - first_child;
    }
    return total;
}

// --- scan a SEGMENT (1-word header) -----------------------------------------
size_t EvChannel::scanSegment(size_t off, int depth, int parent)
{
    SegmentHeader h(&buffer[off]);
    size_t total = 1 + h.length;

    int idx = static_cast<int>(nodes.size());
    nodes.push_back({h.tag, h.type, 0, depth, parent,
                     off + 1, h.length, 0, 0});

    if (IsContainer(h.type))
        scanChildren(off + 1, h.length, h.type, depth + 1, idx);
    return total;
}

// --- scan a TAGSEGMENT (1-word header) --------------------------------------
size_t EvChannel::scanTagSegment(size_t off, int depth, int parent)
{
    TagSegmentHeader h(&buffer[off]);
    size_t total = 1 + h.length;

    int idx = static_cast<int>(nodes.size());
    nodes.push_back({h.tag, h.type, 0, depth, parent,
                     off + 1, h.length, 0, 0});

    if (IsContainer(h.type))
        scanChildren(off + 1, h.length, h.type, depth + 1, idx);
    return total;
}

// --- scan children of a container -------------------------------------------
void EvChannel::scanChildren(size_t off, size_t nwords, uint32_t ptype, int depth, int pidx)
{
    size_t first_child = nodes.size();
    size_t count = 0, pos = 0;

    while (pos < nwords) {
        size_t consumed = 0;
        switch (ptype) {
        case DATA_BANK: case DATA_BANK2:
            consumed = scanBank(off + pos, depth, pidx); break;
        case DATA_SEGMENT: case DATA_SEGMENT2:
            consumed = scanSegment(off + pos, depth, pidx); break;
        case DATA_TAGSEGMENT:
            consumed = scanTagSegment(off + pos, depth, pidx); break;
        default: return;
        }
        if (consumed == 0) break;
        pos += consumed;
        ++count;
    }

    nodes[pidx].child_first = first_child;
    nodes[pidx].child_count = count;
}

// === accessors ==============================================================

std::vector<const EvNode*> EvChannel::FindByTag(uint32_t tag) const
{
    std::vector<const EvNode*> result;
    for (auto &n : nodes)
        if (n.tag == tag) result.push_back(&n);
    return result;
}

const EvNode *EvChannel::FindFirstByTag(uint32_t tag) const
{
    for (auto &n : nodes)
        if (n.tag == tag) return &n;
    return nullptr;
}

const uint8_t *EvChannel::GetCompositePayload(const EvNode &n, size_t &nbytes) const
{
    nbytes = 0;
    if (n.type != DATA_COMPOSITE || n.child_count < 2) return nullptr;
    auto &inner = nodes[n.child_first + 1];
    nbytes = inner.data_words * sizeof(uint32_t);
    return reinterpret_cast<const uint8_t*>(&buffer[inner.data_begin]);
}

// === CODA built trigger bank decoding (spec pages 21, 26, 31) ===============
//
// Physics Event structure:
//   Top-level bank (tag 0xFFXX, type 0x10, num = M events)
//     Built Trigger Bank (tag 0xFF2X, type 0x20 = bank of segments, num = N ROCs)
//       Segment 1 (tag=EB_id, type=0xa ULONG64): first_event_number(64), [avg_ts(64)]xM, [run_info(64)]
//       Segment 2 (tag=EB_id, type=0x05 USHORT16): event_type[0..M-1]
//       ROC segment (tag=roc_id, type=0x01 UINT32): [ts_low, ts_high, misc] per event
//     Data Bank 1 (ROC 1 data)
//     ...
//
// Trigger bank tag encodes content (page 26):
//   bit 0: has timestamps, bit 1: has run#/type, bit 2: no run-specific data

bool EvChannel::decodeTriggerBank(int event_idx, fdec::EventInfo &info) const
{
    // find the built trigger bank at depth 1 (direct child of top-level event)
    const EvNode *tb = nullptr;
    for (auto &n : nodes) {
        if (n.depth == 1 && DaqConfig::is_built_trigger_bank(n.tag)) {
            tb = &n;
            break;
        }
    }
    // fallback: try raw trigger bank (0xFF1X) from ROC raw data
    if (!tb) {
        for (auto &n : nodes) {
            if (n.depth == 1 && DaqConfig::is_raw_trigger_bank(n.tag)) {
                tb = &n;
                break;
            }
        }
    }
    if (!tb || tb->child_count < 1) return false;

    uint32_t tb_tag = tb->tag;
    bool is_built = DaqConfig::is_built_trigger_bank(tb_tag);
    bool has_timestamps = DaqConfig::trigger_bank_has_timestamps(tb_tag);
    bool has_run_info   = DaqConfig::trigger_bank_has_run_info(tb_tag);

    if (is_built) {
        // --- Built trigger bank (0xFF2X) per page 31 ---
        if (tb->child_count < 2) return false;

        // Segment 1: common data (type 0xa = ULONG64)
        // Contains: first_event_number(64-bit), [avg_timestamps(64-bit) x M], [run#_type(64-bit)]
        auto &seg1 = nodes[tb->child_first];
        if (seg1.type == DATA_ULONG64 && seg1.data_words >= 2) {
            const uint32_t *d = GetData(seg1);

            // 64-bit first event number (low word first)
            uint64_t first_event = (uint64_t)d[0] | ((uint64_t)d[1] << 32);
            info.event_number = static_cast<int32_t>(first_event + event_idx);

            size_t pos = 2; // past event number (2 words for 64-bit)

            if (has_timestamps) {
                // M average timestamps, each 64-bit (2 words)
                size_t ts_off = pos + 2 * event_idx;
                if (ts_off + 1 < seg1.data_words) {
                    uint64_t ts_low  = d[ts_off];
                    uint64_t ts_high = d[ts_off + 1];
                    info.timestamp = ts_low | (ts_high << 32);
                }
                pos += 2 * nevents;
            }

            if (has_run_info && pos + 1 < seg1.data_words) {
                // 64-bit value: run# in high 32, run type in low 32
                info.run_number = d[pos + 1];  // high word = run number
                // d[pos] = run type (low word), not stored currently
            }
        }

        // Segment 2: event types (type 0x05 = USHORT16)
        // Packed as unsigned shorts: 2 per 32-bit word, low 16 first
        auto &seg2 = nodes[tb->child_first + 1];
        if (seg2.type == DATA_USHORT16 && seg2.data_words > 0) {
            const uint32_t *d = GetData(seg2);
            int word_idx = event_idx / 2;
            int shift    = (event_idx & 1) * 16;
            if (static_cast<size_t>(word_idx) < seg2.data_words) {
                info.trigger_bits = static_cast<uint8_t>((d[word_idx] >> shift) & 0xFF);
            }
        }

        // ROC segments (type 0x01 = UINT32): per-ROC timestamps
        // These start at child_first + 2 (after the 2 common segments)
        // Each ROC segment has: [ts_low, ts_high, misc] per event
        // We can use these for more precise per-ROC timestamps if avg is missing
        if (!has_timestamps && tb->child_count > 2) {
            auto &roc_seg = nodes[tb->child_first + 2]; // first ROC
            if (roc_seg.type == DATA_UINT32) {
                const uint32_t *d = GetData(roc_seg);
                // words per event: at least 2 (ts_low, ts_high), possibly more (misc)
                size_t words_per_evt = (nevents > 0) ? roc_seg.data_words / nevents : 0;
                if (words_per_evt >= 2) {
                    size_t off = words_per_evt * event_idx;
                    if (off + 1 < roc_seg.data_words) {
                        uint64_t ts_low  = d[off];
                        uint64_t ts_high = d[off + 1] & 0xFFFF; // upper 16 bits of 48-bit ts
                        info.timestamp = ts_low | (ts_high << 32);
                    }
                }
            }
        }

    } else {
        // --- Raw trigger bank (0xFF1X) per page 21 ---
        // Children are segments, one per event:
        //   segment tag = event ID / trigger type
        //   segment data: event_number, [ts_low, ts_high], [misc]
        if (event_idx >= static_cast<int>(tb->child_count)) return false;

        auto &seg = nodes[tb->child_first + event_idx];
        info.trigger_bits = static_cast<uint8_t>(seg.tag); // trigger type in segment tag

        const uint32_t *d = GetData(seg);
        size_t nw = seg.data_words;
        if (nw >= 1)
            info.event_number = static_cast<int32_t>(d[0]);
        if (has_timestamps && nw >= 3) {
            uint64_t ts_low  = d[1];
            uint64_t ts_high = d[2] & 0xFFFF; // bits 47-32
            info.timestamp = ts_low | (ts_high << 32);
        }
    }

    return true;
}

// === JLab TI bank fallback (0xE10A) =========================================
// Used when no CODA trigger bank is present (single-event mode, legacy data)

bool EvChannel::decodeTI(fdec::EventInfo &info) const
{
    // --- TI data bank (0xE10A): trigger type, trigger number, timestamp -----
    auto *ti = FindFirstByTag(config.ti_bank_tag);
    if (ti) {
        const uint32_t *d = GetData(*ti);
        size_t nw = ti->data_words;

        // raw trigger bits
        if (config.ti_trigger_type_word >= 0 &&
            static_cast<size_t>(config.ti_trigger_type_word) < nw)
        {
            info.trigger_bits = (d[config.ti_trigger_type_word] >> config.ti_trigger_type_shift)
                                & config.ti_trigger_type_mask;
        }

        if (config.ti_trigger_word >= 0 &&
            static_cast<size_t>(config.ti_trigger_word) < nw)
            info.trigger_number = static_cast<int32_t>(d[config.ti_trigger_word]);

        int lo = config.ti_time_low_word;
        int hi = config.ti_time_high_word;
        if (lo >= 0 && hi >= 0 &&
            static_cast<size_t>(lo) < nw &&
            static_cast<size_t>(hi) < nw)
        {
            uint64_t time_low  = d[lo];
            uint64_t time_high = (d[hi] & config.ti_time_high_mask);
            if (config.ti_time_high_shift > 0)
                time_high >>= config.ti_time_high_shift;
            info.timestamp = (time_high << 32) | time_low;
        }
    }

    // --- run info bank (0xE10F, in TI master crate) -------------------------
    if (auto *ri = FindFirstByTag(config.run_info_tag)) {
        const uint32_t *d = GetData(*ri);
        size_t nw = ri->data_words;

        if (config.ri_run_number_word >= 0 &&
            static_cast<size_t>(config.ri_run_number_word) < nw)
            info.run_number = d[config.ri_run_number_word];

        if (config.ri_unix_time_word >= 0 &&
            static_cast<size_t>(config.ri_unix_time_word) < nw)
            info.unix_time = d[config.ri_unix_time_word];
    }

    return ti != nullptr;
}

// === DecodeEvent ============================================================

bool EvChannel::DecodeEvent(int i, fdec::EventData &evt) const
{
    evt.clear();
    if (i < 0 || i >= nevents) return false;

    BankHeader evh(&buffer[0]);

    // fill event info
    evt.info.event_tag = evh.tag;
    evt.info.type      = static_cast<uint8_t>(evtype);

    // Try CODA trigger bank first (0xFF2X / 0xFF1X), then fall back to JLab TI
    if (!decodeTriggerBank(i, evt.info))
        decodeTI(evt.info);

    // find all FADC composite banks matching configured tag
    int roc_idx = 0;
    for (size_t ni = 0; ni < nodes.size() && roc_idx < fdec::MAX_ROCS; ++ni) {
        auto &n = nodes[ni];
        if (n.tag != config.fadc_composite_tag || n.type != DATA_COMPOSITE)
            continue;

        size_t nbytes;
        auto *payload = GetCompositePayload(n, nbytes);
        if (!payload) continue;

        // find parent ROC tag
        uint32_t roc_tag = (n.parent >= 0) ? nodes[n.parent].tag : 0;

        fdec::RocData &roc = evt.rocs[roc_idx];
        // evt.clear() already zeroed everything — just set active fields
        roc.present = true;
        roc.tag = roc_tag;

        fdec::Fadc250Decoder::DecodeRoc(payload, nbytes, roc);

        evt.roc_index[roc_idx] = roc_idx;
        roc_idx++;
    }
    evt.nrocs = roc_idx;
    return roc_idx > 0;
}

// === EPICS text extraction ==================================================

std::string EvChannel::ExtractEpicsText() const
{
    // look for the EPICS bank by configured tag
    auto epics_nodes = FindByTag(config.epics_bank_tag);

    // fallback: if no bank with epics_bank_tag, try string-type banks
    // at depth 1 (direct children of the event)
    if (epics_nodes.empty()) {
        for (auto &n : nodes) {
            if (n.depth == 1 &&
                (n.type == DATA_CHARSTAR8 || n.type == DATA_CHAR8) &&
                n.data_words > 0)
            {
                epics_nodes.push_back(&n);
            }
        }
    }

    if (epics_nodes.empty()) return {};

    // extract text from the first matching node
    const EvNode &n = *epics_nodes[0];
    const char *raw = reinterpret_cast<const char*>(&buffer[n.data_begin]);
    size_t max_len = n.data_words * sizeof(uint32_t);

    // find actual string length (may be null-padded)
    size_t len = 0;
    while (len < max_len && raw[len] != '\0') ++len;

    return std::string(raw, len);
}

// === PrintTree ==============================================================
void EvChannel::PrintTree(std::ostream &os) const
{
    for (auto &n : nodes) {
        for (int i = 0; i < n.depth; ++i) os << "  ";

        os << std::setw(6) << std::left << TypeName(n.type) << std::right
           << " tag=0x" << std::hex << n.tag << std::dec << "(" << n.tag << ")"
           << " type=0x" << std::hex << n.type << std::dec
           << " num=" << n.num
           << " data=" << n.data_words << "w";

        if (n.child_count > 0)
            os << " children=" << n.child_count;

        if (n.child_count == 0 && n.data_words > 0 && !IsContainer(n.type) && n.type != DATA_COMPOSITE) {
            os << " |";
            size_t nshow = std::min<size_t>(n.data_words, 4);
            for (size_t i = 0; i < nshow; ++i)
                os << " " << std::hex << std::setw(8) << std::setfill('0')
                   << buffer[n.data_begin + i] << std::setfill(' ') << std::dec;
            if (n.data_words > nshow) os << " ...";
        }

        if ((n.type == DATA_CHARSTAR8 || n.type == DATA_CHAR8) && n.data_words > 0) {
            const char *s = reinterpret_cast<const char*>(&buffer[n.data_begin]);
            size_t maxlen = n.data_words * 4;
            os << " \"";
            for (size_t i = 0; i < maxlen && s[i]; ++i) {
                if (s[i] >= 32 && s[i] < 127) os << s[i];
                else os << '.';
            }
            os << "\"";
        }
        os << "\n";
    }
}
