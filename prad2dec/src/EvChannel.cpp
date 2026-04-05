#include "EvChannel.h"
#include "Fadc250Decoder.h"
#include "Fadc250RawDecoder.h"
#include "Adc1881mDecoder.h"
#include "SspDecoder.h"
#include "evio.h"
#include <cstring>
#include <iostream>
#include <iomanip>
#include <set>

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

// =============================================================================
// Depth-2 data bank decoders
// =============================================================================

// --- 0xC000: CODA trigger bank [event_number, event_tag, reserved] ----------
void EvChannel::decodeTriggerInfo(const EvNode &node, fdec::EventInfo &info) const
{
    const uint32_t *d = GetData(node);
    size_t nw = node.data_words;
    if (config.trig_event_number_word >= 0 &&
        static_cast<size_t>(config.trig_event_number_word) < nw)
        info.event_number = static_cast<int32_t>(d[config.trig_event_number_word]);
}

// --- 0xE10A: TI hardware data (after rol2 block-header stripping) -----------
//
// Slave (4 words, nwords=3):
//   d[0]: event_header  — event_type(8) | 0x01(8) | nwords(16)
//   d[1]: event_number  — 32-bit
//   d[2]: timestamp_low — 32-bit
//   d[3]: ts_high[15:0] | evnum_high[19:16]
//
// Master (7 words, nwords=6):
//   d[0]-d[3]: same as slave
//   d[4]: trigger type byte (often zero)
//   d[5]: 32-bit FP trigger inputs (if tiSetFPInputReadout enabled)
//   d[6]: additional TI flags
//
void EvChannel::decodeTIBank(const EvNode &node, fdec::EventInfo &info,
                             bool is_master) const
{
    const uint32_t *d = GetData(node);
    size_t nw = node.data_words;
    if (nw == 0) return;

    // d[1] event/trigger number
    if (config.ti_trigger_word >= 0 &&
        static_cast<size_t>(config.ti_trigger_word) < nw)
        info.trigger_number = static_cast<int32_t>(d[config.ti_trigger_word]);

    // d[2], d[3]: 48-bit timestamp
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

    // d[5] (TI master only): 32-bit FP trigger input snapshot.
    // Multiple bits can fire simultaneously — tells you which detector
    // signals were active at trigger time. Independent of trigger_type.
    // See database/trigger_bit.json "trigger_bits" section.
    if (is_master && config.ti_trigger_type_word >= 0 &&
        static_cast<size_t>(config.ti_trigger_type_word) < nw)
    {
        info.trigger_bits = (d[config.ti_trigger_type_word]
                             >> config.ti_trigger_type_shift)
                            & config.ti_trigger_type_mask;
    }
}

// --- 0xE10F: Run info bank [hdr, run#, evt_count, unix_time, ...] -----------
void EvChannel::decodeRunInfo(const EvNode &node, fdec::EventInfo &info) const
{
    const uint32_t *d = GetData(node);
    size_t nw = node.data_words;
    if (config.ri_run_number_word >= 0 &&
        static_cast<size_t>(config.ri_run_number_word) < nw)
        info.run_number = d[config.ri_run_number_word];
    if (config.ri_unix_time_word >= 0 &&
        static_cast<size_t>(config.ri_unix_time_word) < nw)
        info.unix_time = d[config.ri_unix_time_word];
}

// =============================================================================
// DecodeEvent — hierarchical dispatch by depth and bank tag
//
// CODA2 single-event structure (see docs/rols/banktags.md):
//
//   [Physics Event: tag = 0x80 + TI_event_type] (depth 0)
//    +-- [0xC000 UINT32]       trigger bank (event#, event_tag)
//    +-- [0x0027 BANK]         TI master crate
//    |    +-- [0xE10A UINT32]  TI data (7 words: trigger#, timestamp, FP bits)
//    |    +-- [0xE10C UINT32]  SSP trigger data
//    |    +-- [0xE10F UINT32]  run info (run#, unix_time)
//    |    +-- [0xE10E STRING]  DAQ config string
//    +-- [0x0080 BANK]         HyCal FADC crate (even tags)
//    |    +-- [0xE10A UINT32]  TI data (4 words)
//    |    +-- [0xE101 COMPOSITE] FADC250 waveforms (physics triggers only)
//    +-- [0x0081 BANK]         TI slave (odd tags, TI only)
//    |    +-- [0xE10A UINT32]  TI data
//    +-- [0x0031 BANK]         GEM VTP/MPD crate (when present)
//    |    +-- [0xE10A UINT32]  TI data
//    |    +-- [0x0DEA UINT32]  MPD strip data (SSP bitfield format)
//    +-- ...
// =============================================================================

// --- Known bank tags from JLab ROLs (rol1.c, rol2.c, vtp1mpd.c) ------------
// Used to provide informative warnings for banks we recognize but don't decode.
struct KnownBankTag {
    uint32_t    tag;
    const char *hardware;
    bool        has_decoder;   // true if prad2dec implements a decoder
};

static const KnownBankTag known_bank_tags[] = {
    // Tags with decoders
    { 0xE10A, "TI/TS",                    true  },
    { 0xE101, "FADC250 composite",         true  },
    { 0xE109, "FADC250 raw",               true  },
    { 0xE120, "ADC1881M (Fastbus)",        true  },
    { 0xE10C, "SSP",                       true  },
    { 0x0DEA, "VTP/MPD (GEM)",             true  },
    { 0xE10F, "Run info",                  true  },
    { 0xE10E, "DAQ config string",         true  },
    { 0xC000, "CODA trigger bank",         true  },
    // Tags from ROLs without decoders
    { 0xE10B, "V1190/V1290 TDC",           false },
    { 0xE141, "FAV3 (FADC v3)",            false },
    { 0xE104, "VSCM",                      false },
    { 0xE105, "DCRB/DC/Vetroc",            false },
    { 0xE115, "DSC2 scaler",               false },
    { 0xE112, "HEAD bank",                 false },
    { 0xE123, "SSP-RICH",                  false },
    { 0xE125, "Per-slot data",             false },
    { 0xE131, "VFTDC",                     false },
    { 0xE133, "Helicity Decoder",          false },
    { 0xE140, "Special (pid=0)",           false },
};
static const int N_KNOWN_TAGS = sizeof(known_bank_tags) / sizeof(known_bank_tags[0]);

static const KnownBankTag *lookupKnownTag(uint32_t tag)
{
    for (int i = 0; i < N_KNOWN_TAGS; ++i)
        if (known_bank_tags[i].tag == tag) return &known_bank_tags[i];
    return nullptr;
}

bool EvChannel::DecodeEvent(int i, fdec::EventData &evt,
                            ssp::SspEventData *ssp_evt) const
{
    evt.clear();
    if (ssp_evt) ssp_evt->clear();
    if (i < 0 || i >= nevents) return false;
    if (nodes.empty()) return false;

    BankHeader evh(&buffer[0]);
    evt.info.event_tag = evh.tag;
    evt.info.type = static_cast<uint8_t>(evtype);

    // Main trigger type: derived from event tag (see docs/rols/banktags.md).
    // event_tag = physics_base + trigger_type (physics_base from daq_config.json).
    // trigger_type identifies WHICH trigger caused this event (single per event).
    // Maps to trigger name via database/trigger_bits.json "trigger_type" section.
    if (config.is_physics(evh.tag) && evh.tag >= config.physics_base)
        evt.info.trigger_type = static_cast<uint8_t>(evh.tag - config.physics_base);

    // --- extract event info: scan all nodes for trigger/TI banks ------------
    // Works for both built (3-level), flat (2-level), and mixed structures.

    // 0xC000: trigger bank → event number
    if (auto *tb = FindFirstByTag(config.trigger_bank_tag))
        decodeTriggerInfo(*tb, evt.info);

    // Scan all TI banks (0xE10A) for event info.
    // The first TI bank found provides trigger_number and timestamp.
    // The longest TI bank (TI master, 7 words) provides FP trigger_bits from d[5].
    // Run info (0xE10F) is extracted from the TI master crate (0x0027).
    bool have_info = false;
    for (auto &n : nodes) {
        if (n.tag != config.ti_bank_tag || n.type != DATA_UINT32) continue;

        // first TI bank: extract trigger_number + timestamp
        if (!have_info) {
            decodeTIBank(n, evt.info, false);
            have_info = true;
        }

        // any TI bank with enough words: extract FP trigger bits from d[5]
        if (evt.info.trigger_bits == 0 && config.ti_trigger_type_word >= 0 &&
            static_cast<size_t>(config.ti_trigger_type_word) < n.data_words)
        {
            const uint32_t *d = GetData(n);
            evt.info.trigger_bits = (d[config.ti_trigger_type_word]
                                     >> config.ti_trigger_type_shift)
                                    & config.ti_trigger_type_mask;
        }
    }

    // TI master crate (0x27): extract run info
    for (auto &n : nodes) {
        if (n.depth == 1 && n.tag == config.ti_master_tag) {
            for (size_t ci = 0; ci < n.child_count; ++ci) {
                auto &child = nodes[n.child_first + ci];
                if (child.tag == config.run_info_tag && child.type == DATA_UINT32)
                    decodeRunInfo(child, evt.info);
            }
            break;
        }
    }

    // --- decode detector data: flat scan all nodes by tag -------------------
    // Works for built (3-level), flat (2-level), and mixed event structures.
    int roc_idx = 0;
    bool ssp_decoded = false;
    static std::set<uint64_t> warned_tags;

    for (size_t ni = 0; ni < nodes.size(); ++ni) {
        auto &n = nodes[ni];

        // === phase 1: skip nodes that are not dispatchable data banks ===
        if (IsContainer(n.type))    continue;  // ROC wrappers, event bank
        if (n.data_words == 0)      continue;  // empty banks
        if (n.parent >= 0 &&                   // composite internals (0x000D, 0x0000)
            nodes[n.parent].type == DATA_COMPOSITE) continue;
        if (n.tag == config.ti_bank_tag)        continue;  // handled above
        if (n.tag == config.trigger_bank_tag)   continue;  // handled above
        if (n.tag == config.run_info_tag)       continue;  // handled above
        if (n.tag == config.daq_config_tag)     continue;  // config string, not data
        if (n.type == DATA_CHARSTAR8 || n.type == DATA_CHAR8) continue;

        // parent ROC tag (for crate_id lookup; 0 if flat/top-level)
        uint32_t roc_tag = (n.parent >= 0) ? nodes[n.parent].tag : 0;

        // === phase 2: dispatch by tag to the appropriate decoder ===

        // FADC250 composite waveforms (0xE101)
        if (n.tag == config.fadc_composite_tag && n.type == DATA_COMPOSITE
            && roc_idx < fdec::MAX_ROCS)
        {
            size_t nbytes;
            auto *payload = GetCompositePayload(n, nbytes);
            if (!payload) continue;
            fdec::RocData &rd = evt.rocs[roc_idx];
            rd.present = true;
            rd.tag = roc_tag;
            fdec::Fadc250Decoder::DecodeRoc(payload, nbytes, rd);
            evt.roc_index[roc_idx] = roc_idx;
            roc_idx++;
            continue;
        }

        // FADC250 raw hardware format (0xE109, fallback)
        if (n.tag == config.fadc_raw_tag && n.type == DATA_UINT32
            && roc_idx < fdec::MAX_ROCS)
        {
            fdec::RocData &rd = evt.rocs[roc_idx];
            rd.present = true;
            rd.tag = roc_tag;
            fdec::Fadc250RawDecoder::DecodeRoc(GetData(n), n.data_words, rd);
            evt.roc_index[roc_idx] = roc_idx;
            roc_idx++;
            continue;
        }

        // SSP/MPD data — GEM (0xE10C, 0x0DEA)
        if (ssp_evt && config.is_ssp_bank(n.tag) && n.type == DATA_UINT32)
        {
            int crate_id = -1;
            for (auto &re : config.roc_tags)
                if (re.tag == roc_tag) { crate_id = re.crate; break; }
            int napvs = ssp::SspDecoder::DecodeRoc(GetData(n), n.data_words,
                                                    crate_id, *ssp_evt);
            if (napvs > 0) ssp_decoded = true;  // only count if actual APV data found
            continue;
        }

        // ADC1881M raw data — PRad legacy
        if (config.adc_format == "adc1881m"
            && n.tag == config.adc1881m_bank_tag
            && roc_idx < fdec::MAX_ROCS)
        {
            int crate_id = -1;
            for (auto &re : config.roc_tags)
                if (re.tag == roc_tag) { crate_id = re.crate; break; }
            fdec::RocData &rd = evt.rocs[roc_idx];
            rd.present = true;
            rd.tag = roc_tag;
            fdec::Adc1881mDecoder::DecodeRoc(GetData(n), n.data_words, rd);

            if (!config.pedestals.empty() && crate_id >= 0) {
                for (int s = 0; s < fdec::MAX_SLOTS; ++s) {
                    auto &slot = rd.slots[s];
                    if (!slot.present) continue;
                    for (int c = 0; c < fdec::MAX_CHANNELS; ++c) {
                        if (!(slot.channel_mask & (1ull << c))) continue;
                        auto &cd = slot.channels[c];
                        if (cd.nsamples != 1) continue;
                        auto *ped = config.get_pedestal(crate_id, s, c);
                        if (!ped) continue;
                        float raw = static_cast<float>(cd.samples[0]);
                        float threshold = ped->mean + config.sparsify_sigma * ped->rms;
                        if (config.sparsify_sigma > 0.f && raw < threshold) {
                            cd.nsamples = 0;
                            slot.channel_mask &= ~(1ull << c);
                            slot.nchannels--;
                        } else {
                            int sub = static_cast<int>(raw) - static_cast<int>(ped->mean + 0.5f);
                            cd.samples[0] = (sub > 0) ? static_cast<uint16_t>(sub) : 0;
                        }
                    }
                }
            }
            evt.roc_index[roc_idx] = roc_idx;
            roc_idx++;
            continue;
        }

        // === phase 3: warn about unhandled data banks (once per tag) ===
        uint64_t key = (uint64_t(roc_tag) << 32) | n.tag;
        if (warned_tags.insert(key).second) {
            auto *known = lookupKnownTag(n.tag);
            if (known && !known->has_decoder)
                std::cerr << "DecodeEvent: skipping " << known->hardware
                          << " bank (0x" << std::hex << n.tag << std::dec
                          << ", " << n.data_words << "w)"
                          << " in ROC 0x" << std::hex << roc_tag << std::dec
                          << " — no decoder implemented\n";
            else if (!known)
                std::cerr << "DecodeEvent: unknown bank"
                          << " tag=0x" << std::hex << n.tag
                          << " type=" << TypeName(n.type)
                          << std::dec << " (" << n.data_words << "w)"
                          << " in ROC 0x" << std::hex << roc_tag << std::dec
                          << " — not in known tags, check ROL source\n";
        }
    }

    evt.nrocs = roc_idx;
    // Return true only when actual detector data was decoded (FADC waveforms
    // or GEM strips).  Event info (event#, timestamp, trigger_type/bits) is
    // always populated in evt.info regardless — callers can use it even when
    // this returns false (e.g. monitoring events with TI data but no waveforms).
    return roc_idx > 0 || ssp_decoded;
}

// === Control event time extraction ==========================================

uint32_t EvChannel::GetControlTime() const
{
    // All CODA control events (Prestart, Go, Sync, End) share the same layout
    if (evtype != EventType::Sync && evtype != EventType::Prestart &&
        evtype != EventType::Go && evtype != EventType::End)
        return 0;

    // Control event layout (after 2-word bank header):
    //   word[0]: [Event Type | 0x01 | 0]  (data word header)
    //   word[1]: unix timestamp
    //   word[2]: A (run number / event counts)
    //   word[3]: B (run type / event counts)
    BankHeader evh(&buffer[0]);
    size_t data_off = BankHeader::size();
    size_t data_words = evh.data_words();
    if (data_words >= 2)
        return buffer[data_off + 1];  // second data word is time
    return 0;
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
