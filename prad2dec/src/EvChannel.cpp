#include "EvChannel.h"
#include "Fadc250Decoder.h"
#include "Fadc250RawDecoder.h"
#include "Adc1881mDecoder.h"
#include "SspDecoder.h"
#include "VtpDecoder.h"
#include "TdcDecoder.h"
#include "Dsc2Decoder.h"
#include "evio.h"
#include <cstring>
#include <iostream>
#include <iomanip>
#include <mutex>
#include <set>

using namespace evc;

// The EVIO C library maintains a global handle table (handleList[]) that is
// NOT thread-safe.  evOpen and evClose both write to this table, so concurrent
// calls from multiple threads corrupt it, causing crashes in evReadRandom,
// SspDecoder, and GemSystem on all threads.  Serialize all evOpen/evClose
// calls with a process-wide mutex.  evRead/evReadRandom on distinct handles
// are safe to call concurrently once the handles are established.
static std::mutex g_evio_open_mutex;

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

status EvChannel::OpenSequential(const std::string &path)
{
    if (fHandle > 0) Close();
    char *cp = strdup(path.c_str()), *cm = strdup("r");
    int st;
    { std::lock_guard<std::mutex> lk(g_evio_open_mutex); st = evOpen(cp, cm, &fHandle); }
    free(cp); free(cm);
    return evio_status(st);
}

void EvChannel::Close()
{
    std::lock_guard<std::mutex> lk(g_evio_open_mutex);
    evClose(fHandle);
    fHandle = -1;
    ra_count = 0;
    ra_pos   = 0;
}

// Unified sequential-style read: in RA mode advances an internal cursor
// via ReadEventByIndex, in sequential mode falls through to evRead.  This
// lets legacy callers (e.g. iterateAll, analysis tools) that follow an
// Open() with a `while (Read() == success)` loop keep working regardless
// of which mode the handle was opened in.
status EvChannel::Read()
{
    if (ra_count > 0) {
        if (ra_pos >= ra_count) return status::eof;
        return ReadEventByIndex(ra_pos);   // bumps ra_pos on success
    }
    return evio_status(evRead(fHandle, buffer.data(), buffer.size()));
}

// --- random-access open (evio "ra" mode) ------------------------------------
status EvChannel::OpenRandomAccess(const std::string &path)
{
    if (fHandle > 0) Close();
    ra_count = 0;
    ra_pos   = 0;
    char *cp = strdup(path.c_str()), *cm = strdup("ra");
    int st;
    { std::lock_guard<std::mutex> lk(g_evio_open_mutex); st = evOpen(cp, cm, &fHandle); }
    free(cp); free(cm);
    if (st != S_SUCCESS) return evio_status(st);

    // Retrieve the event count from the random-access table.  The table
    // itself is owned by evio — we only want the size.  The pointer type
    // differs between evio versions:
    //   evio-4: int evGetRandomAccessTable(int, const uint32_t ***, uint32_t *);
    //   evio-6: int evGetRandomAccessTable(int, uint32_t *** const,  uint32_t *);
    // (the top-level `const` on the 6.0 parameter is cosmetic; the pointee
    // type is what breaks compatibility). EV_VERSION is declared by both
    // headers so we can pick the matching local type.
#if defined(EV_VERSION) && EV_VERSION >= 6
    uint32_t **table = nullptr;
#else
    const uint32_t **table = nullptr;
#endif
    uint32_t n = 0;
    int r = evGetRandomAccessTable(fHandle, &table, &n);
    if (r != S_SUCCESS) {
        Close();
        return evio_status(r);
    }
    ra_count = static_cast<int>(n);
    return status::success;
}

status EvChannel::ReadEventByIndex(int evio_event_index)
{
    if (fHandle <= 0 || ra_count == 0) return status::failure;
    if (evio_event_index < 0 || evio_event_index >= ra_count) return status::failure;

    const uint32_t *pEvent = nullptr;
    uint32_t buflen = 0;
    // evReadRandom is 1-indexed per the evio API; callers use 0-based.
    int r = evReadRandom(fHandle, &pEvent, &buflen,
                          static_cast<uint32_t>(evio_event_index + 1));
    if (r != S_SUCCESS || pEvent == nullptr) return evio_status(r);

    if (buflen > buffer.size()) buffer.resize(buflen);
    std::memcpy(buffer.data(), pEvent, buflen * sizeof(uint32_t));
    ra_pos = evio_event_index + 1;   // so Read() picks up from the next one
    return status::success;
}

// --- auto-detect open: RA first, sequential fallback -----------------------
status EvChannel::OpenAuto(const std::string &path)
{
    status s = OpenRandomAccess(path);
    if (s == status::success) return s;
    // RA failed (e.g. file missing the optional block-index arrays).  Clean
    // up any half-open handle before trying the sequential mode.
    if (fHandle > 0) Close();
    return OpenSequential(path);
}

// === SetConfig ==============================================================
// Back-fill data_banks from legacy tag fields for anything the JSON didn't
// declare, then precompute per-product tag lists consulted by the lazy
// accessors.  Keeps existing configs (no bank_structure section) working.
void EvChannel::SetConfig(const DaqConfig &cfg_in)
{
    config = cfg_in;

    auto ensure = [&](uint32_t tag, const char *mod, const char *prod) {
        if (tag != 0 && config.data_banks.find(tag) == config.data_banks.end())
            config.data_banks[tag] = DaqConfig::DataBankInfo{mod, prod, ""};
    };
    ensure(config.fadc_composite_tag, "fadc250_composite", DaqConfig::product_fadc);
    ensure(config.fadc_raw_tag,       "fadc250_raw",       DaqConfig::product_fadc);
    ensure(config.adc1881m_bank_tag,  "adc1881m",          DaqConfig::product_fadc);
    ensure(config.ti_bank_tag,        "ti",                DaqConfig::product_event_info);
    ensure(config.trigger_bank_tag,   "trigger_bank",      DaqConfig::product_event_info);
    ensure(config.run_info_tag,       "run_info",          DaqConfig::product_event_info);
    ensure(config.daq_config_tag,     "daq_config_string", DaqConfig::product_daq_config);
    ensure(config.tdc_bank_tag,       "v1190_tdc",         DaqConfig::product_tdc);
    ensure(config.epics_bank_tag,     "epics_data",        DaqConfig::product_epics);
    for (auto t : config.ssp_bank_tags)
        ensure(t, "ssp_mpd", DaqConfig::product_gem);
    // VTP ECAL peaks/clusters — no DaqConfig field yet, tag hardcoded.
    ensure(0xE122, "vtp_ecal", DaqConfig::product_vtp);

    fadc_tags = config.banks_for_product(DaqConfig::product_fadc);
    gem_tags  = config.banks_for_product(DaqConfig::product_gem);
    tdc_tags  = config.banks_for_product(DaqConfig::product_tdc);
    vtp_tags  = config.banks_for_product(DaqConfig::product_vtp);
}

// === Scan ===================================================================
bool EvChannel::Scan()
{
    nodes.clear();
    tag_index.clear();
    cached_event_idx = -1;
    clearCache();

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

    // Build tag index consulted by the lazy accessors for O(1) dispatch.
    for (size_t i = 0; i < nodes.size(); ++i)
        tag_index[nodes[i].tag].push_back(static_cast<int>(i));

    // Default-select sub-event 0 so Info()/Fadc()/etc. work without an
    // explicit SelectEvent() call in single-event mode (the common case).
    cached_event_idx = 0;

    // Fresh event → allow Sync() to attempt one decode if this event carries
    // SYNC / control data; the snapshot itself is only updated on a hit, so
    // physics events keep the prior absolute-time reference visible.
    sync_decoded_this_event_ = false;

    return true;
}

void EvChannel::clearCache() const
{
    info_ready = fadc_ready = gem_ready = tdc_ready = vtp_ready = false;
    dsc_ready = epics_ready = false;
}

void EvChannel::SelectEvent(int i) const
{
    if (cached_event_idx != i) {
        cached_event_idx = i;
        clearCache();
    }
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

const std::vector<int> &EvChannel::NodesForTag(uint32_t tag) const
{
    static const std::vector<int> empty;
    auto it = tag_index.find(tag);
    return it != tag_index.end() ? it->second : empty;
}

// === Lazy data-product accessors ============================================
// Each one decodes on first call (for the currently-selected sub-event) and
// returns a cached reference thereafter.  clearCache() invalidates the flags
// on Read()/Scan()/SelectEvent() transitions.

const fdec::EventInfo &EvChannel::Info() const
{
    if (!cache_fadc) cache_fadc = std::make_unique<fdec::EventData>();
    if (!info_ready) {
        decodeInfoInto(cache_fadc->info);
        info_ready = true;
    }
    return cache_fadc->info;
}

const fdec::EventData &EvChannel::Fadc() const
{
    if (!cache_fadc) cache_fadc = std::make_unique<fdec::EventData>();
    if (!fadc_ready) {
        decodeFadcInto(*cache_fadc);   // also refills cache_fadc->info
        fadc_ready = true;
        info_ready = true;
    }
    return *cache_fadc;
}

const ssp::SspEventData &EvChannel::Gem() const
{
    if (!cache_gem) cache_gem = std::make_unique<ssp::SspEventData>();
    if (!gem_ready) {
        decodeGemInto(*cache_gem);
        gem_ready = true;
    }
    return *cache_gem;
}

const tdc::TdcEventData &EvChannel::Tdc() const
{
    if (!cache_tdc) cache_tdc = std::make_unique<tdc::TdcEventData>();
    if (!tdc_ready) {
        decodeTdcInto(*cache_tdc);
        tdc_ready = true;
    }
    return *cache_tdc;
}

const vtp::VtpEventData &EvChannel::Vtp() const
{
    if (!cache_vtp) cache_vtp = std::make_unique<vtp::VtpEventData>();
    if (!vtp_ready) {
        decodeVtpInto(*cache_vtp);
        vtp_ready = true;
    }
    return *cache_vtp;
}

const psync::SyncInfo &EvChannel::Sync() const
{
    // Sync snapshot persists across events; only re-attempt decode once per
    // Scan().  On a SYNC/EPICS or control event the snapshot refreshes;
    // otherwise last_sync_info_ keeps the prior values.
    if (!sync_decoded_this_event_) {
        decodeSyncInto(last_sync_info_);
        sync_decoded_this_event_ = true;
    }
    return last_sync_info_;
}

const dsc::DscEventData &EvChannel::Dsc() const
{
    if (!cache_dsc) cache_dsc = std::make_unique<dsc::DscEventData>();
    if (!dsc_ready) {
        decodeDscInto(*cache_dsc);
        dsc_ready = true;
    }
    return *cache_dsc;
}

const epics::EpicsRecord &EvChannel::Epics() const
{
    if (!cache_epics) cache_epics = std::make_unique<epics::EpicsRecord>();
    if (!epics_ready) {
        decodeEpicsInto(*cache_epics);
        epics_ready = true;
    }
    return *cache_epics;
}

const uint8_t *EvChannel::GetCompositePayload(const EvNode &n, size_t &nbytes) const
{
    nbytes = 0;
    if (n.type != DATA_COMPOSITE || n.child_count < 2) return nullptr;
    auto &inner = nodes[n.child_first + 1];
    if (inner.data_begin + inner.data_words > buffer.size()) return nullptr;
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
// Per-product dispatchers — shared by the lazy cache accessors and the legacy
// DecodeEvent compat wrapper.  Each walks tag_index for the tags belonging to
// its product, then invokes the registered decoder (looked up by module name
// from DaqConfig::data_banks).
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
//    +-- [0x008E BANK]         Tagger crate
//    |    +-- [0xE10A UINT32]  TI data
//    |    +-- [0xE107 UINT32]  V1190 TDC hits
//    +-- [0x0031 BANK]         GEM VTP/MPD crate (when present)
//    |    +-- [0xE10A UINT32]  TI data
//    |    +-- [0x0DE9 UINT32]  MPD strip data (SSP bitfield format)
//    +-- ...
// =============================================================================

void EvChannel::decodeInfoInto(fdec::EventInfo &info) const
{
    info = fdec::EventInfo{};
    info.clear();
    if (nodes.empty()) return;

    BankHeader evh(&buffer[0]);
    info.event_tag = evh.tag;
    info.type = static_cast<uint8_t>(evtype);

    // event_tag = physics_base + trigger_type (see docs/rols/banktags.md).
    if (config.is_physics(evh.tag) && evh.tag >= config.physics_base)
        info.trigger_type = static_cast<uint8_t>(evh.tag - config.physics_base);

    // 0xC000 trigger bank → event number.
    {
        auto tb = tag_index.find(config.trigger_bank_tag);
        if (tb != tag_index.end() && !tb->second.empty())
            decodeTriggerInfo(nodes[tb->second[0]], info);
    }

    // 0xE10A TI banks — the first supplies trigger#/timestamp; any bank with
    // enough words also yields FP trigger_bits (TI master's 7-word variant).
    {
        auto ti = tag_index.find(config.ti_bank_tag);
        if (ti != tag_index.end()) {
            bool have_info = false;
            for (int ni : ti->second) {
                auto &n = nodes[ni];
                if (n.type != DATA_UINT32) continue;
                if (!have_info) { decodeTIBank(n, info, false); have_info = true; }
                if (info.trigger_bits == 0 && config.ti_trigger_type_word >= 0 &&
                    static_cast<size_t>(config.ti_trigger_type_word) < n.data_words)
                {
                    const uint32_t *d = GetData(n);
                    info.trigger_bits = (d[config.ti_trigger_type_word]
                                         >> config.ti_trigger_type_shift)
                                        & config.ti_trigger_type_mask;
                }
            }
        }
    }

    // 0xE10F run info lives inside the TI master crate (0x27).
    {
        auto tm = tag_index.find(config.ti_master_tag);
        if (tm != tag_index.end()) {
            for (int ni : tm->second) {
                auto &n = nodes[ni];
                if (n.depth != 1) continue;
                for (size_t ci = 0; ci < n.child_count; ++ci) {
                    auto &child = nodes[n.child_first + ci];
                    if (child.tag == config.run_info_tag && child.type == DATA_UINT32)
                        decodeRunInfo(child, info);
                }
                break;
            }
        }
    }

    // Track the latest physics event_number AND TI timestamp so slow
    // events (DSC2 SYNCs, EPICS) can stamp themselves with both — the
    // event_number for offline join-by-key, the timestamp so analysis
    // can place the slow row on the run timeline without depending on
    // the events tree containing that physics event (a downstream
    // trigger/error filter may have dropped it before write).
    if (evtype == EventType::Physics && info.event_number != 0) {
        last_physics_event_number_ = info.event_number;
        last_physics_timestamp_    = info.timestamp;
    }
}

void EvChannel::decodeFadcInto(fdec::EventData &evt) const
{
    evt.clear();
    decodeInfoInto(evt.info);

    int roc_idx = 0;
    for (uint32_t tag : fadc_tags) {
        auto it = tag_index.find(tag);
        if (it == tag_index.end()) continue;
        auto *bank_info = config.find_data_bank(tag);
        if (!bank_info) continue;
        const std::string &mod = bank_info->module;

        for (int ni : it->second) {
            if (roc_idx >= fdec::MAX_ROCS) break;
            auto &n = nodes[ni];
            if (n.data_words == 0) continue;
            if (n.parent >= 0 && nodes[n.parent].type == DATA_COMPOSITE) continue;
            uint32_t roc_tag = (n.parent >= 0) ? nodes[n.parent].tag : 0;

            if (mod == "fadc250_composite" && n.type == DATA_COMPOSITE) {
                size_t nbytes;
                auto *payload = GetCompositePayload(n, nbytes);
                if (!payload) continue;
                fdec::RocData &rd = evt.rocs[roc_idx];
                rd.present = true;
                rd.tag = roc_tag;
                fdec::Fadc250Decoder::DecodeRoc(payload, nbytes, rd);
                evt.roc_index[roc_idx] = roc_idx;
                roc_idx++;
            }
            else if (mod == "fadc250_raw" && n.type == DATA_UINT32) {
                if (n.data_begin + n.data_words > buffer.size()) continue;
                fdec::RocData &rd = evt.rocs[roc_idx];
                rd.present = true;
                rd.tag = roc_tag;
                fdec::Fadc250RawDecoder::DecodeRoc(GetData(n), n.data_words, rd);
                evt.roc_index[roc_idx] = roc_idx;
                roc_idx++;
            }
            else if (mod == "adc1881m" && config.adc_format == "adc1881m"
                     && n.type == DATA_UINT32)
            {
                if (n.data_begin + n.data_words > buffer.size()) continue;
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
            }
        }
    }
    evt.nrocs = roc_idx;
}

int EvChannel::decodeGemInto(ssp::SspEventData &ssp_evt) const
{
    ssp_evt.clear();
    int total_apvs = 0;
    for (uint32_t tag : gem_tags) {
        auto it = tag_index.find(tag);
        if (it == tag_index.end()) continue;
        for (int ni : it->second) {
            auto &n = nodes[ni];
            if (n.type != DATA_UINT32 || n.data_words == 0) continue;
            if (n.parent >= 0 && nodes[n.parent].type == DATA_COMPOSITE) continue;
            uint32_t roc_tag = (n.parent >= 0) ? nodes[n.parent].tag : 0;
            int crate_id = -1;
            for (auto &re : config.roc_tags)
                if (re.tag == roc_tag) { crate_id = re.crate; break; }
            // SspDecoder::DecodeRoc is safe for short banks (returns 0 APVs),
            // so the 3-word stub 0xE10C in the TI master crate is a no-op.
            int napvs = ssp::SspDecoder::DecodeRoc(GetData(n), n.data_words,
                                                    crate_id, ssp_evt);
            if (napvs > 0) total_apvs += napvs;
        }
    }
    return total_apvs;
}

void EvChannel::decodeTdcInto(tdc::TdcEventData &tdc_evt) const
{
    tdc_evt.clear();
    for (uint32_t tag : tdc_tags) {
        if (tag == 0) continue;
        auto it = tag_index.find(tag);
        if (it == tag_index.end()) continue;
        for (int ni : it->second) {
            auto &n = nodes[ni];
            if (n.type != DATA_UINT32 || n.data_words == 0) continue;
            if (n.parent >= 0 && nodes[n.parent].type == DATA_COMPOSITE) continue;
            uint32_t roc_tag = (n.parent >= 0) ? nodes[n.parent].tag : 0;
            tdc::TdcDecoder::DecodeRoc(GetData(n), n.data_words, roc_tag, tdc_evt);
        }
    }
}

void EvChannel::decodeVtpInto(vtp::VtpEventData &vtp_evt) const
{
    vtp_evt.clear();
    for (uint32_t tag : vtp_tags) {
        auto it = tag_index.find(tag);
        if (it == tag_index.end()) continue;
        for (int ni : it->second) {
            auto &n = nodes[ni];
            if (n.type != DATA_UINT32 || n.data_words == 0) continue;
            if (n.parent >= 0 && nodes[n.parent].type == DATA_COMPOSITE) continue;
            uint32_t roc_tag = (n.parent >= 0) ? nodes[n.parent].tag : 0;
            // VtpDecoder::DecodeRoc is tolerant of short stub banks.
            vtp::VtpDecoder::DecodeRoc(GetData(n), n.data_words, roc_tag, vtp_evt);
        }
    }
}

// Reads the 0xE112 HEAD bank (SYNC/EPICS events) or, failing that, the first
// UINT32 child of a control event (PRESTART/GO/END).  `out` is updated in
// place only when a source bank is found — physics events leave the prior
// snapshot intact, which is how Sync() provides persistent run-level state.
bool EvChannel::decodeSyncInto(psync::SyncInfo &out) const
{
    if (nodes.empty()) return false;

    auto read_word = [](const uint32_t *d, size_t nw, int off) -> uint32_t {
        return (off >= 0 && static_cast<size_t>(off) < nw) ? d[off] : 0;
    };

    // SYNC / EPICS path: 0xE112 HEAD bank somewhere in the tree.  Updates the
    // fields the HEAD bank actually carries; leaves run_type alone so a
    // prior PRESTART's run_type stays visible across intervening SYNCs.
    {
        auto it = tag_index.find(config.sync_head_tag);
        if (it != tag_index.end()) {
            for (int ni : it->second) {
                const auto &n = nodes[ni];
                if (n.type != DATA_UINT32 || n.data_words == 0) continue;
                const uint32_t *d = GetData(n);
                size_t nw = n.data_words;
                out.run_number   = read_word(d, nw, config.sync_head_run_number_word);
                out.sync_counter = read_word(d, nw, config.sync_head_counter_word);
                out.unix_time    = read_word(d, nw, config.sync_head_unix_time_word);
                out.event_tag    = read_word(d, nw, config.sync_head_event_tag_word);
                // Fall back to the wrapping event bank's tag if the declared
                // event_tag word is blank, so callers can always distinguish
                // PRESTART vs EPICS vs SYNC by event_tag alone.
                if (out.event_tag == 0)
                    out.event_tag = BankHeader(&buffer[0]).tag;
                return true;
            }
        }
    }

    // Control-event path: PRESTART/GO/END carry a 3-word UINT32 payload as
    // the first leaf child of the top-level event bank (same tag as the
    // event in CODA2 convention).  sync_counter is left untouched — it only
    // makes sense for 0xE112 HEAD banks, and preserving the last SYNC's
    // value lets callers still detect "new SYNC arrived" after a control
    // event by diffing counters across events.
    if (evtype == EventType::Prestart || evtype == EventType::Go ||
        evtype == EventType::End)
    {
        const auto &evn = nodes[0];
        for (size_t ci = 0; ci < evn.child_count; ++ci) {
            const auto &child = nodes[evn.child_first + ci];
            if (child.type != DATA_UINT32 || child.data_words == 0) continue;
            const uint32_t *d = GetData(child);
            size_t nw = child.data_words;
            out.unix_time  = read_word(d, nw, config.sync_control_unix_time_word);
            out.run_number = read_word(d, nw, config.sync_control_run_number_word);
            out.run_type   = static_cast<uint8_t>(
                read_word(d, nw, config.sync_control_run_type_word));
            out.event_tag  = evn.tag;          // 0x11 / 0x12 / 0x14
            return true;
        }
    }

    return false;
}

// === DSC2 livetime / EPICS slow-control =====================================

void EvChannel::decodeDscInto(dsc::DscEventData &out) const
{
    out.clear();
    const auto &cfg = config.dsc_scaler;
    if (!cfg.enabled() || cfg.bank_tag < 0) return;

    auto it = tag_index.find(static_cast<uint32_t>(cfg.bank_tag));
    if (it == tag_index.end()) return;
    for (int ni : it->second) {
        const auto &n = nodes[ni];
        if (n.data_words == 0) continue;
        // First bank that matches the configured slot wins.  DecodeBank
        // leaves out.present false on mismatch, so we keep probing.
        if (dsc::Dsc2Decoder::DecodeBank(GetData(n), n.data_words, cfg, out))
            return;
    }
}

void EvChannel::decodeEpicsInto(epics::EpicsRecord &out) const
{
    out.clear();
    if (evtype != EventType::Epics) return;

    const std::string text = ExtractEpicsText();
    if (text.empty()) return;
    epics::ParseEpicsText(text, out);

    // Stamp absolute time / run from the 0xE112 HEAD bank that this EPICS
    // event also carries (Sync() pulls it out and caches it).
    const auto &si = Sync();
    out.unix_time    = si.unix_time;
    out.sync_counter = si.sync_counter;
    out.run_number   = si.run_number;

    // last_physics_{event_number,timestamp}_ are updated by decodeInfoInto()
    // whenever a physics event is decoded.  Stamping both here lets the
    // EPICS row carry its own anchor on the run timeline, so analysis
    // does not have to look the timestamp back up via the events tree
    // (which may not include the anchored physics event if the replay
    // filtered it out before writing).
    out.event_number_at_arrival = last_physics_event_number_;
    out.timestamp_at_arrival    = last_physics_timestamp_;
    out.present = true;
}

// =============================================================================
// Legacy compat wrapper — writes directly into caller-owned structs, bypassing
// the lazy cache.  Semantics identical to the pre-refactor DecodeEvent so
// existing consumers compile and behave unchanged.  New code should prefer
// SelectEvent() + Info()/Fadc()/Gem()/Tdc()/Vtp().
// =============================================================================

bool EvChannel::DecodeEvent(int i, fdec::EventData &evt,
                            ssp::SspEventData *ssp_evt,
                            vtp::VtpEventData *vtp_evt,
                            tdc::TdcEventData *tdc_evt) const
{
    evt.clear();
    if (ssp_evt) ssp_evt->clear();
    if (vtp_evt) vtp_evt->clear();
    if (tdc_evt) tdc_evt->clear();
    if (i < 0 || i >= nevents) return false;
    if (nodes.empty()) return false;

    SelectEvent(i);
    decodeFadcInto(evt);                          // also fills evt.info

    bool ssp_decoded = false;
    if (ssp_evt)
        ssp_decoded = (decodeGemInto(*ssp_evt) > 0);
    if (vtp_evt) decodeVtpInto(*vtp_evt);
    if (tdc_evt) decodeTdcInto(*tdc_evt);

    // Same return convention as the original DecodeEvent: true iff any
    // detector data was decoded (FADC waveforms or GEM strips).  evt.info
    // is always populated regardless.
    return evt.nrocs > 0 || ssp_decoded;
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

// === DAQ config text extraction =============================================
//
// PRESTART events carry a 0xE10E STRING bank with the full concatenated DAQ
// configuration file the run was started with (TI / DSC / FADC / TDC / SSP /
// VTP / TS settings, pedestals & gains, trigger masks).  We pull the text by
// configured tag and strip null padding, mirroring ExtractEpicsText().

std::string EvChannel::ExtractDaqConfigText() const
{
    auto cfg_nodes = FindByTag(config.daq_config_tag);
    if (cfg_nodes.empty()) return {};

    // First STRING-typed bank wins (0xE10E is declared STRING in the
    // bank-structure JSON, but tolerate CHAR8 in case a future format
    // change relabels the type).
    const EvNode *picked = nullptr;
    for (const EvNode *n : cfg_nodes) {
        if (n->type == DATA_CHARSTAR8 || n->type == DATA_CHAR8) {
            picked = n;
            break;
        }
    }
    if (!picked) picked = cfg_nodes.front();
    if (picked->data_words == 0) return {};

    const char *raw = reinterpret_cast<const char*>(&buffer[picked->data_begin]);
    size_t max_len = picked->data_words * sizeof(uint32_t);
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
