#pragma once
//=============================================================================
// EvChannel.h — read evio events, scan bank tree, lazily decode per product
//
// New (lazy) API:
//   EvChannel ch;
//   ch.SetConfig(cfg);             // builds per-product tag lists
//   ch.Open("file.evio");
//   while (ch.Read() == status::success) {
//       if (!ch.Scan()) continue;
//       if (ch.GetEventType() != EventType::Physics) continue;
//       for (int i = 0; i < ch.GetNEvents(); ++i) {
//           ch.SelectEvent(i);
//           const auto &info = ch.Info();   // always cheap
//           const auto &fadc = ch.Fadc();   // decoded on first call
//           // multiple requests for the same product reuse the cached result
//       }
//   }
//
// Legacy DecodeEvent is preserved as a thin compat wrapper so existing
// callers compile and behave unchanged.  The earlier DecodeEventInfo /
// DecodeEventTdc fast paths were removed in favour of SelectEvent() + the
// corresponding Info()/Tdc() accessors, which decode only what's requested
// and cache the result for repeat calls on the same event.
//=============================================================================

#include "EvStruct.h"
#include "Fadc250Data.h"
#include "SspData.h"
#include "VtpData.h"
#include "TdcData.h"
#include "DaqConfig.h"
#include <string>
#include <vector>
#include <unordered_map>
#include <memory>

namespace evc {

enum class status : int { failure = -1, success = 1, incomplete = 2, empty = 3, eof = 4 };

class EvChannel
{
public:
    EvChannel(size_t buflen = 1024 * 2000);
    virtual ~EvChannel() { Close(); }
    EvChannel(const EvChannel &) = delete;
    EvChannel &operator=(const EvChannel &) = delete;

    // --- configuration ------------------------------------------------------
    // Stores the config and precomputes per-product tag lists used by the
    // lazy accessors.  If the config's data_banks map is empty (legacy JSON
    // without a bank_structure section), default entries are synthesised
    // from the legacy bank-tag fields so older configs keep working.
    void SetConfig(const DaqConfig &cfg);
    const DaqConfig &GetConfig() const   { return config; }

    virtual status Open(const std::string &path);
    virtual void   Close();
    virtual status Read();

    // --- scan the current event into a flat tree ----------------------------
    // Rebuilds nodes[] and the tag index; invalidates the per-product cache.
    bool Scan();

    // --- event type (valid after Scan) --------------------------------------
    EventType GetEventType() const { return evtype; }

    // --- tree accessors -----------------------------------------------------
    BankHeader                  GetEvHeader() const { return BankHeader(&buffer[0]); }
    const std::vector<EvNode>  &GetNodes()    const { return nodes; }
    const EvNode               &GetChild(const EvNode &n, size_t i) const { return nodes[n.child_first + i]; }
    std::vector<const EvNode*>  FindByTag(uint32_t tag) const;

    // Find first node with given tag (no allocation). Returns nullptr if not found.
    const EvNode *FindFirstByTag(uint32_t tag) const;

    // O(1) lookup of every node index carrying a given tag in the current
    // event (populated by Scan).  Empty span if the tag is not present.
    const std::vector<int> &NodesForTag(uint32_t tag) const;

    const uint32_t *GetData(const EvNode &n) const { return &buffer[n.data_begin]; }
    const uint8_t  *GetBytes(const EvNode &n) const
    { return reinterpret_cast<const uint8_t*>(&buffer[n.data_begin]); }
    size_t GetDataBytes(const EvNode &n) const { return n.data_words * sizeof(uint32_t); }
    const uint8_t *GetCompositePayload(const EvNode &n, size_t &nbytes) const;

    uint32_t       *GetRawBuffer()       { return buffer.data(); }
    const uint32_t *GetRawBuffer() const { return buffer.data(); }

    // Number of events in this block (from the physics event header num field).
    // For single-event mode this is 1. For multi-event blocks this is M.
    int GetNEvents() const { return nevents; }

    // --- lazy data-product accessors (new API) ------------------------------
    //
    // Choose the sub-event index subsequent Get*() calls refer to.  Clears
    // the product cache if the index changed; for PRad-II single-event data,
    // pass 0 (the default after Scan()).  Safe to call repeatedly.
    void SelectEvent(int i) const;

    // Each accessor decodes on first call after SelectEvent(), then returns
    // a cached reference.  References are invalidated by the next Read(),
    // Scan(), or SelectEvent() call.
    const fdec::EventInfo    &Info() const;  // always cheap
    const fdec::EventData    &Fadc() const;  // FADC250 + ADC1881M waveforms
    const ssp::SspEventData  &Gem()  const;  // SSP/MPD GEM strips
    const tdc::TdcEventData  &Tdc()  const;  // V1190 timing hits
    const vtp::VtpEventData  &Vtp()  const;  // VTP ECAL peaks/clusters

    // --- legacy API (compat, writes directly to caller-owned structs) -------
    //
    // Preserves the original semantics: populate the caller's structs without
    // touching the lazy cache.  Kept so existing consumers compile unchanged;
    // new code should use SelectEvent() + Info()/Fadc()/Gem()/Tdc()/Vtp().
    bool DecodeEvent(int i, fdec::EventData &evt,
                     ssp::SspEventData *ssp_evt = nullptr,
                     vtp::VtpEventData *vtp_evt = nullptr,
                     tdc::TdcEventData *tdc_evt = nullptr) const;

    // --- Control event extraction (Prestart/Go/End) -------------------------

    // Extract unix timestamp from PRESTART or GO event.
    // CODA2 format: data words [time, run_number, run_type]
    // Returns 0 if not a control event or no time found.
    uint32_t GetControlTime() const;

    // --- EPICS extraction (call when GetEventType() == Epics) ---------------

    // Extract raw EPICS text from the current event buffer.
    // Returns the text payload (lines of "value  channel_name").
    // Returns empty string if no EPICS data found.
    std::string ExtractEpicsText() const;

    // debug
    void PrintTree(std::ostream &os) const;

protected:
    DaqConfig config;
    int fHandle;
    std::vector<uint32_t> buffer;
    std::vector<EvNode>   nodes;
    int nevents = 0;
    EventType evtype = EventType::Unknown;

    // tag → every node index in the current event that carries it.
    // Rebuilt by Scan(); consulted by the lazy accessors to avoid re-scanning.
    std::unordered_map<uint32_t, std::vector<int>> tag_index;

    // Per-product tag lists derived from config.data_banks at SetConfig().
    std::vector<uint32_t> fadc_tags;
    std::vector<uint32_t> gem_tags;
    std::vector<uint32_t> tdc_tags;
    std::vector<uint32_t> vtp_tags;

    // --- product cache (populated by Info/Fadc/Gem/Tdc/Vtp) -----------------
    // Cleared on Read/Scan/SelectEvent(i != current).  Heap-allocated on first
    // use — each product struct is sized for the worst-case event (EventData
    // alone is ~8 MB), so keeping them inline would blow the stack when
    // callers declare `EvChannel ch;` as a local.  Marked mutable so the
    // accessors stay const-callable — the cache is an implementation detail.
    mutable int  cached_event_idx = -1;
    mutable bool info_ready = false;
    mutable bool fadc_ready = false;
    mutable bool gem_ready  = false;
    mutable bool tdc_ready  = false;
    mutable bool vtp_ready  = false;
    mutable std::unique_ptr<fdec::EventData>   cache_fadc;   // .info also serves Info()
    mutable std::unique_ptr<ssp::SspEventData> cache_gem;
    mutable std::unique_ptr<tdc::TdcEventData> cache_tdc;
    mutable std::unique_ptr<vtp::VtpEventData> cache_vtp;

    // --- per-bank decoders (shared by legacy and lazy paths) ----------------
    void decodeTriggerInfo(const EvNode &node, fdec::EventInfo &info) const;
    void decodeTIBank(const EvNode &node, fdec::EventInfo &info, bool is_master) const;
    void decodeRunInfo(const EvNode &node, fdec::EventInfo &info) const;

    // --- per-product dispatchers (write into caller-supplied structs) -------
    void decodeInfoInto(fdec::EventInfo &info) const;
    void decodeFadcInto(fdec::EventData &evt) const;    // fills evt.info too
    // Returns the total APV count across every SSP/MPD bank decoded — used by
    // the legacy DecodeEvent compat wrapper to preserve "true iff data found".
    int  decodeGemInto (ssp::SspEventData &ssp) const;
    void decodeTdcInto (tdc::TdcEventData &tdc) const;
    void decodeVtpInto (vtp::VtpEventData &vtp) const;

    // Invalidate all product cache flags.
    void clearCache() const;

    size_t scanBank      (size_t off, int depth, int parent);
    size_t scanSegment   (size_t off, int depth, int parent);
    size_t scanTagSegment(size_t off, int depth, int parent);
    void   scanChildren  (size_t off, size_t nwords, uint32_t ptype, int depth, int pidx);
};

} // namespace evc
