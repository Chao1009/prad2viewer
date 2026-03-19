#pragma once
//=============================================================================
// EvChannel.h — read evio events, scan bank tree, decode FADC data
//
// Usage:
//   DaqConfig cfg;                           // default tags, or load from JSON
//   EvChannel ch;
//   ch.SetConfig(cfg);
//   ch.Open("file.evio");
//   while (ch.Read() == status::success) {
//       if (!ch.Scan()) continue;
//       auto etype = ch.GetEventType();
//
//       if (etype == EventType::Physics) {
//           int nevt = ch.GetNEvents();
//           for (int i = 0; i < nevt; ++i) {
//               ch.DecodeEvent(i, event);
//               // event.info has timestamp, trigger number, event type
//               // event.rocs[r].slots[s].channels[c].samples[]
//           }
//       }
//       else if (etype == EventType::Epics) {
//           std::string text = ch.ExtractEpicsText();
//           // parse "value  channel_name" lines
//       }
//   }
//=============================================================================

#include "EvStruct.h"
#include "Fadc250Data.h"
#include "DaqConfig.h"
#include <string>
#include <vector>

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
    void SetConfig(const DaqConfig &cfg) { config = cfg; }
    const DaqConfig &GetConfig() const   { return config; }

    virtual status Open(const std::string &path);
    virtual void   Close();
    virtual status Read();

    // --- scan the current event into a flat tree ----------------------------
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

    const uint32_t *GetData(const EvNode &n) const { return &buffer[n.data_begin]; }
    const uint8_t  *GetBytes(const EvNode &n) const
    { return reinterpret_cast<const uint8_t*>(&buffer[n.data_begin]); }
    size_t GetDataBytes(const EvNode &n) const { return n.data_words * sizeof(uint32_t); }
    const uint8_t *GetCompositePayload(const EvNode &n, size_t &nbytes) const;

    uint32_t       *GetRawBuffer()       { return buffer.data(); }
    const uint32_t *GetRawBuffer() const { return buffer.data(); }

    // --- event-by-event access (after Scan) ---------------------------------

    // Number of events in this block (from the physics event header num field).
    // For single-event mode this is 1. For multi-event blocks this is M.
    int GetNEvents() const { return nevents; }

    // Decode the i-th event (0-based) into the pre-allocated EventData.
    // Populates EventInfo (type, trigger, timestamp) and FADC ROC data.
    // Returns true on success.
    bool DecodeEvent(int i, fdec::EventData &evt) const;

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

    // CODA built trigger bank decoding (spec pages 21, 26, 31)
    // Extracts event number, timestamps, trigger type, run info from 0xFF2X bank.
    bool decodeTriggerBank(int event_idx, fdec::EventInfo &info) const;

    // JLab TI bank fallback (0xE10A) — used when no CODA trigger bank is present
    bool decodeTI(fdec::EventInfo &info) const;

    size_t scanBank      (size_t off, int depth, int parent);
    size_t scanSegment   (size_t off, int depth, int parent);
    size_t scanTagSegment(size_t off, int depth, int parent);
    void   scanChildren  (size_t off, size_t nwords, uint32_t ptype, int depth, int pidx);
};

} // namespace evc
