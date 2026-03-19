#pragma once
//=============================================================================
// EpicsStore.h — accumulate and index EPICS slow control snapshots
//
// Fed by the decoder's ExtractEpicsText() output. Dynamically discovers
// channels. Provides lookup of the most recent snapshot for any event number.
//
// Usage:
//   EpicsStore epics;
//   // in event loop:
//   if (ch.GetEventType() == EventType::Epics) {
//       epics.Feed(event_number, timestamp, ch.ExtractEpicsText());
//   }
//   // later, for any physics event:
//   float beam_current;
//   if (epics.GetValue(event_number, "beam_current", beam_current)) { ... }
//=============================================================================

#include <string>
#include <vector>
#include <unordered_map>
#include <cstdint>

namespace fdec
{

class EpicsStore
{
public:
    EpicsStore() = default;

    // --- feeding data -------------------------------------------------------

    // Parse raw EPICS text and store a snapshot.
    // Text format: one "value  channel_name" pair per line.
    // event_number: the trigger/event number at the time of this EPICS update.
    // timestamp: 48-bit TI timestamp (0 if unavailable).
    void Feed(int32_t event_number, uint64_t timestamp, const std::string &text);

    // --- querying -----------------------------------------------------------

    // Get the most recent value of a channel at or before the given event number.
    // Returns true if found, false if channel unknown or no snapshot before this event.
    bool GetValue(int32_t event_number, const std::string &channel, float &value) const;

    // Get all channel values from the most recent snapshot at or before event_number.
    // Returns pointer to the values array (indexed by channel id), or nullptr if none.
    // Use GetChannelId() to map names to indices.
    struct Snapshot {
        int32_t              event_number;
        uint64_t             timestamp;
        std::vector<float>   values;     // indexed by channel id
    };

    const Snapshot *FindSnapshot(int32_t event_number) const;

    // --- channel info -------------------------------------------------------

    int  GetChannelCount() const { return static_cast<int>(channel_names_.size()); }
    int  GetChannelId(const std::string &name) const;
    const std::string &GetChannelName(int id) const { return channel_names_[id]; }

    // all known channel names
    const std::vector<std::string> &GetChannelNames() const { return channel_names_; }

    // number of snapshots stored
    int  GetSnapshotCount() const { return static_cast<int>(snapshots_.size()); }

    // --- reset --------------------------------------------------------------
    void Clear();

private:
    int get_or_create_channel(const std::string &name);

    std::vector<std::string>                   channel_names_;
    std::unordered_map<std::string, int>       channel_map_;    // name → id
    std::vector<Snapshot>                      snapshots_;      // sorted by event_number
};

} // namespace fdec
