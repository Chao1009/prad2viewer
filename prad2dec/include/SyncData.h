#pragma once
//=============================================================================
// SyncData.h — absolute-time / run-state snapshot from SYNC and control events
//
// In PRad-II the CODA DAQ writes one of two bank flavors depending on the
// event type.  Both carry a 32-bit absolute unix timestamp that the physics
// stream's 48-bit TI tick counter doesn't:
//
//   * SYNC/EPICS (event tag 0x001F) — 0xE112 "HEAD" bank, 5 UINT32 words:
//         d[0] : reserved (observed 0)
//         d[1] : run number
//         d[2] : monotonic sync counter (increments by 1 per SYNC event)
//         d[3] : absolute unix time (seconds since epoch)
//         d[4] : wrapping event tag (0x1F)
//
//   * PRESTART / GO / END (event tags 0x0011 / 0x0012 / 0x0014) — first
//     child UINT32 bank, 3 words:
//         d[0] : absolute unix time
//         d[1] : run number
//         d[2] : run type (PRESTART only; 0 for GO/END)
//
// `EvChannel::Sync()` returns the latest SyncInfo observed since the channel
// was opened.  It auto-updates whenever `Scan()` lands on a SYNC/EPICS or
// control event and otherwise returns the cached snapshot — so a physics-
// event loop always has an absolute time to anchor its 48-bit TI tick delta
// against.  To detect "a new SYNC just arrived", diff `sync_counter` against
// your last-seen value: it's monotonic for 0xE112 and stays 0 for control
// events (distinguish those via `event_tag`).
//=============================================================================

#include <cstdint>

namespace sync
{

struct SyncInfo {
    uint32_t run_number   = 0;   // 0xE112 d[1] / control d[1]
    uint32_t sync_counter = 0;   // 0xE112 d[2]  (0 for control events)
    uint32_t unix_time    = 0;   // 0xE112 d[3] / control d[0]
    uint32_t event_tag    = 0;   // wrapping event tag (0x1F / 0x11 / 0x12 / 0x14)
    uint8_t  run_type     = 0;   // control-event d[2]  (0 for 0xE112)

    void clear() {
        run_number = sync_counter = unix_time = event_tag = 0;
        run_type = 0;
    }

    bool valid() const { return unix_time != 0; }
};

} // namespace sync
