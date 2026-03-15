#include "EtChannel.h"
#include <iostream>
#include <cstring>

using namespace evc;


static inline status et_status(int code, bool verbose = false)
{
    switch (code) {
    case ET_OK:
        return status::success;
    default:
        if (verbose) {
            std::cerr << "EtChannel Error: " << et_wrap::get_error_str(code) << "\n";
        }
        return status::failure;
    }
}


EtChannel::EtChannel(size_t chunk_buf)
: EvChannel(0), et_id(nullptr), stat_id(ID_NULL), att_id(ID_NULL)
{
    // large enough chunk
    pe.resize(chunk_buf);
    sconf.set_cue(ET_STATION_CUE);
    sconf.set_user(ET_STATION_USER_MULTI);
    sconf.set_restore(ET_STATION_RESTORE_OUT);
    sconf.set_select(ET_STATION_SELECT_ALL);
    sconf.set_block(ET_STATION_NONBLOCKING);
    sconf.set_prescale(1);
}

// Connect a ET system and create a monitor station with pre-settings
status EtChannel::Connect(const std::string &ip, int port, const std::string &et_file)
{
    if (IsETOpen()) {
        std::cout << "EtChannel Warning: a ET system has already been connected.\n";
    }

    // open et system
    et_wrap::OpenConfig conf;
    // use a direct connection to the ET system
    conf.set_cast(ET_DIRECT);
    conf.set_host(ip.c_str());
    conf.set_serverport(port);

    char *fname = strdup(et_file.c_str());
    auto status = et_open(&et_id, fname, conf.configure().get());
    free(fname);
    return et_status(status, true);
}

// create a station and attach to it
status EtChannel::Open(const std::string &station)
{
    if (!IsETOpen()) {
        std::cout << "EtChannel Error: ET System is not opened, cannot open the channel." << std::endl;
        return status::failure;
    }

    // create a station first
    if (stat_id == ID_NULL) {
        char *sname = strdup(station.c_str());
        et_station_create(et_id, &stat_id, sname, sconf.configure().get());
        free(sname);
    }

    return et_status(et_station_attach(et_id, stat_id, &att_id), true);
}

// detach from the station
void EtChannel::Close()
{
    if (IsETOpen() && (att_id != ID_NULL)) {
        et_status(et_station_detach(et_id, att_id), true);
        att_id = ID_NULL;
        et_status(et_station_remove(et_id, stat_id), true);
        stat_id = ID_NULL;
    }
}

void EtChannel::Disconnect()
{
    Close();
    if (IsETOpen()) {
        et_status(et_close(et_id), true);
        et_id = nullptr;
    }
}

// read an event
status EtChannel::Read()
{
    // buffers are not empty, just get the first element from it
    if (!buffers.empty()) {
        buffer = std::move(buffers.front());
        buffers.pop_front();
        return status::success;
    // read from ET system
    } else {
        int nread = 0;
        int chunk = sconf.get_cue();
        if (chunk > static_cast<int>(pe.size())) {
            pe.resize(chunk);
        }
        auto status = et_events_get(et_id, att_id, &pe[0], ET_ASYNC, nullptr, chunk, &nread);

        switch (status) {
        case ET_OK:
            if (!copyEvent(&pe[0], nread)) {
                return status::empty;
            }
            if (et_events_put(et_id, att_id, &pe[0], nread) != ET_OK) {
                std::cerr << "EtChannel Error: failed to put back et_event after reading.\n";
                return status::eof;
            }
            // get an event from the buffers
            buffer = std::move(buffers.front());
            buffers.pop_front();
            break;
        case ET_ERROR_BUSY:
        case ET_ERROR_TIMEOUT:
        case ET_ERROR_WAKEUP:
            std::cout << "EtChannel Warning: " << et_wrap::get_error_str(status) << "\n";
        case ET_ERROR_EMPTY:
            return status::empty;
        // fatal errors
        default:
            std::cerr << "EtChannel Error: " << et_wrap::get_error_str(status) << "\n";
            return status::failure;
        }
    }
    return status::success;
}

// helper functions
template<class Func, class... Args>
inline bool ev_filter(Func beg, Func end, Args&&... args)
{
    for (auto it = beg; it != end; ++it) {
        if (!(*it)(args...)) return false;
    }
    return true;
}

template<class Func>
std::vector<uint32_t> copy_event(const uint32_t *buf, bool swap, Func fil_beg, Func fil_end)
{
    std::vector<uint32_t> event;
    auto header = BankHeader(buf);

    // invalid header
    if (header.length < 1) {
        return event;
    }

    // cannot pass filters
    for (auto it = fil_beg; it != fil_end; ++it) {
        if (!(*it)(header)) {
            return event;
        }
    }

    // good event, copy it!
    event.assign(buf, buf + header.length + 1);
    if (swap) {
        for (auto &val : event) {
            val = ET_SWAP32(val);
        }
    }

    return event;
}

// copy event to the buffer
bool EtChannel::copyEvent(et_event **pe, int nread)
{
    if (nread <= 0) {
        return false;
    }

    void *data;
    size_t len, bytes = sizeof(uint32_t);
    int endian, swap;

    for (int i = 0; i < nread; ++i) {
        // get event data and attributes from ET
        et_event_getdata(pe[i], &data);
        et_event_getlength(pe[i], &len);
        // et_event_getendian(pe[i], &endian);
        et_event_needtoswap(pe[i], &swap);

        // size of the buffer
        len = len / bytes + ((len % bytes) ? 1 : 0);
        uint32_t *dbuf = static_cast<uint32_t*>(data);

        // check if it is a block
        if ((len > 6) && (0xc0da0100 == (swap ? ET_SWAP32(dbuf[7]) : dbuf[7]))) {
            // skip the block header (size 8)
            uint32_t index = 8;
            while (index < dbuf[0]) {
                auto event = copy_event(&dbuf[index], swap, filters.begin(), filters.end());
                index += dbuf[index] + 1;
                if (event.size()) {
                    buffers.emplace_back(std::move(event));
                }
            }
        // a single event
        } else {
            auto event = copy_event(dbuf, swap, filters.begin(), filters.end());
            if (event.size()) {
                buffers.emplace_back(std::move(event));
            }
        }
    }

    return true;
}

