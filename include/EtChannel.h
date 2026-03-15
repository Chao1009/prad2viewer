//=============================================================================
// Class EtChannel                                                           ||
// Read event from ET instead of CODA file                                   ||
//                                                                           ||
// Developer:                                                                ||
// Chao Peng                                                                 ||
// 11/22/2019                                                                ||
//=============================================================================
#pragma once

#include "EtConfigWrapper.h"
#include "EvChannel.h"
#include "EvStruct.h"
#include <functional>
#include <iostream>
#include <chrono>
#include <string>
#include <thread>
#include <list>
#include <vector>


#define ID_NULL -9999

namespace evc {

class EtChannel : public EvChannel
{
public:
    EtChannel(size_t chunk_buf = 2000);
    virtual ~EtChannel() { Disconnect(); }

    EtChannel(const EtChannel &)  = delete;
    void operator =(const EtChannel &)  = delete;

    virtual status Open(const std::string &station);
    virtual void Close();
    virtual status Read();

    status Connect(const std::string &ip, int port, const std::string &et_file);
    void Disconnect();
    bool IsETOpen() const { return (et_id != nullptr) && et_alive(et_id); }
    void AddEvFilter(std::function<bool(const BankHeader &)> &&func) { filters.emplace_back(func); }

    et_wrap::StationConfig &GetConfig() { return sconf; }
    const et_wrap::StationConfig &GetConfig() const { return sconf; }

private:
    bool copyEvent(et_event **pe, int nread);

    et_wrap::StationConfig sconf;
    et_sys_id et_id;
    et_stat_id stat_id;
    et_att_id att_id;
    std::list<std::vector<uint32_t>> buffers;
    std::vector<std::function<bool(const BankHeader &)>> filters;
    std::vector<et_event*> pe;
};

}   // namespace evc
