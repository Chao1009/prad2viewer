//=============================================================================
// EtConfigWrapper                                                           ||
// A C++ wrapper for the configs of ET system/station                        ||
//                                                                           ||
// Developer:                                                                ||
// Chao Peng                                                                 ||
// 11/22/2019                                                                ||
//=============================================================================

#pragma once

#include <vector>
#include <memory>
#include <string>
#include <iostream>
#include <unordered_set>
#include "et.h"

#define ET_VERSION 16

#define SET_BIT(n, i) ( (n) |= (1ULL << i) )
#define CLEAR_BIT(n, i) ( (n) &= ~(1ULL << i) )
#define TEST_BIT(n, i) ( (bool)( n & (1ULL << i) ) )

#define ETCONF_ADD_MEMBER(type, var, flag) \
    public:\
        void set_##var(type val) { var = val; SET_BIT(flag, static_cast<uint32_t>(Flag::var)); } \
        type get_##var() const { return var; } \
    private:\
        type var;

namespace et_wrap {

static std::string get_error_str(int error)
{
    switch(error) {
    case ET_ERROR: return "General error.";
    case ET_ERROR_TOOMANY: return "Too many somethings (stations, attachments, temp events, ET system responses) exist.";
    case ET_ERROR_EXISTS: return "ET system file or station already exists.";
    case ET_ERROR_WAKEUP: return "Sleeping routine woken up by et_wakeup_attachment() or et_wakeup_all().";
    case ET_ERROR_TIMEOUT: return "Timed out.";
    case ET_ERROR_EMPTY: return "No events available in async mode.";
    case ET_ERROR_BUSY: return "Resource is busy.";
    case ET_ERROR_DEAD: return "ET system is dead.";
    case ET_ERROR_READ: return "Network read error.";
    case ET_ERROR_WRITE: return "Network write error,";
    case ET_ERROR_REMOTE: return "Cannot allocate memory in remote client.";
#if ET_VERSION >= 14
    case ET_ERROR_TOOBIG: return "Client is 32 bits & server is 64 (or vice versa) and event is too big for one.";
    case ET_ERROR_NOMEM: return "Cannot allocate memory.";
    case ET_ERROR_BADARG: return "Bad argument given to function.";
    case ET_ERROR_SOCKET: return "Socket option could not be set.";
    case ET_ERROR_NETWORK: return "Host name or address could not be resolved, or cannot connect.";
    case ET_ERROR_CLOSED: return "ET system has been closed by client.";
    case ET_ERROR_JAVASYS: return "C code trying to open Java-based ET system file locally.";
#endif
    default: break;
    }
    return "Unknown error";
}

class OpenConfig
{
#define OPEN_CONFIG_SET(flag, ptr, var) \
    if (TEST_BIT(flag, static_cast<uint32_t>(Flag::var))) { et_open_config_set##var(ptr, get_##var()); }

public:
    enum class Flag : uint32_t {
        wait = 0,
        cast,
        TTL,
        mode,
        debugdefault,
        port,
        serverport,
        policy,
        tcp,
        timeout,
        host,
        interface
    };

    struct tcp_setting {
        int rbuf_size, sbuf_size, no_delay;
    };

public:
    // initialize
    OpenConfig() { flag = 0; }

    // set configuration and return a smart pointer
    std::shared_ptr<void> configure() const
    {
        void *ptr;
        et_open_config_init(&ptr);

        OPEN_CONFIG_SET(flag, ptr, wait);
        OPEN_CONFIG_SET(flag, ptr, cast);
        OPEN_CONFIG_SET(flag, ptr, mode);
        OPEN_CONFIG_SET(flag, ptr, port);
        OPEN_CONFIG_SET(flag, ptr, serverport);
        OPEN_CONFIG_SET(flag, ptr, timeout);

        if (TEST_BIT(flag, static_cast<uint32_t>(Flag::host))) {
            char temp[1024];
            strncpy(temp, host.c_str(), 1024);
            et_open_config_sethost(ptr, temp);
        }

#if ET_VERSION >= 14
        OPEN_CONFIG_SET(flag, ptr, debugdefault);
        OPEN_CONFIG_SET(flag, ptr, TTL);
        OPEN_CONFIG_SET(flag, ptr, policy);

        if (TEST_BIT(flag, static_cast<uint32_t>(Flag::interface))) {
            et_open_config_setinterface(ptr, interface.c_str());
        }
        for (auto &c : multi_casts) {
            et_open_config_addmulticast(ptr, c.c_str());
        }
        for (auto &c : broad_casts) {
            et_open_config_addbroadcast(ptr, c.c_str());
        }
        if (TEST_BIT(flag, static_cast<uint32_t>(Flag::tcp))) {
            et_open_config_settcp(ptr, tcp.rbuf_size, tcp.sbuf_size, tcp.no_delay);
        }
#endif

        return std::shared_ptr<void>(ptr, [] (void *p) { et_open_config_destroy(p); });
    }

    std::unordered_set<std::string> broad_casts, multi_casts;

private:
    uint32_t flag;

    ETCONF_ADD_MEMBER(int, wait, flag);
    ETCONF_ADD_MEMBER(int, cast, flag);
    ETCONF_ADD_MEMBER(int, TTL, flag);
    ETCONF_ADD_MEMBER(int, mode, flag);
    ETCONF_ADD_MEMBER(int, debugdefault, flag);
    ETCONF_ADD_MEMBER(int, port, flag);
    ETCONF_ADD_MEMBER(int, serverport, flag);
    ETCONF_ADD_MEMBER(int, policy, flag);
    ETCONF_ADD_MEMBER(struct tcp_setting, tcp, flag);
    ETCONF_ADD_MEMBER(struct timespec, timeout, flag);
    ETCONF_ADD_MEMBER(std::string, host, flag);
    ETCONF_ADD_MEMBER(std::string, interface, flag);
};

class StationConfig
{
#define STATION_CONFIG_SET(flag, ptr, var) \
    if (TEST_BIT(flag, static_cast<uint32_t>(Flag::var))) { et_station_config_set##var(ptr, get_##var()); }

public:
    enum class Flag : uint32_t {
        block = 0,
        flow,
        select,
        user,
        restore,
        cue,
        prescale,
        selectwords,
        function,
        lib,
        myclass,
    };

public:
    // initialize
    StationConfig() { flag = 0; }

    // set configuration and return a smart pointer
    std::shared_ptr<void> configure() const
    {
        void *ptr;
        et_station_config_init(&ptr);

        STATION_CONFIG_SET(flag, ptr, block);
        STATION_CONFIG_SET(flag, ptr, select);
        STATION_CONFIG_SET(flag, ptr, user);
        STATION_CONFIG_SET(flag, ptr, restore);
        STATION_CONFIG_SET(flag, ptr, cue);
        STATION_CONFIG_SET(flag, ptr, prescale);

        if (!selectwords.empty() && TEST_BIT(flag, static_cast<uint32_t>(Flag::selectwords))) {
            // copy a vector to maintain the const behavior
            auto words = selectwords;
            et_station_config_setselectwords(ptr, &words[0]);
        }

        if (TEST_BIT(flag, static_cast<uint32_t>(Flag::function))) {
            char temp[1024];
            strncpy(temp, function.c_str(), 1024);
            if (et_station_config_setfunction(ptr, temp) != ET_OK) {
                std::cerr << "Could not set function \"" << function << "\" for station config." << std::endl;
            }
        }
        if (TEST_BIT(flag, static_cast<uint32_t>(Flag::lib))) {
            char temp[1024];
            strncpy(temp, lib.c_str(), 1024);
            if (et_station_config_setlib(ptr, temp) != ET_OK) {
                std::cerr << "Could not set library \"" << lib << "\" for station config." << std::endl;
            }
        }

#if ET_VERSION >= 14
        STATION_CONFIG_SET(flag, ptr, flow);
        if (TEST_BIT(flag, static_cast<uint32_t>(Flag::myclass))) {
            if (et_station_config_setclass(ptr, myclass.c_str()) != ET_OK) {
                std::cerr << "Could not set class \"" << myclass << "\" for station config." << std::endl;
            }
        }
#endif

        return std::shared_ptr<void>(ptr, [] (void *p) { et_station_config_destroy(p); });
    }

    std::unordered_set<std::string> broad_casts, multi_casts;

private:
    uint32_t flag;

    ETCONF_ADD_MEMBER(int, block, flag);
    ETCONF_ADD_MEMBER(int, flow, flag);
    ETCONF_ADD_MEMBER(int, select, flag);
    ETCONF_ADD_MEMBER(int, user, flag);
    ETCONF_ADD_MEMBER(int, restore, flag);
    ETCONF_ADD_MEMBER(int, cue, flag);
    ETCONF_ADD_MEMBER(int, prescale, flag);
    ETCONF_ADD_MEMBER(std::vector<int>, selectwords, flag);
    ETCONF_ADD_MEMBER(std::string, function, flag);
    ETCONF_ADD_MEMBER(std::string, lib, flag);
    ETCONF_ADD_MEMBER(std::string, myclass, flag);
};

}; // namespace et_wrapper

