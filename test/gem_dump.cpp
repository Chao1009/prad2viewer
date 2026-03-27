// test/gem_dump.cpp — GEM data diagnostic tool
//
// Decodes SSP/MPD banks from EVIO data, optionally runs GEM reconstruction,
// and prints diagnostic output at each stage of the pipeline.
//
// Usage:
//   gem_dump <evio_file> [-D <daq_config.json>] [options]
//
// Modes (default: summary):
//   -m raw        Dump raw SSP-decoded APV data (strips × time samples)
//   -m hits       Process through GemSystem → show strip hits per plane
//   -m clusters   Full reconstruction → show clusters and 2D GEM hits
//   -m summary    Statistics: MPDs, APVs, strips, hits, clusters per event
//   -m evdump     Dump event(s) with 2D hits to JSON (-n K for first K)
//
// Options:
//   -D <file>     DAQ configuration (auto-searches daq_config.json if omitted)
//   -G <file>     GEM map file (default: gem_map.json from DAQ config dir)
//   -P <file>     GEM pedestal file (optional, required for good hit finding)
//   -n <N>        Max physics events to process (default: 10, 0=all)
//   -t <bit>      Trigger bit filter (default: -1 = accept all)
//   -e <N>        Dump only event N (1-based physics event number)

#include "EvChannel.h"
#include "Fadc250Data.h"
#include "SspData.h"
#include "load_daq_config.h"
#include "GemSystem.h"
#include "GemCluster.h"

#include <nlohmann/json.hpp>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <map>
#include <set>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <getopt.h>
#include <memory>

using namespace evc;

// -------------------------------------------------------------------------
// Helpers
// -------------------------------------------------------------------------
static std::string hex(uint32_t v)
{
    char buf[16];
    snprintf(buf, sizeof(buf), "0x%04X", v);
    return buf;
}

// -------------------------------------------------------------------------
// Event filter for evdump: field=val[,val]:min_dets
//   e.g.  pos=10,11:3   match=+Y,-Y:2   plane=X:4
// -------------------------------------------------------------------------
struct EvdumpFilter {
    std::string field;
    std::vector<std::string> values;
    int min_dets = 1;
};

static bool parseEvdumpFilter(const std::string &expr, EvdumpFilter &filt)
{
    // split at ':'
    auto colon = expr.rfind(':');
    std::string lhs = (colon != std::string::npos) ? expr.substr(0, colon) : expr;
    filt.min_dets = (colon != std::string::npos) ? std::atoi(expr.substr(colon + 1).c_str()) : 1;

    // split at '='
    auto eq = lhs.find('=');
    if (eq == std::string::npos) {
        std::cerr << "Error: filter must be field=val[,val]:min_dets\n";
        return false;
    }
    filt.field = lhs.substr(0, eq);
    std::string vals = lhs.substr(eq + 1);

    // split values by ','
    filt.values.clear();
    size_t start = 0;
    while (start < vals.size()) {
        auto comma = vals.find(',', start);
        if (comma == std::string::npos) comma = vals.size();
        filt.values.push_back(vals.substr(start, comma - start));
        start = comma + 1;
    }
    return !filt.values.empty();
}

static bool apvMatchesFilter(const gem::ApvConfig &cfg, const EvdumpFilter &filt)
{
    std::string val;
    if (filt.field == "pos")
        val = std::to_string(cfg.plane_index);
    else if (filt.field == "plane")
        val = (cfg.plane_type == 0) ? "X" : "Y";
    else if (filt.field == "match")
        val = cfg.match;
    else if (filt.field == "orient")
        val = std::to_string(cfg.orient);
    else if (filt.field == "det")
        val = std::to_string(cfg.det_id);
    else
        return false;

    for (auto &v : filt.values)
        if (v == val) return true;
    return false;
}

static bool checkEvdumpFilter(const EvdumpFilter &filt,
                               gem::GemSystem &sys,
                               const ssp::SspEventData &ssp)
{
    std::set<int> matching_dets;
    for (int m = 0; m < ssp.nmpds; ++m) {
        auto &mpd = ssp.mpds[m];
        if (!mpd.present) continue;
        for (int a = 0; a < ssp::MAX_APVS_PER_MPD; ++a) {
            if (!mpd.apvs[a].present) continue;
            int idx = sys.FindApvIndex(mpd.crate_id, mpd.mpd_id, a);
            if (idx < 0) continue;
            if (!sys.HasApvZsHits(idx)) continue;

            auto &cfg = sys.GetApvConfig(idx);
            if (apvMatchesFilter(cfg, filt))
                matching_dets.insert(cfg.det_id);
        }
    }
    return static_cast<int>(matching_dets.size()) >= filt.min_dets;
}

// -------------------------------------------------------------------------
// Mode: raw — dump decoded SSP APV data
// -------------------------------------------------------------------------
static void dumpRawSsp(const ssp::SspEventData &ssp, int phys_ev)
{
    std::cout << "--- Event " << phys_ev << ": "
              << ssp.nmpds << " MPD(s) ---\n";

    for (int m = 0; m < ssp.nmpds; ++m) {
        auto &mpd = ssp.mpds[m];
        if (!mpd.present) continue;

        for (int a = 0; a < ssp::MAX_APVS_PER_MPD; ++a) {
            auto &apv = mpd.apvs[a];
            if (!apv.present) continue;

            std::cout << "  MPD crate=" << mpd.crate_id
                      << " fiber=" << mpd.mpd_id
                      << " APV=" << a
                      << " strips=" << apv.nstrips;
            if (apv.has_online_cm) {
                std::cout << " online_cm=[";
                for (int t = 0; t < ssp::SSP_TIME_SAMPLES; ++t) {
                    if (t) std::cout << ",";
                    std::cout << apv.online_cm[t];
                }
                std::cout << "]";
            }
            std::cout << "\n";

            // print first 8 strips with data (or all if ≤ 8)
            int printed = 0;
            for (int s = 0; s < ssp::APV_STRIP_SIZE && printed < 8; ++s) {
                if (!apv.hasStrip(s)) continue;
                printed++;
                std::cout << "    ch[" << std::setw(3) << s << "] =";
                for (int t = 0; t < ssp::SSP_TIME_SAMPLES; ++t)
                    std::cout << " " << std::setw(6) << apv.strips[s][t];
                std::cout << "\n";
            }
            if (apv.nstrips > 8)
                std::cout << "    ... (" << apv.nstrips - 8 << " more strips)\n";
        }
    }
    std::cout << "\n";
}

// -------------------------------------------------------------------------
// Mode: hits — show strip hits after GemSystem processing
// -------------------------------------------------------------------------
static void dumpHits(gem::GemSystem &sys, int phys_ev)
{
    int total_hits = 0;
    for (int d = 0; d < sys.GetNDetectors(); ++d)
        for (int p = 0; p < 2; ++p)
            total_hits += (int)sys.GetPlaneHits(d, p).size();

    std::cout << "--- Event " << phys_ev << ": "
              << total_hits << " strip hit(s) ---\n";

    for (int d = 0; d < sys.GetNDetectors(); ++d) {
        auto &det = sys.GetDetectors()[d];
        for (int p = 0; p < 2; ++p) {
            auto &hits = sys.GetPlaneHits(d, p);
            if (hits.empty()) continue;

            std::cout << "  " << det.name << " " << (p == 0 ? "X" : "Y")
                      << ": " << hits.size() << " hit(s)\n";

            int show = std::min((int)hits.size(), 12);
            for (int i = 0; i < show; ++i) {
                auto &h = hits[i];
                std::cout << "    strip=" << std::setw(4) << h.strip
                          << " pos=" << std::fixed << std::setprecision(2)
                          << std::setw(8) << h.position << "mm"
                          << " charge=" << std::setw(8) << std::setprecision(1)
                          << h.charge
                          << " tbin=" << h.max_timebin;
                if (h.cross_talk) std::cout << " [xtalk]";
                // show time samples
                if (!h.ts_adc.empty()) {
                    std::cout << "  ts=[";
                    for (size_t t = 0; t < h.ts_adc.size(); ++t) {
                        if (t) std::cout << ",";
                        std::cout << std::setprecision(0) << h.ts_adc[t];
                    }
                    std::cout << "]";
                }
                std::cout << "\n";
            }
            if ((int)hits.size() > show)
                std::cout << "    ... (" << hits.size() - show << " more)\n";
        }
    }
    std::cout << "\n";
}

// -------------------------------------------------------------------------
// Mode: clusters — show clusters and 2D reconstructed hits
// -------------------------------------------------------------------------
static void dumpClusters(gem::GemSystem &sys, int phys_ev)
{
    auto &all_hits = sys.GetAllHits();
    std::cout << "--- Event " << phys_ev << ": "
              << all_hits.size() << " reconstructed 2D hit(s) ---\n";

    for (int d = 0; d < sys.GetNDetectors(); ++d) {
        auto &det = sys.GetDetectors()[d];

        // show 1D clusters per plane
        for (int p = 0; p < 2; ++p) {
            auto &clusters = sys.GetPlaneClusters(d, p);
            if (clusters.empty()) continue;

            std::cout << "  " << det.name << " " << (p == 0 ? "X" : "Y")
                      << ": " << clusters.size() << " cluster(s)\n";

            int show = std::min((int)clusters.size(), 8);
            for (int i = 0; i < show; ++i) {
                auto &cl = clusters[i];
                std::cout << "    pos=" << std::fixed << std::setprecision(2)
                          << std::setw(8) << cl.position << "mm"
                          << " peak=" << std::setprecision(1) << std::setw(8) << cl.peak_charge
                          << " total=" << std::setw(8) << cl.total_charge
                          << " size=" << cl.hits.size()
                          << " tbin=" << cl.max_timebin;
                if (cl.cross_talk) std::cout << " [xtalk]";
                std::cout << "\n";
            }
            if ((int)clusters.size() > show)
                std::cout << "    ... (" << clusters.size() - show << " more)\n";
        }

        // show 2D hits
        auto &hits2d = sys.GetHits(d);
        if (!hits2d.empty()) {
            std::cout << "  " << det.name << " 2D hits: " << hits2d.size() << "\n";
            int show = std::min((int)hits2d.size(), 8);
            for (int i = 0; i < show; ++i) {
                auto &h = hits2d[i];
                std::cout << "    (" << std::fixed << std::setprecision(2)
                          << std::setw(8) << h.x << ", "
                          << std::setw(8) << h.y << ") mm"
                          << "  Qx=" << std::setprecision(0) << h.x_charge
                          << " Qy=" << h.y_charge
                          << " Nx=" << h.x_size << " Ny=" << h.y_size
                          << "\n";
            }
            if ((int)hits2d.size() > show)
                std::cout << "    ... (" << hits2d.size() - show << " more)\n";
        }
    }
    std::cout << "\n";
}

// -------------------------------------------------------------------------
// Mode: evdump — dump single event data to JSON file
// -------------------------------------------------------------------------
static int dumpEventJson(const ssp::SspEventData &ssp,
                         gem::GemSystem &sys,
                         int phys_ev,
                         int32_t trigger_num,
                         uint32_t trigger_bits,
                         const std::string &output_file)
{
    using json = nlohmann::json;
    auto r1 = [](float v) -> double { return std::round(v * 10.) / 10.; };
    auto r2 = [](float v) -> double { return std::round(v * 100.) / 100.; };

    json root;
    root["event_number"]  = phys_ev;
    root["trigger_number"] = trigger_num;
    root["trigger_bits"]   = trigger_bits;

    // --- raw APV data (before pedestal/CM/zero-sup) ---
    json raw_arr = json::array();
    for (int m = 0; m < ssp.nmpds; ++m) {
        auto &mpd = ssp.mpds[m];
        if (!mpd.present) continue;
        for (int a = 0; a < ssp::MAX_APVS_PER_MPD; ++a) {
            auto &apv = mpd.apvs[a];
            if (!apv.present) continue;

            json aj;
            aj["crate"] = mpd.crate_id;
            aj["mpd"]   = mpd.mpd_id;
            aj["adc"]   = a;

            int idx = sys.FindApvIndex(mpd.crate_id, mpd.mpd_id, a);
            if (idx >= 0) {
                auto &cfg = sys.GetApvConfig(idx);
                aj["det"]   = cfg.det_id;
                aj["plane"] = (cfg.plane_type == 0) ? "X" : "Y";
                aj["pos"]   = cfg.plane_index;
            }

            json ch_obj = json::object();
            for (int s = 0; s < ssp::APV_STRIP_SIZE; ++s) {
                if (!apv.hasStrip(s)) continue;
                json samples = json::array();
                for (int t = 0; t < ssp::SSP_TIME_SAMPLES; ++t)
                    samples.push_back(apv.strips[s][t]);
                ch_obj[std::to_string(s)] = samples;
            }
            aj["channels"] = ch_obj;
            raw_arr.push_back(aj);
        }
    }
    root["raw_apvs"] = raw_arr;

    // --- zero-suppressed APV channels (APV address preserved) ---
    json zs_arr = json::array();
    float zs_thres = sys.GetZeroSupThreshold();
    float xt_thres = sys.GetCrossTalkThreshold();
    for (int m = 0; m < ssp.nmpds; ++m) {
        auto &mpd = ssp.mpds[m];
        if (!mpd.present) continue;
        for (int a = 0; a < ssp::MAX_APVS_PER_MPD; ++a) {
            auto &apv = mpd.apvs[a];
            if (!apv.present) continue;
            int idx = sys.FindApvIndex(mpd.crate_id, mpd.mpd_id, a);
            if (idx < 0) continue;

            json ch_obj = json::object();
            auto &cfg = sys.GetApvConfig(idx);
            for (int ch = 0; ch < ssp::APV_STRIP_SIZE; ++ch) {
                if (!sys.IsChannelHit(idx, ch)) continue;

                float max_charge = -1e9f;
                short max_tb = 0;
                json ts = json::array();
                for (int t = 0; t < ssp::SSP_TIME_SAMPLES; ++t) {
                    float val = sys.GetProcessedAdc(idx, ch, t);
                    ts.push_back(r1(val));
                    if (val > max_charge) { max_charge = val; max_tb = static_cast<short>(t); }
                }
                bool xtalk = (max_charge < cfg.pedestal[ch].noise * xt_thres)
                          && (max_charge > cfg.pedestal[ch].noise * zs_thres);
                ch_obj[std::to_string(ch)] = {
                    {"charge", r1(max_charge)}, {"max_timebin", max_tb},
                    {"cross_talk", xtalk}, {"ts_adc", ts}
                };
            }
            if (!ch_obj.empty()) {
                zs_arr.push_back({
                    {"crate", mpd.crate_id}, {"mpd", mpd.mpd_id}, {"adc", a},
                    {"channels", ch_obj}
                });
            }
        }
    }
    root["zs_apvs"] = zs_arr;

    // --- per-detector: clusters, 2D hits ---
    json det_arr = json::array();
    auto &dets = sys.GetDetectors();
    for (int d = 0; d < sys.GetNDetectors(); ++d) {
        json dj;
        dj["id"]       = d;
        dj["name"]     = dets[d].name;
        dj["x_pitch"]  = dets[d].planes[0].pitch;
        dj["y_pitch"]  = dets[d].planes[1].pitch;
        dj["x_strips"] = dets[d].planes[0].n_apvs * 128;
        dj["y_strips"] = dets[d].planes[1].n_apvs * 128;

        for (int p = 0; p < 2; ++p) {
            std::string pre = (p == 0) ? "x" : "y";

            // 1D clusters
            auto &clusters = sys.GetPlaneClusters(d, p);
            json cl_arr = json::array();
            for (auto &cl : clusters) {
                json cj;
                cj["position"]     = r2(cl.position);
                cj["peak_charge"]  = r1(cl.peak_charge);
                cj["total_charge"] = r1(cl.total_charge);
                cj["max_timebin"]  = cl.max_timebin;
                cj["cross_talk"]   = cl.cross_talk;
                cj["size"]         = static_cast<int>(cl.hits.size());
                json hs = json::array();
                for (auto &sh : cl.hits) hs.push_back(sh.strip);
                cj["hit_strips"] = hs;
                cl_arr.push_back(cj);
            }
            dj[pre + "_clusters"] = cl_arr;
        }

        // 2D reconstructed hits
        auto &h2d = sys.GetHits(d);
        json h2d_arr = json::array();
        for (auto &h : h2d) {
            json hj;
            hj["x"] = r2(h.x);  hj["y"] = r2(h.y);
            hj["x_charge"] = r1(h.x_charge);
            hj["y_charge"] = r1(h.y_charge);
            hj["x_peak"]   = r1(h.x_peak);
            hj["y_peak"]   = r1(h.y_peak);
            hj["x_size"]   = h.x_size;
            hj["y_size"]   = h.y_size;
            h2d_arr.push_back(hj);
        }
        dj["hits_2d"] = h2d_arr;

        det_arr.push_back(dj);
    }
    root["detectors"] = det_arr;

    // write (binary mode to avoid BOM/encoding issues on Windows)
    std::ofstream of(output_file, std::ios::binary);
    if (!of.is_open()) {
        std::cerr << "Error: cannot write " << output_file << "\n";
        return 1;
    }
    of << root.dump(2) << "\n";
    of.close();
    std::cerr << "Written: " << output_file << "\n";
    return 0;
}

// -------------------------------------------------------------------------
// Summary accumulator
// -------------------------------------------------------------------------
struct EventStats {
    int nmpds       = 0;
    int napvs       = 0;
    int nstrips     = 0;
    int nhits_x     = 0;
    int nhits_y     = 0;
    int nclusters_x = 0;
    int nclusters_y = 0;
    int nhits_2d    = 0;
};

static void accumulateStats(const ssp::SspEventData &ssp,
                            gem::GemSystem *sys,
                            EventStats &st)
{
    for (int m = 0; m < ssp.nmpds; ++m) {
        if (!ssp.mpds[m].present) continue;
        st.nmpds++;
        for (int a = 0; a < ssp::MAX_APVS_PER_MPD; ++a) {
            if (!ssp.mpds[m].apvs[a].present) continue;
            st.napvs++;
            st.nstrips += ssp.mpds[m].apvs[a].nstrips;
        }
    }

    if (sys) {
        for (int d = 0; d < sys->GetNDetectors(); ++d) {
            st.nhits_x     += (int)sys->GetPlaneHits(d, 0).size();
            st.nhits_y     += (int)sys->GetPlaneHits(d, 1).size();
            st.nclusters_x += (int)sys->GetPlaneClusters(d, 0).size();
            st.nclusters_y += (int)sys->GetPlaneClusters(d, 1).size();
            st.nhits_2d    += (int)sys->GetHits(d).size();
        }
    }
}

// -------------------------------------------------------------------------
// Mode: ped — compute per-strip pedestals from raw SSP data
// -------------------------------------------------------------------------
struct StripAccum {
    double sum  = 0.;
    double sum2 = 0.;
    int    count = 0;

    void add(double v) { sum += v; sum2 += v * v; ++count; }
    double mean() const { return count > 0 ? sum / count : 0.; }
    double rms() const {
        if (count < 2) return 0.;
        double m = mean();
        double var = sum2 / count - m * m;
        return var > 0 ? std::sqrt(var) : 0.;
    }
};

// key: pack (crate, mpd, apv, strip) into uint64
static uint64_t packPedKey(int crate, int mpd, int apv, int strip)
{
    return (static_cast<uint64_t>(crate & 0xFFFF) << 48) |
           (static_cast<uint64_t>(mpd   & 0xFFFF) << 32) |
           (static_cast<uint64_t>(apv   & 0xFFFF) << 16) |
           static_cast<uint64_t>(strip  & 0xFFFF);
}

// MPD pedestal algorithm:
//   Per event, per APV:
//     1) For each time sample: sort 128 strip ADCs, discard highest 28 and
//        lowest 28, average the middle 72 → commonMode[ts].
//        Subtract commonMode[ts] from each strip.
//     2) For each strip: average the 6 CM-corrected time samples → value.
//        Accumulate value into per-strip mean/RMS.
//   After all events: offset = mean, noise = RMS.

static constexpr int CM_DISCARD = 28;  // discard top/bottom 28 of 128

static void accumulatePedestals(const ssp::SspEventData &ssp,
                                std::map<uint64_t, StripAccum> &accum)
{
    float sorted[ssp::APV_STRIP_SIZE];
    float cm_corrected[ssp::APV_STRIP_SIZE][ssp::SSP_TIME_SAMPLES];

    for (int m = 0; m < ssp.nmpds; ++m) {
        auto &mpd = ssp.mpds[m];
        if (!mpd.present) continue;

        for (int a = 0; a < ssp::MAX_APVS_PER_MPD; ++a) {
            auto &apv = mpd.apvs[a];
            if (!apv.present || apv.nstrips == 0) continue;

            // Step 1.1: common mode subtraction per time sample
            for (int t = 0; t < ssp::SSP_TIME_SAMPLES; ++t) {
                // collect all strip values for this time sample
                for (int s = 0; s < ssp::APV_STRIP_SIZE; ++s)
                    sorted[s] = apv.hasStrip(s) ? (float)apv.strips[s][t] : 0.f;

                // sort to find middle 72
                std::sort(sorted, sorted + ssp::APV_STRIP_SIZE);

                // average middle strips (skip lowest CM_DISCARD and highest CM_DISCARD)
                double cm_sum = 0.;
                int cm_count = ssp::APV_STRIP_SIZE - 2 * CM_DISCARD;
                for (int s = CM_DISCARD; s < ssp::APV_STRIP_SIZE - CM_DISCARD; ++s)
                    cm_sum += sorted[s];
                double common_mode = cm_sum / cm_count;

                // subtract common mode
                for (int s = 0; s < ssp::APV_STRIP_SIZE; ++s)
                    cm_corrected[s][t] = (apv.hasStrip(s) ? (float)apv.strips[s][t] : 0.f) - common_mode;
            }

            // Step 1.2: per-strip average across time samples, accumulate
            for (int s = 0; s < ssp::APV_STRIP_SIZE; ++s) {
                if (!apv.hasStrip(s)) continue;

                double avg = 0.;
                for (int t = 0; t < ssp::SSP_TIME_SAMPLES; ++t)
                    avg += cm_corrected[s][t];
                avg /= ssp::SSP_TIME_SAMPLES;

                uint64_t key = packPedKey(mpd.crate_id, mpd.mpd_id, a, s);
                accum[key].add(avg);
            }
        }
    }
}

static int writePedestals(const std::map<uint64_t, StripAccum> &accum,
                          const std::string &output_file)
{
    std::ofstream of(output_file, std::ios::binary);
    if (!of.is_open()) {
        std::cerr << "Error: cannot write " << output_file << "\n";
        return 1;
    }

    // group strips by APV: (crate, mpd, adc) → strip data
    struct ApvKey { int crate, mpd, adc; };
    std::map<uint64_t, ApvKey> apv_keys;
    std::map<uint64_t, std::vector<std::pair<int, const StripAccum*>>> apv_strips;

    for (auto &[key, acc] : accum) {
        int crate = (key >> 48) & 0xFFFF;
        int mpd   = (key >> 32) & 0xFFFF;
        int apv   = (key >> 16) & 0xFFFF;
        int strip = key & 0xFFFF;
        uint64_t akey = (static_cast<uint64_t>(crate) << 32) |
                        (static_cast<uint64_t>(mpd) << 16) | apv;
        apv_keys[akey] = {crate, mpd, apv};
        apv_strips[akey].emplace_back(strip, &acc);
    }

    nlohmann::json arr = nlohmann::json::array();
    int napvs = 0;
    for (auto &[akey, strips] : apv_strips) {
        auto &ak = apv_keys[akey];
        // sort by strip number
        std::sort(strips.begin(), strips.end());

        nlohmann::json offsets = nlohmann::json::array();
        nlohmann::json noises  = nlohmann::json::array();
        for (auto &[s, acc] : strips) {
            offsets.push_back(std::round(acc->mean() * 1000.) / 1000.);
            noises.push_back(std::round(acc->rms() * 10000.) / 10000.);
        }
        arr.push_back({
            {"crate", ak.crate}, {"mpd", ak.mpd}, {"adc", ak.adc},
            {"offset", offsets}, {"noise", noises}
        });
        napvs++;
    }
    of << arr.dump(2) << "\n";
    of.close();

    std::cerr << "Written: " << output_file << " (" << napvs << " APVs, "
              << accum.size() << " strips)\n";
    return 0;
}

// -------------------------------------------------------------------------
// Main
// -------------------------------------------------------------------------
static void usage(const char *prog)
{
    std::cerr
        << "GEM data diagnostic tool\n\n"
        << "Usage:\n"
        << "  " << prog << " <evio_file> -D <daq_config.json> [options]\n\n"
        << "Modes (default: summary):\n"
        << "  -m raw        Dump raw SSP-decoded APV data\n"
        << "  -m hits       Strip hits after pedestal/CM/zero-sup\n"
        << "  -m clusters   Full reconstruction: clusters + 2D hits\n"
        << "  -m summary    Per-event statistics table\n"
        << "  -m ped        Compute per-strip pedestals → output file\n"
        << "  -m evdump     Dump event(s) with 2D hits to JSON (-n K for first K)\n\n"
        << "Options:\n"
        << "  -D <file>     DAQ configuration (auto-searches daq_config.json if omitted)\n"
        << "  -G <file>     GEM map file (default: gem_map.json)\n"
        << "  -P <file>     GEM pedestal file (input, for hits/clusters)\n"
        << "  -o <file>     Output file (ped mode, default: gem_ped.json)\n"
        << "  -n <N>        Max physics events (default: 10, 0=all for ped)\n"
        << "  -t <bit>      Trigger bit filter (-1=all, default)\n"
        << "  -e <N>        Dump only physics event N (1-based)\n"
        << "  -z <sigma>    Override zero-suppression threshold (default: from gem_map)\n"
        << "  -f <filter>   Event filter for evdump: field=val[,val]:min_dets\n"
        << "                Fields: pos, plane (X/Y), match (+Y/-Y), orient, det\n"
        << "                Example: -f pos=10,11:3 (ZS hits from pos 10/11 in >=3 dets)\n";
}

int main(int argc, char *argv[])
{
    if (argc < 2) { usage(argv[0]); return 1; }

    std::string daq_config_file;
    std::string gem_map_file;
    std::string gem_ped_file;
    std::string output_file = "gem_ped.json";
    std::string mode = "summary";
    std::string filter_expr;
    int max_events  = 10;
    int trigger_bit = -1;   // -1 = accept all
    int target_event = 0;   // 0 = disabled
    float zerosup_override = -1.f;  // <0 = use gem_map default

    int opt;
    while ((opt = getopt(argc, argv, "D:G:P:o:m:n:t:e:z:f:h")) != -1) {
        switch (opt) {
        case 'D': daq_config_file = optarg; break;
        case 'G': gem_map_file = optarg; break;
        case 'P': gem_ped_file = optarg; break;
        case 'o': output_file = optarg; break;
        case 'm': mode = optarg; break;
        case 'n': max_events = std::atoi(optarg); break;
        case 't': trigger_bit = std::atoi(optarg); break;
        case 'e': target_event = std::atoi(optarg); max_events = 0; break;
        case 'f': filter_expr = optarg; break;
        case 'z': zerosup_override = std::atof(optarg); break;
        default:  usage(argv[0]); return 1;
        }
    }
    if (optind >= argc) { usage(argv[0]); return 1; }
    std::string evio_file = argv[optind];

    // ped mode defaults to all events
    if (mode == "ped" && max_events == 10)
        max_events = 0;

    // evdump mode: dump first N events passing filter
    //   default: 1 event;  -n K: K events;  -n 0: all matching events
    //   -e N overrides: dump that specific event regardless of filter
    //   -f <filter>: custom APV-based filter (default: require 2D hits)
    int evdump_limit = 1;
    int evdump_count = 0;
    EvdumpFilter evdump_filter;
    bool has_evdump_filter = false;
    if (mode == "evdump") {
        if (target_event == 0) {
            evdump_limit = (max_events == 10) ? 1 : max_events;
            max_events = 0;   // let evdump_limit control the loop
        }
        if (output_file == "gem_ped.json")
            output_file = "gem_event.json";
        if (!filter_expr.empty()) {
            if (!parseEvdumpFilter(filter_expr, evdump_filter))
                return 1;
            has_evdump_filter = true;
            std::cerr << "Filter   : " << evdump_filter.field << "=";
            for (size_t i = 0; i < evdump_filter.values.size(); ++i)
                std::cerr << (i ? "," : "") << evdump_filter.values[i];
            std::cerr << " in >=" << evdump_filter.min_dets << " dets\n";
        }
    }

    // validate mode
    bool need_gem = (mode == "hits" || mode == "clusters" || mode == "summary" || mode == "evdump");
    bool need_cluster = (mode == "clusters" || mode == "summary" || mode == "evdump");

    // auto-search for daq_config.json if not specified
    if (daq_config_file.empty()) {
        for (auto p : {"daq_config.json", "database/daq_config.json", "../database/daq_config.json"}) {
            std::ifstream f(p);
            if (f.good()) { daq_config_file = p; break; }
        }
    }

    // load DAQ config
    DaqConfig daq_cfg;
    if (daq_config_file.empty() || !load_daq_config(daq_config_file, daq_cfg)) {
        std::cerr << "Error: failed to load DAQ config"
                  << (daq_config_file.empty() ? " (not found)" : ": " + daq_config_file)
                  << "\n";
        return 1;
    }
    std::cerr << "DAQ config: " << daq_config_file
              << " (adc_format=" << daq_cfg.adc_format << ")\n";

    // resolve GEM map file (default: gem_map.json next to DAQ config)
    if (gem_map_file.empty() && need_gem) {
        // try same directory as daq config
        auto pos = daq_config_file.rfind('/');
        if (pos == std::string::npos) pos = daq_config_file.rfind('\\');
        std::string dir = (pos != std::string::npos) ? daq_config_file.substr(0, pos + 1) : "";
        gem_map_file = dir + "gem_map.json";
    }

    // initialize GEM system
    std::unique_ptr<gem::GemSystem> gem_sys;
    std::unique_ptr<gem::GemCluster> gem_clusterer;

    if (need_gem) {
        gem_sys = std::make_unique<gem::GemSystem>();
        gem_sys->Init(gem_map_file);
        std::cerr << "GEM map  : " << gem_map_file
                  << " (" << gem_sys->GetNDetectors() << " detectors)\n";

        if (!gem_ped_file.empty()) {
            gem_sys->LoadPedestals(gem_ped_file);
            std::cerr << "GEM peds : " << gem_ped_file << "\n";
        }

        if (zerosup_override >= 0.f) {
            gem_sys->SetZeroSupThreshold(zerosup_override);
            std::cerr << "Zero-sup : " << zerosup_override << " sigma (override)\n";
        }

        if (need_cluster)
            gem_clusterer = std::make_unique<gem::GemCluster>();
    }

    // open EVIO file
    EvChannel ch;
    ch.SetConfig(daq_cfg);
    if (ch.Open(evio_file) != status::success) {
        std::cerr << "Error: cannot open " << evio_file << "\n";
        return 1;
    }
    std::cerr << "File     : " << evio_file << "\n\n";

    // trigger filter
    uint32_t trigger_mask = 0;
    if (trigger_bit >= 0) {
        trigger_mask = 1u << trigger_bit;
        std::cerr << "Trigger  : bit " << trigger_bit
                  << " (mask 0x" << std::hex << trigger_mask << std::dec << ")\n";
    }

    // summary mode header
    if (mode == "summary") {
        std::cout << std::setw(6) << "ev#"
                  << std::setw(10) << "trigger#"
                  << std::setw(6) << "MPDs"
                  << std::setw(6) << "APVs"
                  << std::setw(8) << "strips"
                  << std::setw(8) << "hits_X"
                  << std::setw(8) << "hits_Y"
                  << std::setw(8) << "clus_X"
                  << std::setw(8) << "clus_Y"
                  << std::setw(8) << "2D_hits"
                  << "\n";
        std::cout << std::string(76, '-') << "\n";
    }

    // event loop
    auto event_ptr = std::make_unique<fdec::EventData>();
    auto ssp_ptr   = std::make_unique<ssp::SspEventData>();
    auto &event    = *event_ptr;
    auto &ssp_evt  = *ssp_ptr;

    int phys_count = 0;
    int ssp_events = 0;

    // totals for summary
    EventStats totals;

    // pedestal accumulator
    std::map<uint64_t, StripAccum> ped_accum;

    while (ch.Read() == status::success) {
        if (!ch.Scan()) continue;
        if (ch.GetEventType() != EventType::Physics) continue;

        for (int i = 0; i < ch.GetNEvents(); ++i) {
            ssp_evt.clear();
            if (!ch.DecodeEvent(i, event, &ssp_evt)) continue;
            phys_count++;

            // trigger filter
            if (trigger_mask && !(event.info.trigger_bits & trigger_mask))
                continue;

            // target event filter
            if (target_event > 0 && phys_count != target_event)
                continue;

            // skip events with no SSP data
            if (ssp_evt.nmpds == 0) {
                if (target_event > 0) {
                    std::cout << "Event " << phys_count << ": no SSP data\n";
                    goto done;
                }
                continue;
            }

            ssp_events++;

            // GEM processing
            if (gem_sys) {
                gem_sys->Clear();
                gem_sys->ProcessEvent(ssp_evt);
                if (gem_clusterer)
                    gem_sys->Reconstruct(*gem_clusterer);
            }

            // output based on mode
            if (mode == "ped") {
                accumulatePedestals(ssp_evt, ped_accum);
                if (ssp_events % 1000 == 0)
                    std::cerr << "  " << ssp_events << " events...\r" << std::flush;
            }
            else if (mode == "raw") {
                std::cout << "[physics event " << phys_count
                          << " trigger#=" << event.info.trigger_number
                          << " bits=0x" << std::hex << event.info.trigger_bits
                          << std::dec << "]\n";
                dumpRawSsp(ssp_evt, phys_count);
            }
            else if (mode == "hits") {
                std::cout << "[physics event " << phys_count
                          << " trigger#=" << event.info.trigger_number
                          << " bits=0x" << std::hex << event.info.trigger_bits
                          << std::dec << "]\n";
                dumpHits(*gem_sys, phys_count);
            }
            else if (mode == "clusters") {
                std::cout << "[physics event " << phys_count
                          << " trigger#=" << event.info.trigger_number
                          << " bits=0x" << std::hex << event.info.trigger_bits
                          << std::dec << "]\n";
                dumpClusters(*gem_sys, phys_count);
            }
            else if (mode == "evdump") {
                // apply filter (unless -e targets a specific event)
                if (target_event == 0) {
                    if (has_evdump_filter) {
                        if (!checkEvdumpFilter(evdump_filter, *gem_sys, ssp_evt))
                            continue;
                    } else {
                        if (gem_sys->GetAllHits().empty())
                            continue;
                    }
                }

                // multi-event: append event number to filename
                std::string out = output_file;
                if (target_event == 0 && evdump_limit != 1) {
                    auto dot = output_file.rfind('.');
                    std::string suffix = "_" + std::to_string(phys_count);
                    if (dot != std::string::npos)
                        out = output_file.substr(0, dot) + suffix + output_file.substr(dot);
                    else
                        out = output_file + suffix;
                }

                if (dumpEventJson(ssp_evt, *gem_sys, phys_count,
                                  event.info.trigger_number,
                                  event.info.trigger_bits, out) != 0)
                    return 1;

                if (evdump_limit > 0 && ++evdump_count >= evdump_limit)
                    goto done;
            }
            else { // summary
                EventStats st;
                accumulateStats(ssp_evt, gem_sys.get(), st);

                std::cout << std::setw(6) << phys_count
                          << std::setw(10) << event.info.trigger_number
                          << std::setw(6) << st.nmpds
                          << std::setw(6) << st.napvs
                          << std::setw(8) << st.nstrips
                          << std::setw(8) << st.nhits_x
                          << std::setw(8) << st.nhits_y
                          << std::setw(8) << st.nclusters_x
                          << std::setw(8) << st.nclusters_y
                          << std::setw(8) << st.nhits_2d
                          << "\n";

                // accumulate totals
                totals.nmpds       += st.nmpds;
                totals.napvs       += st.napvs;
                totals.nstrips     += st.nstrips;
                totals.nhits_x     += st.nhits_x;
                totals.nhits_y     += st.nhits_y;
                totals.nclusters_x += st.nclusters_x;
                totals.nclusters_y += st.nclusters_y;
                totals.nhits_2d    += st.nhits_2d;
            }

            if (target_event > 0)
                goto done;
            if (max_events > 0 && ssp_events >= max_events)
                goto done;
        }
        if (max_events > 0 && ssp_events >= max_events)
            break;
    }

done:
    ch.Close();

    // pedestal output
    if (mode == "ped" && !ped_accum.empty()) {
        std::cerr << "Computed pedestals from " << ssp_events << " events, "
                  << ped_accum.size() << " strips\n";
        return writePedestals(ped_accum, output_file);
    }

    // summary footer
    if (mode == "summary" && ssp_events > 0) {
        std::cout << std::string(76, '-') << "\n";
        std::cout << "Totals: " << ssp_events << " events with SSP data"
                  << " (of " << phys_count << " physics events)\n";
        if (ssp_events > 0) {
            std::cout << "  Avg per event:"
                      << " MPDs=" << std::fixed << std::setprecision(1)
                      << (float)totals.nmpds / ssp_events
                      << " APVs=" << (float)totals.napvs / ssp_events
                      << " strips=" << (float)totals.nstrips / ssp_events
                      << "\n";
            if (gem_sys) {
                std::cout << "  Avg per event:"
                          << " hits_X=" << (float)totals.nhits_x / ssp_events
                          << " hits_Y=" << (float)totals.nhits_y / ssp_events
                          << " clus_X=" << (float)totals.nclusters_x / ssp_events
                          << " clus_Y=" << (float)totals.nclusters_y / ssp_events
                          << " 2D=" << (float)totals.nhits_2d / ssp_events
                          << "\n";
            }
        }
    }
    else if (ssp_events == 0) {
        std::cerr << "No events with SSP data found in " << phys_count
                  << " physics events.\n";
    }

    std::cerr << "Done: parsed " << phys_count << " physics events, "
              << ssp_events << " with SSP data.\n";
    return 0;
}
