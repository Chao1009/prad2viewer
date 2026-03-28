#include "GemSystem.h"
#include "GemCluster.h"
#include "SspData.h"
#include <nlohmann/json.hpp>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <sstream>

using namespace gem;

//=============================================================================
// Construction / destruction
//=============================================================================

GemSystem::GemSystem() = default;
GemSystem::~GemSystem() = default;

//=============================================================================
// Init — load GEM map from JSON
//=============================================================================

void GemSystem::Init(const std::string &map_file)
{
    std::ifstream f(map_file);
    if (!f.is_open()) {
        std::cerr << "GemSystem::Init: cannot open " << map_file << std::endl;
        return;
    }

    nlohmann::json j;
    try { j = nlohmann::json::parse(f, nullptr, true, true); }
    catch (const nlohmann::json::parse_error &e) {
        std::cerr << "GemSystem::Init: parse error: " << e.what() << std::endl;
        return;
    }

    // --- parse layers/detectors ---
    detectors_.clear();
    if (j.contains("layers")) {
        for (auto &layer : j["layers"]) {
            DetectorConfig det;
            det.id   = layer.value("id", static_cast<int>(detectors_.size()));
            det.name = layer.value("name", "GEM" + std::to_string(det.id));
            det.type = layer.value("type", "PRADGEM");

            det.planes[0].type   = 0; // X
            det.planes[0].n_apvs = layer.value("x_apvs", 12);
            det.planes[0].pitch  = layer.value("x_pitch", 0.4f);
            det.planes[0].size   = det.planes[0].n_apvs * APV_STRIP_SIZE * det.planes[0].pitch;

            det.planes[1].type   = 1; // Y
            det.planes[1].n_apvs = layer.value("y_apvs", 24);
            det.planes[1].pitch  = layer.value("y_pitch", 0.4f);
            det.planes[1].size   = det.planes[1].n_apvs * APV_STRIP_SIZE * det.planes[1].pitch;

            detectors_.push_back(det);
        }
    }

    // --- parse global parameters ---
    apv_channels_    = j.value("apv_channels", 128);
    readout_center_  = j.value("readout_center", 32);
    common_thres_    = j.value("common_mode_threshold", 20.f);
    zerosup_thres_   = j.value("zero_suppression_threshold", 5.f);
    crosstalk_thres_ = j.value("cross_talk_threshold", 8.f);
    online_zero_sup_ = j.value("online_zero_suppression", false);
    position_res_    = j.value("position_resolution", 0.08f);

    // strip-level cuts
    reject_first_timebin_ = j.value("reject_first_timebin", true);
    reject_last_timebin_  = j.value("reject_last_timebin", true);
    min_peak_adc_         = j.value("min_peak_adc", 0.f);
    min_sum_adc_          = j.value("min_sum_adc", 0.f);

    // XY matching cuts (stored here, applied via GemCluster)
    match_adc_asymmetry_  = j.value("match_adc_asymmetry", 0.8f);
    match_time_diff_      = j.value("match_time_diff", 50.f);
    match_ts_period_      = j.value("match_ts_period", 25.f);

    // --- parse APV entries ---
    apvs_.clear();
    apv_map_.clear();
    if (j.contains("apvs")) {
        for (auto &entry : j["apvs"]) {
            ApvConfig apv;
            apv.crate_id    = entry.value("crate", -1);
            apv.mpd_id      = entry.value("mpd", -1);
            apv.adc_ch      = entry.value("adc", -1);
            apv.det_id      = entry.value("det", 0);

            std::string plane_str = entry.value("plane", "X");
            apv.plane_type = (plane_str == "Y" || plane_str == "1") ? 1 : 0;

            apv.orient      = entry.value("orient", 0);
            apv.plane_index = entry.value("pos", 0);
            apv.det_pos     = entry.value("det_pos", 0);

            apv.pin_rotate   = entry.value("pin_rotate", 0);
            apv.shared_pos   = entry.value("shared_pos", -1);
            apv.hybrid_board = entry.value("hybrid_board", true);
            apv.match        = entry.value("match", "");

            int idx = static_cast<int>(apvs_.size());
            apv_map_[packApvKey(apv.crate_id, apv.mpd_id, apv.adc_ch)] = idx;
            apvs_.push_back(apv);
        }
    }

    // --- allocate working data ---
    apv_work_.resize(apvs_.size());
    for (size_t i = 0; i < apvs_.size(); ++i)
        buildStripMap(static_cast<int>(i));

    // --- allocate per-plane and per-detector storage ---
    plane_data_.resize(detectors_.size());
    det_hits_.resize(detectors_.size());
}

//=============================================================================
// LoadPedestals — load per-strip pedestal from JSON file
//
// Format: [{"crate":49,"mpd":0,"adc":0,"offset":[...128...],"noise":[...128...]}, ...]
//=============================================================================

void GemSystem::LoadPedestals(const std::string &ped_file)
{
    std::ifstream f(ped_file);
    if (!f.is_open()) {
        std::cerr << "GemSystem::LoadPedestals: cannot open " << ped_file << std::endl;
        return;
    }

    nlohmann::json j;
    try { j = nlohmann::json::parse(f, nullptr, true, true); }
    catch (const nlohmann::json::parse_error &e) {
        std::cerr << "GemSystem::LoadPedestals: parse error: " << e.what() << std::endl;
        return;
    }

    int loaded = 0;
    for (auto &entry : j) {
        int crate = entry.value("crate", -1);
        int mpd   = entry.value("mpd", -1);
        int adc   = entry.value("adc", -1);

        int idx = FindApvIndex(crate, mpd, adc);
        if (idx < 0) continue;

        auto &offsets = entry["offset"];
        auto &noises  = entry["noise"];
        int nstrips = std::min({(int)offsets.size(), (int)noises.size(), APV_STRIP_SIZE});

        for (int s = 0; s < nstrips; ++s) {
            apvs_[idx].pedestal[s].offset = offsets[s].get<float>();
            apvs_[idx].pedestal[s].noise  = noises[s].get<float>();
        }
        loaded += nstrips;
    }
    std::cerr << "GemSystem: loaded " << loaded << " pedestal entries from " << ped_file << "\n";
}

//=============================================================================
// LoadCommonModeRange — load per-APV common mode range
//
// Format: crate  mpd  adc  cm_min  cm_max
//=============================================================================

void GemSystem::LoadCommonModeRange(const std::string &cm_file)
{
    std::ifstream f(cm_file);
    if (!f.is_open()) {
        std::cerr << "GemSystem::LoadCommonModeRange: cannot open " << cm_file << std::endl;
        return;
    }

    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        int crate, mpd, adc;
        float cm_min, cm_max;
        if (!(iss >> crate >> mpd >> adc >> cm_min >> cm_max)) continue;

        int idx = FindApvIndex(crate, mpd, adc);
        if (idx < 0) continue;

        apvs_[idx].cm_range_min = cm_min;
        apvs_[idx].cm_range_max = cm_max;
    }
}

//=============================================================================
// Clear — reset per-event data
//=============================================================================

void GemSystem::Clear()
{
    for (auto &w : apv_work_)
        std::memset(w.hit_pos, 0, sizeof(w.hit_pos));
    for (auto &pd : plane_data_) {
        pd[0].hits.clear();
        pd[0].clusters.clear();
        pd[1].hits.clear();
        pd[1].clusters.clear();
    }
    for (auto &dh : det_hits_)
        dh.clear();
    all_hits_.clear();
}

//=============================================================================
// ProcessEvent — decode SSP data → strip hits
//=============================================================================

void GemSystem::ProcessEvent(const ssp::SspEventData &evt)
{
    for (int mi = 0; mi < evt.nmpds; ++mi) {
        auto &mpd = evt.mpds[mi];
        if (!mpd.present) continue;

        for (int ai = 0; ai < ssp::MAX_APVS_PER_MPD; ++ai) {
            auto &apv = mpd.apvs[ai];
            if (!apv.present) continue;

            int idx = FindApvIndex(apv.addr.crate_id, apv.addr.mpd_id, apv.addr.adc_ch);
            if (idx < 0) continue;

            processApv(idx, apv);
        }
    }
}

//=============================================================================
// Reconstruct — run clustering on all planes, then 2D matching
//=============================================================================

void GemSystem::Reconstruct(GemCluster &clusterer)
{
    // apply XY matching cuts from gem_map config
    auto cfg = clusterer.GetConfig();
    cfg.match_adc_asymmetry = match_adc_asymmetry_;
    cfg.match_time_diff     = match_time_diff_;
    cfg.ts_period           = match_ts_period_;
    clusterer.SetConfig(cfg);

    for (int d = 0; d < static_cast<int>(detectors_.size()); ++d) {
        // cluster X and Y planes
        for (int p = 0; p < 2; ++p) {
            auto &pd = plane_data_[d][p];
            clusterer.FormClusters(pd.hits, pd.clusters);
        }

        // Cartesian reconstruction: match X and Y clusters
        auto &xc = plane_data_[d][0].clusters;
        auto &yc = plane_data_[d][1].clusters;
        clusterer.CartesianReconstruct(xc, yc, det_hits_[d], d, position_res_);

        // accumulate all hits
        all_hits_.insert(all_hits_.end(), det_hits_[d].begin(), det_hits_[d].end());
    }
}

//=============================================================================
// FindApvIndex — O(1) lookup
//=============================================================================

int GemSystem::FindApvIndex(int crate, int mpd, int adc) const
{
    auto it = apv_map_.find(packApvKey(crate, mpd, adc));
    return (it != apv_map_.end()) ? it->second : -1;
}

//=============================================================================
// Accessors
//=============================================================================

const std::vector<StripHit>& GemSystem::GetPlaneHits(int det, int plane) const
{
    return plane_data_[det][plane].hits;
}

const std::vector<StripCluster>& GemSystem::GetPlaneClusters(int det, int plane) const
{
    return plane_data_[det][plane].clusters;
}

const std::vector<GEMHit>& GemSystem::GetHits(int det) const
{
    return det_hits_[det];
}

//=============================================================================
// processApv — per-APV: pedestal subtraction, common mode, zero suppression
//=============================================================================

// Macro to compute data index: raw[ch + ts * (APV_STRIP_SIZE + 1)]
// The +1 accounts for the APV header word per time sample (MPD_APV_TS_LEN = 129)
// For our flat buffer, we use a simpler layout: raw[ts * APV_STRIP_SIZE + ch]
#define RAW_IDX(ch, ts) ((ts) * APV_STRIP_SIZE + (ch))

void GemSystem::processApv(int apv_idx, const ssp::ApvData &data)
{
    auto &cfg = apvs_[apv_idx];
    auto &work = apv_work_[apv_idx];

    // --- copy raw data into working buffer ---
    for (int ch = 0; ch < APV_STRIP_SIZE; ++ch) {
        for (int ts = 0; ts < SSP_TIME_SAMPLES; ++ts) {
            work.raw[RAW_IDX(ch, ts)] = static_cast<float>(data.strips[ch][ts]);
        }
    }

    // --- common mode correction for each time sample ---
    for (int ts = 0; ts < SSP_TIME_SAMPLES; ++ts) {
        float *buf = &work.raw[ts * APV_STRIP_SIZE];

        // subtract pedestal offset
        if (!online_zero_sup_) {
            for (int ch = 0; ch < APV_STRIP_SIZE; ++ch)
                buf[ch] -= cfg.pedestal[ch].offset;
        }

        // compute and subtract common mode (sorting algorithm)
        float cm = commonModeSorting(buf, APV_STRIP_SIZE, apv_idx);
        for (int ch = 0; ch < APV_STRIP_SIZE; ++ch)
            buf[ch] -= cm;
    }

    // --- zero suppression ---
    for (int ch = 0; ch < APV_STRIP_SIZE; ++ch) {
        float avg = 0.f;
        for (int ts = 0; ts < SSP_TIME_SAMPLES; ++ts)
            avg += work.raw[RAW_IDX(ch, ts)];
        avg /= SSP_TIME_SAMPLES;

        work.hit_pos[ch] = (avg > cfg.pedestal[ch].noise * zerosup_thres_);
    }

    // --- collect hits to plane ---
    collectHits(apv_idx);
}

//=============================================================================
// commonModeSorting — MPD version: remove top N high-ADC strips from average
//=============================================================================

float GemSystem::commonModeSorting(float *buf, int size, [[maybe_unused]] int apv_idx)
{
    float sum = 0.f;
    int count = 0;

    // Track the top NUM_HIGH_STRIPS highest values to exclude
    std::vector<float> high_vals(NUM_HIGH_STRIPS, -9999.f);

    for (int i = 0; i < size; ++i) {
        sum += buf[i];
        count++;

        // Maintain sorted list of top N highest values
        if (buf[i] > high_vals[0]) {
            // Find insertion point and shift
            int pos = 0;
            while (pos < NUM_HIGH_STRIPS - 1 && high_vals[pos + 1] < buf[i])
                pos++;
            for (int j = 0; j < pos; ++j)
                high_vals[j] = high_vals[j + 1];
            high_vals[pos] = buf[i];
        }
    }

    // Subtract the top N values from sum
    for (int i = 0; i < NUM_HIGH_STRIPS && count > 1; ++i) {
        sum -= high_vals[i];
        count--;
    }

    return (count > 0) ? sum / static_cast<float>(count) : 0.f;
}

//=============================================================================
// commonModeDanning — Danning algorithm with common mode range
//=============================================================================

float GemSystem::commonModeDanning(float *buf, int size, int apv_idx)
{
    auto &cfg = apvs_[apv_idx];

    // Step 1: average A — only values within common mode range
    float avgA = 0.f;
    int countA = 0;
    for (int i = 0; i < size; ++i) {
        if (buf[i] >= cfg.cm_range_min && buf[i] <= cfg.cm_range_max) {
            avgA += buf[i];
            countA++;
        }
    }
    if (countA == 0) return 0.f;
    avgA /= static_cast<float>(countA);

    // Step 2: average B — values below avgA + RMS_THRESHOLD * noise
    static constexpr float RMS_THRESHOLD = 3.f;
    float avgB = 0.f;
    int countB = 0;
    for (int i = 0; i < size; ++i) {
        if (buf[i] < avgA + RMS_THRESHOLD * cfg.pedestal[i].noise) {
            avgB += buf[i];
            countB++;
        }
    }

    return (countB > 0) ? avgB / static_cast<float>(countB) : 0.f;
}

//=============================================================================
// collectHits — gather zero-suppressed hits into plane data
//=============================================================================

void GemSystem::collectHits(int apv_idx)
{
    auto &cfg = apvs_[apv_idx];
    auto &work = apv_work_[apv_idx];

    if (cfg.det_id < 0 || cfg.det_id >= static_cast<int>(detectors_.size()))
        return;
    if (cfg.plane_type < 0 || cfg.plane_type > 1)
        return;

    auto &det = detectors_[cfg.det_id];
    auto &plane = det.planes[cfg.plane_type];
    auto &hits = plane_data_[cfg.det_id][cfg.plane_type].hits;

    for (int ch = 0; ch < APV_STRIP_SIZE; ++ch) {
        if (!work.hit_pos[ch]) continue;

        int plane_strip = work.strip_map[ch];
        if (plane_strip < 0) continue;

        // Find max charge and max timebin
        float max_charge = -1e9f;
        float sum_adc = 0.f;
        short max_tb = 0;
        std::vector<float> ts_adc(SSP_TIME_SAMPLES);
        for (int ts = 0; ts < SSP_TIME_SAMPLES; ++ts) {
            float val = work.raw[RAW_IDX(ch, ts)];
            ts_adc[ts] = val;
            sum_adc += val;
            if (val > max_charge) {
                max_charge = val;
                max_tb = static_cast<short>(ts);
            }
        }

        // Strip-level cuts
        if (reject_first_timebin_ && max_tb == 0) continue;
        if (reject_last_timebin_  && max_tb == SSP_TIME_SAMPLES - 1) continue;
        if (min_peak_adc_ > 0.f   && max_charge < min_peak_adc_) continue;
        if (min_sum_adc_  > 0.f   && sum_adc < min_sum_adc_) continue;

        // Calculate physical position
        float pos = static_cast<float>(plane_strip) * plane.pitch
                    - plane.size * 0.5f + plane.pitch * 0.5f;

        // Check cross-talk
        bool xtalk = (max_charge < cfg.pedestal[ch].noise * crosstalk_thres_)
                     && (max_charge > cfg.pedestal[ch].noise * zerosup_thres_);

        StripHit hit;
        hit.strip       = plane_strip;
        hit.charge      = max_charge;
        hit.max_timebin = max_tb;
        hit.position    = pos;
        hit.cross_talk  = xtalk;
        hit.ts_adc      = std::move(ts_adc);
        hits.push_back(std::move(hit));
    }
}

#undef RAW_IDX

//=============================================================================
// buildStripMap — compute APV channel → plane strip mapping
//
// Implements the full MapStrip pipeline (from PRadAnalyzer/mpd_gem_view_ssp):
//   1. APV25 internal channel mapping (chip wiring, universal)
//   2. Hybrid board pin conversion (MPD electronics only)
//   3. Readout strip scaling (configurable offset: 32 normal, 48 for special APVs)
//   4. 7-bit mask
//   5. Orient flip
//   6. Plane-wide strip number with configurable offset
//=============================================================================

void GemSystem::buildStripMap(int apv_idx)
{
    auto &cfg = apvs_[apv_idx];
    auto &work = apv_work_[apv_idx];

    const int N = apv_channels_;
    const int center = readout_center_;

    // derive from physical parameters
    int readout_off = center + cfg.pin_rotate;
    int eff_pos = (cfg.shared_pos >= 0) ? cfg.shared_pos : cfg.plane_index;
    int plane_shift = (eff_pos - cfg.plane_index) * N - cfg.pin_rotate;

    for (int ch = 0; ch < N; ++ch) {
        // Step 1: APV25 internal channel mapping (chip wiring, universal)
        int strip = 32 * (ch % 4) + 8 * (ch / 4) - 31 * (ch / 16);

        // Step 2: hybrid board pin conversion (MPD electronics)
        if (cfg.hybrid_board)
            strip = strip + 1 + strip % 4 - 5 * ((strip / 4) % 2);

        // Step 3: readout strip mapping (odd/even fan-out around center)
        // readout_off > 0: apply mapping; 0: skip (steps 1+2 give final strip)
        if (readout_off > 0) {
            if (strip & 1)
                strip = readout_off - (strip + 1) / 2;
            else
                strip = readout_off + strip / 2;
        }

        // Step 4: channel mask
        strip &= (N - 1);

        // Step 5: orient flip
        if (cfg.orient == 1)
            strip = (N - 1) - strip;

        // Step 6: plane-wide strip number
        strip += plane_shift + cfg.plane_index * N;

        work.strip_map[ch] = strip;
    }
}
