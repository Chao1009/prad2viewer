#include "PulseTemplateStore.h"

#include <nlohmann/json.hpp>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>
#include <cstdio>
#include <cstdlib>

using namespace fdec;
using nlohmann::json;

namespace {

// Pull "<roc>_<slot>_<channel>" out of the entry's `channel_id` string.
// Returns true on success, false if the format doesn't match.  We tolerate
// unmapped channels (where the channel_id is itself the key, e.g. for
// LMS lines that the daq_map doesn't resolve) — same parser handles both.
bool parse_channel_id(const std::string &s, int &roc, int &slot, int &ch)
{
    int a = 0, b = 0, c = 0;
    if (std::sscanf(s.c_str(), "%d_%d_%d", &a, &b, &c) != 3) return false;
    roc = a; slot = b; ch = c;
    return true;
}

// Read a {"median": <num>, "mad": <num>} sub-object's median field.  The
// fitter writes NaN when no pulses contributed; we propagate NaN.
double read_median(const json &j, const char *key)
{
    if (!j.contains(key)) return std::nan("");
    const auto &sub = j[key];
    if (!sub.is_object() || !sub.contains("median")) return std::nan("");
    const auto &m = sub["median"];
    if (m.is_number()) return m.get<double>();
    return std::nan("");
}

} // anon

void PulseTemplateStore::Clear()
{
    by_csc_.clear();
    global_      = {};
    valid_       = false;
    has_global_  = false;
    n_channels_loaded_ = 0;
    n_good_      = 0;
}

namespace {

// Fill `tmpl.grid` with the unit-amplitude two-tau template on the
// hardwired 8× oversampled time axis t_i = i · (clk_ns / GRID_OVERSAMPLE).
// Sets tmpl.grid_clk_ns so consumers can detect "grid filled" via
// grid_clk_ns > 0.  See PulseTemplate doc for why the oversample factor
// isn't a runtime knob.
void fill_grid(PulseTemplate &tmpl, float clk_ns)
{
    const float dt = clk_ns / static_cast<float>(PulseTemplate::GRID_OVERSAMPLE);
    tmpl.grid_clk_ns = dt;
    const float tr = tmpl.tau_r_ns;
    const float tf = tmpl.tau_f_ns;
    for (int i = 0; i < PulseTemplate::GRID_N; ++i) {
        const float t = i * dt;
        tmpl.grid[i] = (1.0f - std::exp(-t / tr)) * std::exp(-t / tf);
    }
    // i=0 evaluates to 0 (1 - 1)·exp(0) = 0, exact.
}

} // anon

bool PulseTemplateStore::LoadFromFile(const std::string &path,
                                      const WaveConfig &cfg)
{
    Clear();

    std::ifstream f(path);
    if (!f.is_open()) {
        std::cerr << "[PulseTemplateStore] WARN: cannot open " << path
                  << " — deconv will fall back to non-deconv mode.\n";
        return false;
    }
    json j;
    try { j = json::parse(f, nullptr, true, true); }   // allow comments
    catch (const json::parse_error &e) {
        std::cerr << "[PulseTemplateStore] WARN: parse error in "
                  << path << " (" << e.what()
                  << ") — deconv will fall back to non-deconv mode.\n";
        return false;
    }

    const auto &dcfg = cfg.nnls_deconv;
    const float clk_ns = (cfg.clk_mhz > 0.0f) ? (1000.0f / cfg.clk_mhz) : 4.0f;

    // Effective gates.  min_pulses / chi2_max come from the fitter's
    // _meta block (those are properties of the fit, not the analyzer);
    // the τ ranges come from the analyzer config.
    int   min_pulses = 50;
    float chi2_max   = 3.0f;
    if (j.contains("_meta")) {
        const auto &meta = j["_meta"];
        if (meta.contains("min_pulses") && meta["min_pulses"].is_number())
            min_pulses = meta["min_pulses"].get<int>();
        if (meta.contains("chi2_max") && meta["chi2_max"].is_number())
            chi2_max = meta["chi2_max"].get<float>();
    }
    const float tr_lo = dcfg.tau_r_min_ns, tr_hi = dcfg.tau_r_max_ns;
    const float tf_lo = dcfg.tau_f_min_ns, tf_hi = dcfg.tau_f_max_ns;

    std::vector<float> good_taur, good_tauf;

    for (auto it = j.begin(); it != j.end(); ++it) {
        const std::string &name = it.key();
        if (!name.empty() && name[0] == '_') continue;     // skip _meta etc.
        const auto &rec = it.value();
        if (!rec.is_object())                  continue;

        const double tr   = read_median(rec, "tau_r_ns");
        const double tf   = read_median(rec, "tau_f_ns");
        const double chi2 = read_median(rec, "chi2_per_dof");
        int n_used = 0;
        if (rec.contains("n_pulses_used") && rec["n_pulses_used"].is_number())
            n_used = rec["n_pulses_used"].get<int>();

        if (!std::isfinite(tr) || !std::isfinite(tf)) continue;

        // Resolve (roc, slot, channel) from channel_id field; if absent
        // OR malformed (e.g. unmapped channel keyed only by name), skip
        // this entry — the analyzer can't look it up by name anyway.
        if (!rec.contains("channel_id") || !rec["channel_id"].is_string())
            continue;
        int roc = -1, slot = -1, ch = -1;
        if (!parse_channel_id(rec["channel_id"].get<std::string>(),
                              roc, slot, ch))
            continue;

        Entry e{};
        e.tmpl.tau_r_ns    = static_cast<float>(tr);
        e.tmpl.tau_f_ns    = static_cast<float>(tf);
        e.tmpl.is_global   = false;
        e.tmpl.grid_clk_ns = 0.0f;          // filled below if gates pass

        const bool gates_ok =
            (n_used >= min_pulses) &&
            std::isfinite(chi2) && (chi2 < chi2_max) &&
            (tr >= tr_lo && tr <= tr_hi) &&
            (tf >= tf_lo && tf <= tf_hi);
        e.is_good = gates_ok;
        if (gates_ok) {
            fill_grid(e.tmpl, clk_ns);              // precompute to skip exp() in hot loop
            ++n_good_;
            good_taur.push_back(e.tmpl.tau_r_ns);
            good_tauf.push_back(e.tmpl.tau_f_ns);
        }

        by_csc_.emplace(pack_key(roc, slot, ch), e);
        ++n_channels_loaded_;
    }

    if (n_channels_loaded_ == 0) {
        std::cerr << "[PulseTemplateStore] WARN: " << path
                  << " contained no usable channel templates"
                  << " — deconv will fall back to non-deconv mode.\n";
        return false;
    }

    // Synthesise the global-median template if at least a handful of
    // channels survived gates — bar set low so a partial run is still
    // useful as a fallback.
    if (good_taur.size() >= 5) {
        std::sort(good_taur.begin(), good_taur.end());
        std::sort(good_tauf.begin(), good_tauf.end());
        const auto med = [](const std::vector<float> &v) {
            const size_t n = v.size();
            return (n % 2 == 1) ? v[n/2]
                                : 0.5f * (v[n/2 - 1] + v[n/2]);
        };
        global_.tau_r_ns  = med(good_taur);
        global_.tau_f_ns  = med(good_tauf);
        global_.is_global = true;
        fill_grid(global_, clk_ns);
        has_global_       = true;
    }

    valid_ = true;
    std::cerr << "[PulseTemplateStore] loaded " << path
              << ": " << n_channels_loaded_ << " channels, "
              << n_good_ << " pass gates";
    if (has_global_) {
        std::cerr << "; global τ_r=" << global_.tau_r_ns
                  << "ns τ_f=" << global_.tau_f_ns << "ns";
    } else {
        std::cerr << "; no global template (need ≥5 good channels)";
    }
    std::cerr << "\n";
    return true;
}

const PulseTemplate *
PulseTemplateStore::Lookup(int roc_tag, int slot, int channel,
                           bool fallback_to_global) const
{
    if (!valid_) return nullptr;

    auto it = by_csc_.find(pack_key(roc_tag, slot, channel));
    if (it != by_csc_.end() && it->second.is_good)
        return &it->second.tmpl;

    if (fallback_to_global && has_global_)
        return &global_;

    return nullptr;
}
