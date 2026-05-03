#include "PulseTemplateStore.h"

#include <nlohmann/json.hpp>
#include <fstream>
#include <iostream>
#include <cmath>
#include <cstdio>
#include <cstdlib>

using namespace fdec;
using nlohmann::json;

namespace {

// Pull "<roc>_<slot>_<channel>" out of the entry's `channel_id` string.
// Returns true on success, false if the format doesn't match.
bool parse_channel_id(const std::string &s, int &roc, int &slot, int &ch)
{
    int a = 0, b = 0, c = 0;
    if (std::sscanf(s.c_str(), "%d_%d_%d", &a, &b, &c) != 3) return false;
    roc = a; slot = b; ch = c;
    return true;
}

// Read a {"median": <num>, "mad": <num>} sub-object's median field.
// The fitter writes NaN when no pulses contributed; we propagate NaN.
double read_median(const json &j, const char *key)
{
    if (!j.contains(key)) return std::nan("");
    const auto &sub = j[key];
    if (!sub.is_object() || !sub.contains("median")) return std::nan("");
    const auto &m = sub["median"];
    if (m.is_number()) return m.get<double>();
    return std::nan("");
}

// Fill `tmpl.grid` with the unit-amplitude two-tau template on the
// hardwired 8× oversampled time axis t_i = i · (clk_ns / GRID_OVERSAMPLE).
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

void PulseTemplateStore::Clear()
{
    channel_type_.clear();
    by_type_.clear();
    valid_ = false;
}

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
    const float tr_lo = dcfg.tau_r_min_ns, tr_hi = dcfg.tau_r_max_ns;
    const float tf_lo = dcfg.tau_f_min_ns, tf_hi = dcfg.tau_f_max_ns;

    // ---- channel_id → module_type lookup --------------------------------
    // We only need each channel's category; per-channel τ values are
    // intentionally ignored (the deconvolver uses the per-type aggregate).
    for (auto it = j.begin(); it != j.end(); ++it) {
        const std::string &name = it.key();
        if (!name.empty() && name[0] == '_') continue;     // skip _meta etc.
        const auto &rec = it.value();
        if (!rec.is_object())                  continue;

        if (!rec.contains("channel_id") || !rec["channel_id"].is_string())
            continue;
        int roc = -1, slot = -1, ch = -1;
        if (!parse_channel_id(rec["channel_id"].get<std::string>(),
                              roc, slot, ch))
            continue;

        if (!rec.contains("module_type") || !rec["module_type"].is_string())
            continue;
        std::string type_name = rec["module_type"].get<std::string>();
        if (type_name.empty() || type_name == "Unknown")    continue;

        channel_type_.emplace(pack_key(roc, slot, ch), std::move(type_name));
    }

    // ---- per-type templates from `_by_type` block ----------------------
    if (j.contains("_by_type") && j["_by_type"].is_object()) {
        for (auto bt = j["_by_type"].begin(); bt != j["_by_type"].end(); ++bt) {
            const std::string &type_name = bt.key();
            const auto &rec = bt.value();
            if (!rec.is_object()) continue;
            const double tr = read_median(rec, "tau_r_ns");
            const double tf = read_median(rec, "tau_f_ns");
            if (!std::isfinite(tr) || !std::isfinite(tf)) continue;
            if (tr < tr_lo || tr > tr_hi || tf < tf_lo || tf > tf_hi) continue;
            PulseTemplate t{};
            t.tau_r_ns    = static_cast<float>(tr);
            t.tau_f_ns    = static_cast<float>(tf);
            t.is_global   = true;     // a category aggregate, not a per-channel fit
            fill_grid(t, clk_ns);
            by_type_.emplace(type_name, t);
        }
    }

    if (by_type_.empty()) {
        std::cerr << "[PulseTemplateStore] WARN: " << path
                  << " yielded no per-type templates"
                  << " — deconv will fall back to non-deconv mode.\n";
        return false;
    }

    valid_ = true;
    std::cerr << "[PulseTemplateStore] loaded " << path
              << ": " << channel_type_.size() << " channels typed,"
              << " per-type:";
    for (const auto &kv : by_type_)
        std::cerr << " " << kv.first
                  << "(τ_r=" << kv.second.tau_r_ns
                  << ",τ_f=" << kv.second.tau_f_ns << ")";
    std::cerr << "\n";
    return true;
}

const PulseTemplate *
PulseTemplateStore::type_template(const std::string &type_name) const
{
    auto it = by_type_.find(type_name);
    return (it != by_type_.end()) ? &it->second : nullptr;
}

const PulseTemplate *
PulseTemplateStore::Lookup(int roc_tag, int slot, int channel) const
{
    if (!valid_) return nullptr;

    auto it = channel_type_.find(pack_key(roc_tag, slot, channel));
    if (it == channel_type_.end()) return nullptr;

    auto bt = by_type_.find(it->second);
    return (bt != by_type_.end()) ? &bt->second : nullptr;
}
