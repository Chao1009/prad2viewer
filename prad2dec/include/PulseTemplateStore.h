#pragma once
//=============================================================================
// PulseTemplateStore.h — per-type pulse-template lookup for the pile-up
// deconvolver.
//
// Loads the JSON written by `analysis/pyscripts/fit_pulse_template.py`
// and exposes the per-type aggregates (PbGlass / PbWO4 / LMS / Veto)
// from its `_by_type` block.  Per-channel τ_r / τ_f entries in the JSON
// are read only to learn each channel's module type — the per-channel
// shapes themselves are no longer used.  The deconvolver gets one
// shape per category, which keeps the inputs well-conditioned (a
// well-fit median across many channels) and avoids amplifying noise
// from low-statistics single-channel fits.
//
// The store is owned by the application layer (Replay, viewer servers,
// gain tools) and bound to a WaveAnalyzer via
// `WaveAnalyzer::SetTemplateStore(&store)` once at setup; the analyzer
// itself never touches a JSON file.
//
// Failure modes are deliberately quiet so production code can keep
// running without templates:
//
//   * file missing / parse error → warning to stderr, store stays
//     invalid()  ⇒ Lookup() returns nullptr always
//   * channel's type is unknown  → Lookup() returns nullptr
//   * per-type entry rejected at load time (τ outside config range) →
//     Lookup() returns nullptr for that type's channels
//=============================================================================

#include "Fadc250Data.h"
#include "WaveAnalyzer.h"   // for WaveConfig::NnlsDeconvConfig

#include <cstdint>
#include <string>
#include <unordered_map>

namespace fdec
{

// Thread safety: after a successful LoadFromFile() the store is
// effectively immutable and Lookup() is safe to call concurrently from
// any number of threads.  LoadFromFile() and Clear() are NOT
// concurrent-safe — invoke them only from a single thread when no
// WaveAnalyzer is using the store.
class PulseTemplateStore
{
public:
    PulseTemplateStore() = default;

    // Load the per-type templates from the JSON written by
    // `fit_pulse_template.py`.  `cfg.nnls_deconv` provides the τ-range
    // gates each per-type entry must satisfy to be accepted, and
    // `cfg.clk_mhz` sets the precomputed grid sample period so the
    // deconv hot loop can skip per-sample exp().  Returns true on
    // success — `valid()` will then be true and at least one per-type
    // template is present.  False on file-not-found, parse error, or
    // empty contents; logs a one-line warning to stderr in any failure
    // path.
    //
    // Re-loading is allowed; the existing contents are dropped first.
    bool LoadFromFile(const std::string &path, const WaveConfig &cfg);

    // Look up the deconv template for a channel identified by its
    // (roc_tag, slot, channel) triple — the same triple the fitter
    // records as `channel_id` "<roc_tag>_<slot>_<channel>".  The store
    // resolves the channel's module type from its loaded per-channel
    // metadata and returns the matching per-type template, or nullptr
    // if the channel is unknown / its type lacks a usable per-type
    // entry.
    //
    // The returned pointer is owned by the store and remains valid
    // until the next LoadFromFile() / store destruction.
    const PulseTemplate *Lookup(int roc_tag, int slot, int channel) const;

    // Pointer to the per-type template (e.g. type_name = "PbGlass" /
    // "PbWO4" / "LMS" / "Veto"), or nullptr if the file's `_by_type`
    // block did not include this type / it failed gates.  Exposed for
    // diagnostics and Python scripts that want the per-type template
    // directly without going through (roc, slot, channel).
    const PulseTemplate *type_template(const std::string &type_name) const;

    // True iff LoadFromFile() succeeded and at least one per-type
    // template is loaded.  Callers can short-circuit deconv when this
    // returns false.
    bool valid() const { return valid_; }

    // Diagnostics.
    int  n_channels_known() const { return static_cast<int>(channel_type_.size()); }
    int  n_types_loaded()   const { return static_cast<int>(by_type_.size()); }

    // Reset to empty state (no templates, invalid).  Useful in tests
    // and for "disable deconv at runtime" code paths.
    void Clear();

private:
    static uint64_t pack_key(int roc_tag, int slot, int channel)
    {
        return (static_cast<uint64_t>(roc_tag) << 40) |
               (static_cast<uint64_t>(slot)    << 20) |
               static_cast<uint64_t>(channel);
    }

    // (roc_tag, slot, channel) → "PbGlass" / "PbWO4" / "LMS" / "Veto".
    // Built from each per-channel entry's `module_type` field at load
    // time; per-channel τ values are intentionally not stored.
    std::unordered_map<uint64_t, std::string>      channel_type_;

    // Per-type templates from the file's `_by_type` block — keyed by
    // the type name string the fitter wrote ("PbGlass", "PbWO4",
    // "LMS", "Veto").  Each entry only present when the JSON carried
    // both τ_r and τ_f for that type AND they passed the analyzer's
    // tau-range gates.
    std::unordered_map<std::string, PulseTemplate> by_type_;

    bool valid_ = false;
};

} // namespace fdec
