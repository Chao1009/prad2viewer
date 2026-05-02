#pragma once
//=============================================================================
// PulseTemplateStore.h — per-channel pulse-template lookup for the NNLS
// pile-up deconvolver.
//
// Loads the JSON written by `analysis/pyscripts/fit_pulse_template.py`,
// validates each per-channel entry against the same gates the analyzer
// would apply (min_pulses / chi2_max / τ_r / τ_f range), synthesises a
// global-median fallback template, and answers (roc_tag, slot, channel)
// → const PulseTemplate* queries.
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
//   * channel absent from JSON   → Lookup() returns nullptr
//   * per-channel gates failed   → Lookup() returns either the global
//     template (if `fallback_to_global_template` AND a global template
//     exists) or nullptr
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

    // Load the per-channel JSON (output of fit_pulse_template.py).
    // Validates each channel against `cfg.nnls_deconv` τ-range gates
    // (and the JSON _meta block's min_pulses / chi2_max, no overrides
    // at this layer); precomputes a unit-amplitude template grid for
    // each good channel at sample period 1000/cfg.clk_mhz so the
    // deconv hot loop can skip per-sample exp().  Returns true on
    // success — `valid()` will then be true and at least one good
    // per-channel template is present.  False on file-not-found,
    // parse error, or empty contents; logs a one-line warning to
    // stderr in any failure path.
    //
    // Re-loading is allowed; the existing contents are dropped first.
    bool LoadFromFile(const std::string &path, const WaveConfig &cfg);

    // Look up the template for a channel identified by its (roc_tag,
    // slot, channel) triple — the same triple the fitter records as
    // `channel_id` "<roc_tag>_<slot>_<channel>".
    //
    // Returns:
    //   * the per-channel template if present AND it passed gates
    //   * the global-median fallback if the per-channel template is
    //     missing / failed gates AND `fallback_to_global` is true AND
    //     a global template exists
    //   * nullptr otherwise (caller should skip deconv)
    //
    // The returned pointer is owned by the store and remains valid
    // until the next LoadFromFile() / store destruction.
    const PulseTemplate *Lookup(int roc_tag, int slot, int channel,
                                bool fallback_to_global) const;

    // True iff LoadFromFile() succeeded and at least one good per-channel
    // template is loaded.  Caller can use this to short-circuit deconv
    // when the store is empty.
    bool valid()      const { return valid_; }
    bool has_global() const { return has_global_; }

    // Diagnostics — counts after LoadFromFile.
    int n_channels_loaded() const { return n_channels_loaded_; }
    int n_good()            const { return n_good_; }

    // Pointer to the synthesised global template, or nullptr if none.
    const PulseTemplate *global_template() const
    {
        return has_global_ ? &global_ : nullptr;
    }

    // Reset to empty state (no templates, invalid).  Useful in tests
    // and for "disable deconv at runtime" code paths.
    void Clear();

private:
    struct Entry {
        PulseTemplate tmpl;
        bool          is_good;
    };

    static uint64_t pack_key(int roc_tag, int slot, int channel)
    {
        return (static_cast<uint64_t>(roc_tag) << 40) |
               (static_cast<uint64_t>(slot)    << 20) |
               static_cast<uint64_t>(channel);
    }

    std::unordered_map<uint64_t, Entry> by_csc_;
    PulseTemplate global_{};
    bool          valid_      = false;
    bool          has_global_ = false;
    int           n_channels_loaded_ = 0;
    int           n_good_     = 0;
};

} // namespace fdec
