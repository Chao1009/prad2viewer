#include "app_state.h"
#include "WaveAnalyzer.h"
#include "HyCalCluster.h"

using json = nlohmann::json;

//=============================================================================
// Filters
//=============================================================================

void AppState::resolveFilterKeys()
{
    filter_wf_keys.clear();
    for (auto &name : waveform_filter.modules) {
        const auto *mod = hycal.module_by_name(name);
        if (!mod || mod->daq.crate < 0) continue;
        for (auto &[roc_tag, crate_id] : roc_to_crate) {
            if (crate_id == mod->daq.crate) {
                filter_wf_keys.insert(std::to_string(roc_tag) + "_"
                    + std::to_string(mod->daq.slot) + "_"
                    + std::to_string(mod->daq.channel));
                break;
            }
        }
    }

    filter_cl_includes.clear();
    for (auto &name : cluster_filter.includes_modules) {
        const auto *mod = hycal.module_by_name(name);
        if (mod) filter_cl_includes.insert(mod->index);
    }

    filter_cl_centers.clear();
    for (auto &name : cluster_filter.center_modules) {
        const auto *mod = hycal.module_by_name(name);
        if (mod) filter_cl_centers.insert(mod->id);
    }
}

std::string AppState::loadFilter(const json &j)
{
    trigger_type_filter = {};
    waveform_filter = {};
    cluster_filter = {};

    // trigger type filter
    if (j.contains("trigger_type")) {
        auto &tt = j["trigger_type"];
        if (tt.contains("enable")) trigger_type_filter.enable = tt["enable"];
        if (tt.contains("accept") && tt["accept"].is_array()) {
            for (auto &v : tt["accept"])
                trigger_type_filter.accept.push_back(static_cast<uint8_t>(v.get<int>()));
        }
    }

    if (j.contains("waveform")) {
        auto &w = j["waveform"];
        auto &f = waveform_filter;
        if (w.contains("enable"))       f.enable       = w["enable"];
        if (w.contains("n_peaks_min"))  f.n_peaks_min  = w["n_peaks_min"];
        if (w.contains("n_peaks_max"))  f.n_peaks_max  = w["n_peaks_max"];
        if (w.contains("time_min"))     f.time_min     = w["time_min"];
        if (w.contains("time_max"))     f.time_max     = w["time_max"];
        if (w.contains("integral_min")) f.integral_min = w["integral_min"];
        if (w.contains("integral_max")) f.integral_max = w["integral_max"];
        if (w.contains("height_min"))   f.height_min   = w["height_min"];
        if (w.contains("height_max"))   f.height_max   = w["height_max"];
        if (w.contains("modules"))
            for (auto &m : w["modules"])
                f.modules.push_back(m.get<std::string>());
    }

    if (j.contains("clustering")) {
        auto &c = j["clustering"];
        auto &f = cluster_filter;
        if (c.contains("enable"))           f.enable       = c["enable"];
        if (c.contains("n_min"))            f.n_min        = c["n_min"];
        if (c.contains("n_max"))            f.n_max        = c["n_max"];
        if (c.contains("energy_min"))       f.energy_min   = c["energy_min"];
        if (c.contains("energy_max"))       f.energy_max   = c["energy_max"];
        if (c.contains("size_min"))         f.size_min     = c["size_min"];
        if (c.contains("size_max"))         f.size_max     = c["size_max"];
        if (c.contains("includes_min"))     f.includes_min = c["includes_min"];
        if (c.contains("includes_modules"))
            for (auto &m : c["includes_modules"])
                f.includes_modules.push_back(m.get<std::string>());
        if (c.contains("center_modules"))
            for (auto &m : c["center_modules"])
                f.center_modules.push_back(m.get<std::string>());
    }

    resolveFilterKeys();

    std::cerr << "Filter loaded: waveform "
              << (waveform_filter.enable ? "ON" : "off")
              << " (" << filter_wf_keys.size() << " modules)"
              << ", cluster "
              << (cluster_filter.enable ? "ON" : "off") << "\n";
    return "";
}

void AppState::unloadFilter()
{
    trigger_type_filter = {};
    waveform_filter = {};
    cluster_filter = {};
    filter_wf_keys.clear();
    filter_cl_includes.clear();
    filter_cl_centers.clear();
    std::cerr << "Filters unloaded\n";
}

json AppState::filterToJson() const
{
    json r;
    r["active"] = filterActive();
    {
        auto &f = waveform_filter;
        json w;
        w["enable"] = f.enable;
        if (!f.modules.empty()) w["modules"] = f.modules;
        w["n_peaks_min"] = f.n_peaks_min;
        w["n_peaks_max"] = f.n_peaks_max;
        if (f.time_min > -1e20f)     w["time_min"]     = f.time_min;
        if (f.time_max <  1e20f)     w["time_max"]     = f.time_max;
        if (f.integral_min > -1e20f) w["integral_min"] = f.integral_min;
        if (f.integral_max <  1e20f) w["integral_max"] = f.integral_max;
        if (f.height_min > -1e20f)   w["height_min"]   = f.height_min;
        if (f.height_max <  1e20f)   w["height_max"]   = f.height_max;
        r["waveform"] = w;
    }
    {
        auto &f = cluster_filter;
        json c;
        c["enable"] = f.enable;
        c["n_min"] = f.n_min;
        c["n_max"] = f.n_max;
        if (f.energy_min > 0)    c["energy_min"] = f.energy_min;
        if (f.energy_max < 1e20f) c["energy_max"] = f.energy_max;
        c["size_min"] = f.size_min;
        c["size_max"] = f.size_max;
        if (!f.includes_modules.empty()) {
            c["includes_modules"] = f.includes_modules;
            c["includes_min"] = f.includes_min;
        }
        if (!f.center_modules.empty()) c["center_modules"] = f.center_modules;
        r["clustering"] = c;
    }
    return r;
}

bool AppState::evaluateFilter(fdec::EventData &event,
                              ssp::SspEventData *ssp) const
{
    if (!filterActive()) return true;

    // --- trigger type filter (fast, check first) ---
    if (trigger_type_filter.enable && !trigger_type_filter(event.info.trigger_type))
        return false;

    bool is_adc1881m = (daq_cfg.adc_format == "adc1881m");

    // --- waveform filter ---
    if (waveform_filter.enable) {
        if (filter_wf_keys.empty()) return false;  // enabled but no modules resolved

        bool any_module_pass = false;
        fdec::WaveAnalyzer ana(daq_cfg.wave_cfg);
        ana.cfg.min_peak_ratio = hist_cfg.min_peak_ratio;
        ana.SetTemplateStore(&template_store);
        fdec::WaveResult wres;

        for (int r = 0; r < event.nrocs && !any_module_pass; ++r) {
            auto &roc = event.rocs[r];
            if (!roc.present) continue;
            for (int s = 0; s < fdec::MAX_SLOTS && !any_module_pass; ++s) {
                if (!roc.slots[s].present) continue;
                auto &slot = roc.slots[s];
                for (int c = 0; c < fdec::MAX_CHANNELS && !any_module_pass; ++c) {
                    if (!(slot.channel_mask & (1ull << c))) continue;
                    auto &cd = slot.channels[c];
                    if (cd.nsamples <= 0) continue;

                    std::string key = std::to_string(roc.tag) + "_"
                                    + std::to_string(s) + "_" + std::to_string(c);
                    if (filter_wf_keys.find(key) == filter_wf_keys.end()) continue;

                    if (is_adc1881m) continue;  // no peak analysis for ADC1881M

                    ana.SetChannelKey(roc.tag, s, c);
                    ana.Analyze(cd.samples, cd.nsamples, wres);
                    int n_qual = 0;
                    for (int p = 0; p < wres.npeaks; ++p) {
                        auto &pk = wres.peaks[p];
                        if (pk.height < waveform_filter.height_min ||
                            pk.height > waveform_filter.height_max) continue;
                        if (pk.time < waveform_filter.time_min ||
                            pk.time > waveform_filter.time_max) continue;
                        if (pk.integral < waveform_filter.integral_min ||
                            pk.integral > waveform_filter.integral_max) continue;
                        n_qual++;
                    }
                    if (n_qual >= waveform_filter.n_peaks_min &&
                        n_qual <= waveform_filter.n_peaks_max)
                        any_module_pass = true;
                }
            }
        }
        if (!any_module_pass) return false;
    }

    // --- cluster filter ---
    if (cluster_filter.enable) {
        // build clusters from this event (local clusterer, no side effects)
        fdec::HyCalCluster clusterer(hycal);
        clusterer.SetConfig(cluster_cfg);

        fdec::WaveAnalyzer ana(daq_cfg.wave_cfg);
        ana.cfg.min_peak_ratio = hist_cfg.min_peak_ratio;
        ana.SetTemplateStore(&template_store);
        fdec::WaveResult wres;

        for (int r = 0; r < event.nrocs; ++r) {
            auto &roc = event.rocs[r];
            if (!roc.present) continue;
            auto cit = roc_to_crate.find(roc.tag);
            if (cit == roc_to_crate.end()) continue;
            int crate = cit->second;

            for (int s = 0; s < fdec::MAX_SLOTS; ++s) {
                if (!roc.slots[s].present) continue;
                auto &slot = roc.slots[s];
                for (int c = 0; c < fdec::MAX_CHANNELS; ++c) {
                    if (!(slot.channel_mask & (1ull << c))) continue;
                    auto &cd = slot.channels[c];
                    if (cd.nsamples <= 0) continue;

                    const auto *mod = hycal.module_by_daq(crate, s, c);
                    if (!mod || !mod->is_hycal()) continue;

                    float adc_val;
                    if (is_adc1881m) {
                        adc_val = (float)cd.samples[0];
                    } else {
                        ana.SetChannelKey(roc.tag, s, c);
                        ana.Analyze(cd.samples, cd.nsamples, wres);
                        adc_val = bestPeakInWindow(wres, hist_cfg.threshold,
                                                   hist_cfg.time_min, hist_cfg.time_max);
                    }
                    if (adc_val > 0) {
                        float energy = (mod->cal_factor > 0.)
                            ? static_cast<float>(mod->energize(adc_val))
                            : adc_val * adc_to_mev;
                        clusterer.AddHit(mod->index, energy, 0.f);
                    }
                }
            }
        }

        clusterer.FormClusters();
        std::vector<fdec::HyCalCluster::RecoResult> results;
        clusterer.ReconstructMatched(results);

        int n_qual = 0;
        for (auto &rr : results) {
            auto &hit = rr.hit;
            auto *cl = rr.cluster;

            // energy range
            if (hit.energy < cluster_filter.energy_min ||
                hit.energy > cluster_filter.energy_max) continue;
            // size range
            if (hit.nblocks < cluster_filter.size_min ||
                hit.nblocks > cluster_filter.size_max) continue;
            // includes check
            if (!filter_cl_includes.empty()) {
                int count = 0;
                for (auto &h : cl->hits)
                    if (filter_cl_includes.count(h.index)) count++;
                if (count < cluster_filter.includes_min) continue;
            }
            // center check
            if (!filter_cl_centers.empty()) {
                if (!filter_cl_centers.count(hit.center_id)) continue;
            }
            n_qual++;
        }
        if (n_qual < cluster_filter.n_min || n_qual > cluster_filter.n_max)
            return false;
    }

    return true;
}
