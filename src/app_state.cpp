#include "app_state.h"
#include "data_source.h"
#include "load_daq_config.h"

#include <fstream>
#include <iostream>
#include <cmath>
#include <cstdlib>

using json = nlohmann::json;

static json histToJson(const Histogram &h, float mn, float mx, float st)
{
    if (h.bins.empty())
        return {{"bins", json::array()}, {"underflow", 0}, {"overflow", 0},
                {"min", mn}, {"max", mx}, {"step", st}};
    return {{"bins", h.bins}, {"underflow", h.underflow}, {"overflow", h.overflow},
            {"min", mn}, {"max", mx}, {"step", st}};
}

//=============================================================================
// Per-event processing
//=============================================================================

// Encode peak array for one channel.
static json encodePeaks(const fdec::WaveResult &wres)
{
    json parr = json::array();
    for (int p = 0; p < wres.npeaks; ++p) {
        auto &pk = wres.peaks[p];
        parr.push_back({
            {"p", pk.pos}, {"t", std::round(pk.time * 10) / 10},
            {"h", std::round(pk.height * 10) / 10},
            {"i", std::round(pk.integral * 10) / 10},
            {"l", pk.left}, {"r", pk.right},
            {"o", pk.overflow ? 1 : 0},
        });
    }
    return parr;
}

json AppState::encodeEventJson(fdec::EventData &event, int ev_id,
                               fdec::WaveAnalyzer &ana, fdec::WaveResult &wres,
                               bool include_samples)
{
    json channels = json::object();
    for (int r = 0; r < event.nrocs; ++r) {
        auto &roc = event.rocs[r];
        if (!roc.present) continue;
        for (int s = 0; s < fdec::MAX_SLOTS; ++s) {
            if (!roc.slots[s].present) continue;
            auto &slot = roc.slots[s];
            for (int c = 0; c < fdec::MAX_CHANNELS; ++c) {
                if (!(slot.channel_mask & (1ull << c))) continue;
                auto &cd = slot.channels[c];
                if (cd.nsamples <= 0) continue;

                ana.Analyze(cd.samples, cd.nsamples, wres);
                std::string key = std::to_string(roc.tag) + "_"
                                + std::to_string(s) + "_" + std::to_string(c);

                json ch_j = {
                    {"pm", std::round(wres.ped.mean * 10) / 10},
                    {"pr", std::round(wres.ped.rms * 10) / 10},
                    {"pk", encodePeaks(wres)},
                };
                if (include_samples) {
                    json sarr = json::array();
                    for (int j = 0; j < cd.nsamples; ++j) sarr.push_back(cd.samples[j]);
                    ch_j["s"] = std::move(sarr);
                }
                channels[key] = std::move(ch_j);
            }
        }
    }
    return {{"event", ev_id}, {"channels", channels},
            {"event_number", event.info.event_number},
            {"trigger_type", event.info.trigger_type},
            {"trigger_bits", event.info.trigger_bits}};
}

json AppState::encodeWaveformJson(fdec::EventData &event, const std::string &chan_key,
                                  fdec::WaveAnalyzer &ana, fdec::WaveResult &wres)
{
    // parse "roc_slot_ch" key
    int roc_tag = 0, sl = 0, ch = 0;
    if (std::sscanf(chan_key.c_str(), "%d_%d_%d", &roc_tag, &sl, &ch) != 3)
        return {{"error", "invalid channel key"}};

    // find the channel in the event
    for (int r = 0; r < event.nrocs; ++r) {
        auto &roc = event.rocs[r];
        if (!roc.present || roc.tag != roc_tag) continue;
        if (!roc.slots[sl].present) break;
        if (!(roc.slots[sl].channel_mask & (1ull << ch))) break;
        auto &cd = roc.slots[sl].channels[ch];
        if (cd.nsamples <= 0) break;

        ana.Analyze(cd.samples, cd.nsamples, wres);

        json sarr = json::array();
        for (int j = 0; j < cd.nsamples; ++j) sarr.push_back(cd.samples[j]);

        return {{"key", chan_key}, {"s", sarr},
                {"pm", std::round(wres.ped.mean * 10) / 10},
                {"pr", std::round(wres.ped.rms * 10) / 10},
                {"pk", encodePeaks(wres)}};
    }
    return {{"error", "channel not found"}};
}

void AppState::projectToHyCalLocal(float Gx, float Gy, float Gz,
                                   float &px, float &py) const
{
    // Transform target and source points into HyCal-local frame, then linearly
    // interpolate along the line to the local z=0 plane.  Equivalent to
    // intersecting the lab-frame line with the tilted HyCal plane, but cleaner
    // because the math reduces to a 1D parametric solve.
    float Tx, Ty, Tz, gx, gy, gz;
    hycal_transform.labToLocal(target_x, target_y, target_z, Tx, Ty, Tz);
    hycal_transform.labToLocal(Gx, Gy, Gz, gx, gy, gz);
    float dz = gz - Tz;
    if (std::abs(dz) < 1e-6f) { px = gx; py = gy; return; }
    float s = -Tz / dz;
    px = Tx + s * (gx - Tx);
    py = Ty + s * (gy - Ty);
}

json AppState::computeClustersJson(fdec::EventData &event, int ev_id,
                                   fdec::WaveAnalyzer &ana, fdec::WaveResult &wres)
{
    if (!cluster_trigger(event.info.trigger_bits))
        return {{"event", ev_id}, {"hits", json::object()}, {"clusters", json::array()},
                {"info", "trigger filtered"}};

    bool is_adc1881m = (daq_cfg.adc_format == "adc1881m");
    fdec::HyCalCluster clusterer(hycal);
    clusterer.SetConfig(cluster_cfg);

    int nmod = hycal.module_count();
    std::vector<float> mod_energy(nmod, 0.f);

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

                float adc_val = 0;
                if (is_adc1881m) {
                    adc_val = cd.samples[0];
                } else {
                    ana.Analyze(cd.samples, cd.nsamples, wres);
                    adc_val = bestPeakInWindow(wres, hist_cfg.threshold,
                                               hist_cfg.time_min, hist_cfg.time_max);
                }
                if (adc_val <= 0) continue;

                float energy = (mod->cal_factor > 0.)
                    ? static_cast<float>(mod->energize(adc_val))
                    : adc_val * adc_to_mev;
                mod_energy[mod->index] = energy;
                clusterer.AddHit(mod->index, energy);
            }
        }
    }

    clusterer.FormClusters();

    json hits_j = json::object();
    for (int i = 0; i < nmod; ++i)
        if (mod_energy[i] > 0.f)
            hits_j[std::to_string(i)] = std::round(mod_energy[i] * 100) / 100;

    std::vector<fdec::HyCalCluster::RecoResult> reco;
    clusterer.ReconstructMatched(reco);

    json cl_arr = json::array();
    for (auto &r : reco) {
        auto &cmod = hycal.module(r.cluster->center.index);
        json indices = json::array();
        for (auto &h : r.cluster->hits) indices.push_back(h.index);
        cl_arr.push_back({
            {"id", static_cast<int>(cl_arr.size())},
            {"center", cmod.name}, {"center_id", cmod.id},
            {"x", std::round(r.hit.x * 10) / 10},
            {"y", std::round(r.hit.y * 10) / 10},
            {"energy", std::round(r.hit.energy * 10) / 10},
            {"nblocks", r.hit.nblocks}, {"npos", r.hit.npos},
            {"modules", indices},
        });
    }

    return {{"event", ev_id}, {"hits", hits_j}, {"clusters", cl_arr}};
}

void AppState::recordSyncTime(uint32_t unix_time, uint64_t last_ti_ts)
{
    if (unix_time == 0) return;
    std::lock_guard<std::mutex> lk(lms_mtx);
    if (sync_unix != 0) return;   // already have a sync reference

    if (lms_first_ts == 0) {
        // No LMS events yet — stash for later.
        // Will be applied when the first LMS event sets lms_first_ts.
        pending_sync_unix = unix_time;
        pending_sync_ti = last_ti_ts;
        return;
    }

    sync_unix = unix_time;
    sync_rel_sec = (last_ti_ts != 0)
        ? static_cast<double>(last_ti_ts - lms_first_ts) * TI_TICK_SEC
        : 0.;
}

void AppState::processEvent(fdec::EventData &event,
                            fdec::WaveAnalyzer &ana, fdec::WaveResult &wres)
{
    // --- check which consumers need this event ---
    uint32_t tb = event.info.trigger_bits;
    bool do_hist    = waveform_trigger(tb);
    bool do_cluster = cluster_trigger(tb);
    bool do_lms     = lms_trigger.accept != 0 && lms_trigger(tb);
    bool do_alpha   = alpha_trigger.accept != 0 && alpha_trigger(tb);

    if (!do_hist && !do_cluster && !do_lms && !do_alpha) {
        std::lock_guard<std::mutex> lk(data_mtx);
        events_processed++;
        return;
    }

    bool is_adc1881m = (daq_cfg.adc_format == "adc1881m");

    // clustering setup (stack-allocated, per-event)
    fdec::HyCalCluster clusterer(hycal);
    if (do_cluster) clusterer.SetConfig(cluster_cfg);

    // LMS timing
    double lms_time = 0;

    // acquire both locks for the merged pass
    std::unique_lock<std::mutex> lk1(data_mtx, std::defer_lock);
    std::unique_lock<std::mutex> lk2(lms_mtx, std::defer_lock);
    std::lock(lk1, lk2);

    if (do_lms) {
        if (lms_first_ts == 0) {
            lms_first_ts = event.info.timestamp;
            // apply stashed sync time from a control event that arrived before LMS data
            if (pending_sync_unix != 0 && sync_unix == 0) {
                sync_unix = pending_sync_unix;
                // PRESTART/GO arrives before physics events, so pending_sync_ti is
                // typically 0. In that case sync_rel_sec = 0 (run start = LMS start).
                sync_rel_sec = (pending_sync_ti != 0)
                    ? static_cast<double>(pending_sync_ti - lms_first_ts) * TI_TICK_SEC
                    : 0.;
            }
        }
        lms_time = static_cast<double>(event.info.timestamp - lms_first_ts) * TI_TICK_SEC;
    }

    // --- single pass: analyze once per channel, feed all consumers ---
    for (int r = 0; r < event.nrocs; ++r) {
        auto &roc = event.rocs[r];
        if (!roc.present) continue;

        // crate lookup (needed by cluster + LMS + Alpha consumers)
        int crate = -1;
        if (do_cluster || do_lms || do_alpha) {
            auto cit = roc_to_crate.find(roc.tag);
            if (cit != roc_to_crate.end()) crate = cit->second;
        }

        for (int s = 0; s < fdec::MAX_SLOTS; ++s) {
            if (!roc.slots[s].present) continue;
            auto &slot = roc.slots[s];
            for (int c = 0; c < fdec::MAX_CHANNELS; ++c) {
                if (!(slot.channel_mask & (1ull << c))) continue;
                auto &cd = slot.channels[c];
                if (cd.nsamples <= 0) continue;

                // ── analyze ONCE ──
                float peak_in_window = -1;
                if (!is_adc1881m) {
                    ana.Analyze(cd.samples, cd.nsamples, wres);
                    peak_in_window = bestPeakInWindow(wres, hist_cfg.threshold,
                                                       hist_cfg.time_min, hist_cfg.time_max);
                } else {
                    wres.npeaks = 0;
                    peak_in_window = cd.samples[0];
                }

                // ── histogram consumer ──
                if (do_hist && !is_adc1881m) {
                    std::string key = std::to_string(roc.tag) + "_"
                                   + std::to_string(s) + "_" + std::to_string(c);
                    bool has_peak = false, has_peak_tcut = false;
                    float best = -1, best_height = -1;
                    for (int p = 0; p < wres.npeaks; ++p) {
                        auto &pk = wres.peaks[p];
                        if (pk.height < hist_cfg.threshold) continue;
                        has_peak = true;
                        if (pk.time >= hist_cfg.time_min && pk.time <= hist_cfg.time_max) {
                            has_peak_tcut = true;
                            if (pk.integral > best) { best = pk.integral; best_height = pk.height; }
                        }
                    }
                    if (best >= 0) {
                        auto &h = histograms[key];
                        if (h.bins.empty()) h.init(hist_nbins);
                        h.fill(best, hist_cfg.bin_min, hist_cfg.bin_step);
                        auto &hh = height_histograms[key];
                        if (hh.bins.empty()) hh.init(height_nbins);
                        hh.fill(best_height, hist_cfg.height_min, hist_cfg.height_step);
                    }
                    for (int p = 0; p < wres.npeaks; ++p) {
                        auto &pk = wres.peaks[p];
                        if (pk.height < hist_cfg.threshold) continue;
                        auto &ph = pos_histograms[key];
                        if (ph.bins.empty()) ph.init(pos_nbins);
                        ph.fill(pk.time, hist_cfg.pos_min, hist_cfg.pos_step);
                    }
                    if (has_peak)      occupancy[key]++;
                    if (has_peak_tcut) occupancy_tcut[key]++;
                }

                // ── cluster consumer ──
                if (do_cluster && crate >= 0) {
                    const auto *mod = hycal.module_by_daq(crate, s, c);
                    if (mod && mod->is_hycal()) {
                        float adc_val = is_adc1881m ? (float)cd.samples[0] : peak_in_window;
                        if (adc_val > 0) {
                            float energy = (mod->cal_factor > 0.)
                                ? static_cast<float>(mod->energize(adc_val))
                                : adc_val * adc_to_mev;
                            clusterer.AddHit(mod->index, energy);
                        }
                    }
                }

                // ── LMS consumer ──
                if (do_lms && crate >= 0) {
                    const auto *mod = hycal.module_by_daq(crate, s, c);
                    if (mod) {
                        float val = is_adc1881m ? (float)cd.samples[0] : peak_in_window;
                        if (val > 0) {
                            auto &hist = lms_history[mod->index];
                            if (static_cast<int>(hist.size()) < lms_max_history)
                                hist.push_back({lms_time, val});
                            // Always track the latest reading, even after history saturates,
                            // so the LMS/Alpha ref correction stays current.
                            latest_lms_integral[mod->index] = val;
                        }
                    }
                }

                // ── Alpha consumer (Am-241 reference) ──
                if (do_alpha && crate >= 0) {
                    const auto *mod = hycal.module_by_daq(crate, s, c);
                    if (mod) {
                        float val = is_adc1881m ? (float)cd.samples[0] : peak_in_window;
                        if (val > 0) latest_alpha_integral[mod->index] = val;
                    }
                }
            }
        }
    }

    // --- post-loop: clustering + physics histograms ---
    if (do_cluster) {
        clusterer.FormClusters();
        std::vector<fdec::ClusterHit> reco_hits;
        clusterer.ReconstructHits(reco_hits);

        struct ClusterInfo { float lx, ly, lz, theta; };
        std::vector<ClusterInfo> cinfo(reco_hits.size());
        for (size_t i = 0; i < reco_hits.size(); ++i) {
            auto &rh = reco_hits[i];
            auto &ci = cinfo[i];
            hycal_transform.toLab(rh.x, rh.y, ci.lx, ci.ly, ci.lz);
            float dx = ci.lx - target_x, dy = ci.ly - target_y, dz = ci.lz - target_z;
            float rv = std::sqrt(dx*dx + dy*dy);
            ci.theta = std::atan2(rv, dz) * (180.f / 3.14159265f);
        }

        // Per-Ncl bucket index (-1 if Ncl falls outside the nclusters_hist
        // range, in which case the bucketed hists get no fill — same
        // semantics as Histogram::fill underflow/overflow).
        int ncl_bucket = -1;
        {
            float fb = ((float)reco_hits.size() - nclusters_hist_min) / nclusters_hist_step;
            if (fb >= 0.f) {
                int b = (int)fb;
                if (b < (int)cluster_energy_hist_by_ncl.size())
                    ncl_bucket = b;
            }
        }
        for (size_t i = 0; i < reco_hits.size(); ++i) {
            cluster_energy_hist.fill(reco_hits[i].energy, cl_hist_min, cl_hist_step);
            nblocks_hist.fill(reco_hits[i].nblocks, nblocks_hist_min, nblocks_hist_step);
            if (ncl_bucket >= 0) {
                cluster_energy_hist_by_ncl[ncl_bucket].fill(
                    reco_hits[i].energy, cl_hist_min, cl_hist_step);
                nblocks_hist_by_ncl[ncl_bucket].fill(
                    reco_hits[i].nblocks, nblocks_hist_min, nblocks_hist_step);
            }
        }
        nclusters_hist.fill(reco_hits.size(), nclusters_hist_min, nclusters_hist_step);
        cluster_events_processed++;

        bool physics_accept = physics_trigger(tb);
        if (physics_accept) {
            float Eb = beam_energy.load();
            for (size_t i = 0; i < reco_hits.size(); ++i) {
                energy_angle_hist.fill(cinfo[i].theta, reco_hits[i].energy,
                    ea_angle_min, ea_angle_step, ea_energy_min, ea_energy_step);
            }
            if (reco_hits.size() == 2 && Eb > 0) {
                float esum = reco_hits[0].energy + reco_hits[1].energy;
                bool energy_ok = std::abs(esum - Eb) < moller_energy_tol * Eb;
                bool angle_ok = false;
                for (int j = 0; j < 2; ++j)
                    if (cinfo[j].theta >= moller_angle_min && cinfo[j].theta <= moller_angle_max)
                        angle_ok = true;
                if (energy_ok && angle_ok) {
                    moller_events++;
                    for (int j = 0; j < 2; ++j) {
                        moller_xy_hist.fill(cinfo[j].lx, cinfo[j].ly,
                            moller_xy_x_min, moller_xy_x_step, moller_xy_y_min, moller_xy_y_step);
                    }
                }
            }
            // HyCal cluster-hit XY: single-cluster ep-elastic candidates.
            // Same gate is reused (when configured) for GEM↔HyCal residuals.
            bool ep_cand = false;
            if ((int)reco_hits.size() == hxy_n_clusters && Eb > 0) {
                const auto &cl = reco_hits[0];
                bool nb_ok = cl.nblocks >= hxy_nblocks_min && cl.nblocks <= hxy_nblocks_max;
                bool e_ok  = cl.energy  >= hxy_energy_frac_min * Eb;
                if (nb_ok && e_ok) {
                    ep_cand = true;
                    hycal_xy_hist.fill(cinfo[0].lx, cinfo[0].ly,
                        hxy_x_min, hxy_x_step, hxy_y_min, hxy_y_step);
                    hycal_xy_events++;
                }
            }
            // GEM↔HyCal matching residuals.  Reference is the FIRST cluster's
            // HyCal-local xy — for ep candidates that's the only cluster, for
            // multi-cluster events it's the leading reconstructed hit.
            if (gem_enabled && (ep_cand || !gem_match_require_ep) && !reco_hits.empty()) {
                const float ref_x = reco_hits[0].x, ref_y = reco_hits[0].y;
                const int n_dets = std::min<int>(gem_sys.GetNDetectors(),
                                                 (int)gem_dx_hist.size());
                for (int d = 0; d < n_dets; ++d) {
                    auto &xform = gem_transforms[d];
                    for (auto &h : gem_sys.GetHits(d)) {
                        float lx, ly, lz;
                        xform.toLab(h.x, h.y, lx, ly, lz);
                        float px, py;
                        projectToHyCalLocal(lx, ly, lz, px, py);
                        float dxr = px - ref_x, dyr = py - ref_y;
                        if (std::sqrt(dxr*dxr + dyr*dyr) < gem_match_window_mm) {
                            gem_dx_hist[d].fill(dxr, gem_resid_min, gem_resid_step);
                            gem_dy_hist[d].fill(dyr, gem_resid_min, gem_resid_step);
                            gem_match_hits[d]++;
                        }
                    }
                }
                gem_match_events++;
            }
        }
    }

    events_processed++;
    if (do_lms) lms_events++;
}

void AppState::processReconEvent(const ReconEventData &recon)
{
    uint32_t tb = recon.trigger_bits;
    bool do_cluster = cluster_trigger(tb);
    bool do_physics = physics_trigger(tb);

    std::lock_guard<std::mutex> lk(data_mtx);
    events_processed++;

    if (do_cluster && !recon.clusters.empty()) {
        int ncl_bucket = -1;
        {
            float fb = ((float)recon.clusters.size() - nclusters_hist_min)
                       / nclusters_hist_step;
            if (fb >= 0.f) {
                int b = (int)fb;
                if (b < (int)cluster_energy_hist_by_ncl.size())
                    ncl_bucket = b;
            }
        }
        for (auto &cl : recon.clusters) {
            cluster_energy_hist.fill(cl.energy, cl_hist_min, cl_hist_step);
            nblocks_hist.fill(cl.nblocks, nblocks_hist_min, nblocks_hist_step);
            if (ncl_bucket >= 0) {
                cluster_energy_hist_by_ncl[ncl_bucket].fill(
                    cl.energy, cl_hist_min, cl_hist_step);
                nblocks_hist_by_ncl[ncl_bucket].fill(
                    cl.nblocks, nblocks_hist_min, nblocks_hist_step);
            }
        }
        nclusters_hist.fill(recon.clusters.size(), nclusters_hist_min, nclusters_hist_step);
        cluster_events_processed++;
    }

    if (do_physics && !recon.clusters.empty()) {
        struct CI { float lx, ly, lz, theta; };
        std::vector<CI> cinfo(recon.clusters.size());
        for (size_t i = 0; i < recon.clusters.size(); ++i) {
            auto &cl = recon.clusters[i];
            auto &ci = cinfo[i];
            hycal_transform.toLab(cl.x, cl.y, ci.lx, ci.ly, ci.lz);
            float dx = ci.lx - target_x, dy = ci.ly - target_y, dz = ci.lz - target_z;
            float r = std::sqrt(dx*dx + dy*dy);
            ci.theta = std::atan2(r, dz) * (180.f / 3.14159265f);
        }
        for (size_t i = 0; i < recon.clusters.size(); ++i)
            energy_angle_hist.fill(cinfo[i].theta, recon.clusters[i].energy,
                ea_angle_min, ea_angle_step, ea_energy_min, ea_energy_step);

        float Eb = beam_energy.load();
        if (recon.clusters.size() == 2 && Eb > 0) {
            float esum = recon.clusters[0].energy + recon.clusters[1].energy;
            bool energy_ok = std::abs(esum - Eb) < moller_energy_tol * Eb;
            bool angle_ok = false;
            for (int j = 0; j < 2; ++j)
                if (cinfo[j].theta >= moller_angle_min && cinfo[j].theta <= moller_angle_max)
                    angle_ok = true;
            if (energy_ok && angle_ok) {
                moller_events++;
                for (int j = 0; j < 2; ++j) {
                    moller_xy_hist.fill(cinfo[j].lx, cinfo[j].ly,
                        moller_xy_x_min, moller_xy_x_step, moller_xy_y_min, moller_xy_y_step);
                }
            }
        }
        // HyCal cluster-hit XY: single-cluster ep-elastic candidates.
        bool ep_cand = false;
        if ((int)recon.clusters.size() == hxy_n_clusters && Eb > 0) {
            const auto &cl = recon.clusters[0];
            bool nb_ok = cl.nblocks >= hxy_nblocks_min && cl.nblocks <= hxy_nblocks_max;
            bool e_ok  = cl.energy  >= hxy_energy_frac_min * Eb;
            if (nb_ok && e_ok) {
                ep_cand = true;
                hycal_xy_hist.fill(cinfo[0].lx, cinfo[0].ly,
                    hxy_x_min, hxy_x_step, hxy_y_min, hxy_y_step);
                hycal_xy_events++;
            }
        }
        // GEM↔HyCal matching residuals (ROOT recon path uses recon.gem_hits,
        // which carry detector-local x,y just like the live gem_sys hits).
        if (gem_enabled && (ep_cand || !gem_match_require_ep) && !recon.clusters.empty()) {
            const float ref_x = recon.clusters[0].x, ref_y = recon.clusters[0].y;
            const int n_dets = (int)gem_dx_hist.size();
            for (auto &gh : recon.gem_hits) {
                if (gh.det_id < 0 || gh.det_id >= n_dets) continue;
                if (gh.det_id >= (int)gem_transforms.size()) continue;
                auto &xform = gem_transforms[gh.det_id];
                float lx, ly, lz;
                xform.toLab(gh.x, gh.y, lx, ly, lz);
                float px, py;
                projectToHyCalLocal(lx, ly, lz, px, py);
                float dxr = px - ref_x, dyr = py - ref_y;
                if (std::sqrt(dxr*dxr + dyr*dyr) < gem_match_window_mm) {
                    gem_dx_hist[gh.det_id].fill(dxr, gem_resid_min, gem_resid_step);
                    gem_dy_hist[gh.det_id].fill(dyr, gem_resid_min, gem_resid_step);
                    gem_match_hits[gh.det_id]++;
                }
            }
            gem_match_events++;
        }
    }
}

json AppState::encodeReconClustersJson(const ReconEventData &recon, int ev_id)
{
    json hits_j = json::object();
    json cl_arr = json::array();

    for (size_t i = 0; i < recon.clusters.size(); ++i) {
        auto &cl = recon.clusters[i];
        std::string center_name;
        if (cl.center_id >= 0 && cl.center_id < hycal.module_count())
            center_name = hycal.module(cl.center_id).name;
        hits_j[std::to_string(cl.center_id)] =
            std::round(cl.energy * 100) / 100;
        cl_arr.push_back({
            {"id", (int)i}, {"center", center_name},
            {"center_id", cl.center_id},
            {"x", std::round(cl.x * 10) / 10},
            {"y", std::round(cl.y * 10) / 10},
            {"energy", std::round(cl.energy * 10) / 10},
            {"nblocks", cl.nblocks}, {"npos", 0},
            {"modules", json::array({cl.center_id})},
        });
    }
    return {{"event", ev_id}, {"hits", hits_j}, {"clusters", cl_arr}};
}

void AppState::prepareGemForView(const ssp::SspEventData &ssp_evt)
{
    if (!gem_enabled || ssp_evt.nmpds == 0) return;
    gem_sys.Clear();
    gem_sys.ProcessEvent(ssp_evt);
    gem_sys.Reconstruct(gem_clusterer);
}

void AppState::processGemEvent(const ssp::SspEventData &ssp_evt)
{
    if (!gem_enabled || ssp_evt.nmpds == 0) return;
    prepareGemForView(ssp_evt);

    // accumulate occupancy + histograms in a single pass
    std::lock_guard<std::mutex> lk(data_mtx);
    int total_clusters = 0;
    const int n_dets = std::min<int>(gem_sys.GetNDetectors(),
                                     (int)gem_transforms.size());
    for (int d = 0; d < n_dets; ++d) {
        auto &det = gem_sys.GetDetectors()[d];
        float xSize = det.planes[0].size;
        float ySize = det.planes[1].size;
        float xStep = xSize / GEM_OCC_NX;
        float yStep = ySize / GEM_OCC_NY;
        auto &xform = gem_transforms[d];
        auto &hits = gem_sys.GetHits(d);
        total_clusters += static_cast<int>(hits.size());
        for (auto &h : hits) {
            // Occupancy is a strip-level diagnostic — fill with raw local
            // detector coords so bin edges line up with the readout grid.
            // Lab-frame orientation lives in the matching plot, not here.
            gem_occupancy[d].fill(h.x, h.y, -xSize/2, xStep, -ySize/2, yStep);
            // theta from the target vertex — same convention as the HyCal
            // cluster theta a few hundred lines up (subtract target_*).
            float lx, ly, lz;
            xform.toLab(h.x, h.y, lx, ly, lz);
            float dx = lx - target_x, dy = ly - target_y, dz = lz - target_z;
            float r = std::sqrt(dx*dx + dy*dy);
            // dz must be > 0 for a forward-going particle; clamp to avoid
            // atan2(r, ≤0) pushing every fill into overflow when the GEM
            // z is mis-configured.
            if (dz <= 0.f) continue;
            float theta = std::atan2(r, dz) * (180.f / 3.14159265f);
            gem_theta_hist.fill(theta, gem_theta_min, gem_theta_step);
        }
    }
    gem_nclusters_hist.fill(static_cast<float>(total_clusters),
                            static_cast<float>(gem_ncl_min),
                            static_cast<float>(gem_ncl_step));
}

//=============================================================================
// GEM API builders
//=============================================================================

nlohmann::json AppState::apiGemHits() const
{
    json result = json::object();
    result["enabled"] = gem_enabled;
    if (!gem_enabled) return result;

    result["n_detectors"] = gem_sys.GetNDetectors();
    json detectors = json::array();
    json all_hits = json::array();
    for (int d = 0; d < gem_sys.GetNDetectors(); ++d) {
        auto &det = gem_sys.GetDetectors()[d];
        json dj;
        dj["name"] = det.name;
        dj["id"] = det.id;

        // 1D clusters per plane
        for (int p = 0; p < 2; ++p) {
            std::string pname = (p == 0) ? "x_clusters" : "y_clusters";
            json clusters = json::array();
            for (auto &cl : gem_sys.GetPlaneClusters(d, p)) {
                clusters.push_back({
                    {"position", cl.position},
                    {"peak_charge", cl.peak_charge},
                    {"total_charge", cl.total_charge},
                    {"size", (int)cl.hits.size()},
                    {"max_timebin", cl.max_timebin}
                });
            }
            dj[pname] = clusters;
        }

        // 2D hits (transformed to lab frame) — build per-det and all_hits in one pass.
        // proj_x/proj_y are the lab→target line projected onto the HyCal local
        // plane; the cluster tab overlays these on the geo view.
        auto &xform = gem_transforms[d];
        json hits = json::array();
        for (auto &h : gem_sys.GetHits(d)) {
            float lx, ly, lz;
            xform.toLab(h.x, h.y, lx, ly, lz);
            float px, py;
            projectToHyCalLocal(lx, ly, lz, px, py);
            hits.push_back({
                {"x", lx}, {"y", ly},
                {"proj_x", px}, {"proj_y", py},
                {"x_charge", h.x_charge}, {"y_charge", h.y_charge},
                {"x_size", h.x_size}, {"y_size", h.y_size}
            });
            all_hits.push_back({
                {"x", lx}, {"y", ly},
                {"proj_x", px}, {"proj_y", py}, {"det", d},
                {"x_charge", h.x_charge}, {"y_charge", h.y_charge}
            });
        }
        dj["hits_2d"] = hits;
        detectors.push_back(dj);
    }
    result["detectors"] = detectors;
    result["all_hits"] = all_hits;
    return result;
}

nlohmann::json AppState::apiGemConfig() const
{
    json result = json::object();
    result["enabled"] = gem_enabled;
    if (!gem_enabled) return result;

    result["n_detectors"] = gem_sys.GetNDetectors();
    json layers = json::array();
    for (int d = 0; d < gem_sys.GetNDetectors(); ++d) {
        auto &det = gem_sys.GetDetectors()[d];
        json lj = {
            {"id", det.id},
            {"name", det.name},
            {"type", det.type},
            {"x_pitch", det.planes[0].pitch},
            {"y_pitch", det.planes[1].pitch},
            {"x_apvs", det.planes[0].n_apvs},
            {"y_apvs", det.planes[1].n_apvs},
            {"x_size", det.planes[0].size},
            {"y_size", det.planes[1].size}
        };
        auto &t = gem_transforms[d];
        lj["position"] = {t.x, t.y, t.z};
        lj["tilting"]  = {t.rx, t.ry, t.rz};
        layers.push_back(lj);
    }
    result["layers"] = layers;
    result["occ_nx"] = GEM_OCC_NX;
    result["occ_ny"] = GEM_OCC_NY;
    return result;
}

nlohmann::json AppState::apiGemOccupancy() const
{
    json result = json::object();
    result["enabled"] = gem_enabled;
    if (!gem_enabled) return result;

    std::lock_guard<std::mutex> lk(data_mtx);
    json dets = json::array();
    for (int d = 0; d < gem_sys.GetNDetectors(); ++d) {
        auto &det = gem_sys.GetDetectors()[d];
        json dj;
        dj["id"] = det.id;
        dj["name"] = det.name;
        dj["x_size"] = det.planes[0].size;
        dj["y_size"] = det.planes[1].size;
        dj["nx"] = GEM_OCC_NX;
        dj["ny"] = GEM_OCC_NY;
        dj["bins"] = gem_occupancy[d].bins;
        dets.push_back(dj);
    }
    result["detectors"] = dets;
    result["total"] = events_processed.load();
    return result;
}

nlohmann::json AppState::apiGemHist() const
{
    std::lock_guard<std::mutex> lk(data_mtx);
    return {
        {"nclusters", histToJson(gem_nclusters_hist, (float)gem_ncl_min, (float)gem_ncl_max, (float)gem_ncl_step)},
        {"theta",     histToJson(gem_theta_hist, gem_theta_min, gem_theta_max, gem_theta_step)}
    };
}

nlohmann::json AppState::apiGemApv(const ssp::SspEventData &ssp_evt, int evnum) const
{
    json result = json::object();
    result["enabled"] = gem_enabled;
    result["event"]   = evnum;
    if (!gem_enabled) {
        result["detectors"] = json::array();
        result["apvs"]      = json::array();
        return result;
    }
    // Global software N-sigma multiplier and pedestal calibration revision.
    // The frontend pairs this zs_sigma with the per-APV noise[] from
    // /api/gem/calib to draw the threshold band; gem_calib_rev lets it
    // detect when the cached calib payload is stale and needs re-fetching.
    result["zs_sigma"]  = gem_sys.GetZeroSupThreshold();
    result["calib_rev"] = gem_calib_rev.load();

    // Detector summary — det_id → name + APV count, used by the frontend
    // to lay out one section per GEM with consistent ordering.
    auto &dets = gem_sys.GetDetectors();
    json det_arr = json::array();
    std::vector<int> apv_counts(dets.size(), 0);
    for (int i = 0; i < gem_sys.GetNApvs(); ++i) {
        auto &cfg = gem_sys.GetApvConfig(i);
        if (cfg.det_id >= 0 && cfg.det_id < (int)apv_counts.size())
            apv_counts[cfg.det_id]++;
    }
    for (size_t d = 0; d < dets.size(); ++d) {
        det_arr.push_back({
            {"id",      dets[d].id},
            {"name",    dets[d].name},
            {"n_apvs",  apv_counts[d]},
        });
    }
    result["detectors"] = det_arr;

    // Per-APV dump.  Each APV carries:
    //   raw[128][6]        — int16 firmware samples (0 if APV not in event)
    //   processed[128][6]  — pedestal + CM corrected float (0 if not present)
    //   hits[128]          — software ZS survivor (post-cut), 0/1
    //   fw_hits[128]       — firmware survivor (pre-software-cut), 0/1
    //   cm[6] | null       — firmware online common mode per time sample
    //   no_hit_fr          — firmware full-readout (nstrips==128) but no survivors
    //   present            — APV showed up in this event's SSP data
    // Per-APV pedestal noise lives on /api/gem/calib (one-shot, cached
    // by the frontend until calib_rev changes).
    constexpr int N_STRIPS = 128;
    constexpr int N_TS     = 6;
    json apvs = json::array();
    for (int i = 0; i < gem_sys.GetNApvs(); ++i) {
        auto &cfg = gem_sys.GetApvConfig(i);
        if (cfg.crate_id < 0 || cfg.mpd_id < 0 || cfg.adc_ch < 0)
            continue;   // unmapped slot

        const ssp::ApvData *raw = ssp_evt.findApv(cfg.crate_id, cfg.mpd_id, cfg.adc_ch);
        bool present = (raw != nullptr) && raw->present;

        json raw_arr = json::array();
        json proc_arr = json::array();
        json hit_arr = json::array();
        json fw_hit_arr = json::array();
        bool any_hit = false;

        for (int s = 0; s < N_STRIPS; ++s) {
            json raw_row  = json::array();
            json proc_row = json::array();
            for (int t = 0; t < N_TS; ++t) {
                if (present)
                    raw_row.push_back(static_cast<int>(raw->strips[s][t]));
                else
                    raw_row.push_back(0);
                // 1-decimal rounding keeps payload tight without losing
                // anything visible at the panel's sub-pixel resolution.
                float v = present ? gem_sys.GetProcessedAdc(i, s, t) : 0.f;
                proc_row.push_back(std::round(v * 10.f) / 10.f);
            }
            raw_arr.push_back(std::move(raw_row));
            proc_arr.push_back(std::move(proc_row));
            bool hit = present && gem_sys.IsChannelHit(i, s);
            if (hit) any_hit = true;
            hit_arr.push_back(hit ? 1 : 0);
            fw_hit_arr.push_back(present && raw->hasStrip(s) ? 1 : 0);
        }

        // Firmware online CM (6 samples) — only present when the MPD emitted
        // type-0xD debug-header words; otherwise null so the frontend can
        // skip the overlay rather than draw zeros.
        json cm_val = nullptr;
        if (present && raw->has_online_cm) {
            json cm_arr = json::array();
            for (int t = 0; t < N_TS; ++t)
                cm_arr.push_back(static_cast<int>(raw->online_cm[t]));
            cm_val = std::move(cm_arr);
        }

        std::string det_name;
        if (cfg.det_id >= 0 && cfg.det_id < (int)dets.size())
            det_name = dets[cfg.det_id].name;
        const char *plane = (cfg.plane_type == 0) ? "X"
                          : (cfg.plane_type == 1) ? "Y" : "?";

        // Firmware full-readout warning: APV reported every strip but none
        // survived ZS.  Mirrors the gem_event_viewer.py "no hits" badge.
        bool no_hit_fr = present && raw->nstrips >= N_STRIPS && !any_hit;

        apvs.push_back({
            {"id",         i},
            {"det_id",     cfg.det_id},
            {"det_name",   det_name},
            {"plane",      plane},
            {"det_pos",    cfg.det_pos},
            {"crate",      cfg.crate_id},
            {"mpd",        cfg.mpd_id},
            {"adc",        cfg.adc_ch},
            {"present",    present},
            {"no_hit_fr",  no_hit_fr},
            {"raw",        std::move(raw_arr)},
            {"processed",  std::move(proc_arr)},
            {"hits",       std::move(hit_arr)},
            {"fw_hits",    std::move(fw_hit_arr)},
            {"cm",         std::move(cm_val)},
        });
    }
    result["apvs"] = apvs;
    return result;
}

nlohmann::json AppState::apiGemCalib() const
{
    json result = json::object();
    result["enabled"]   = gem_enabled;
    result["rev"]       = gem_calib_rev.load();
    result["zs_sigma"]  = gem_enabled ? gem_sys.GetZeroSupThreshold() : 0.f;
    json apvs = json::array();
    if (gem_enabled) {
        constexpr int N_STRIPS = 128;
        for (int i = 0; i < gem_sys.GetNApvs(); ++i) {
            auto &cfg = gem_sys.GetApvConfig(i);
            if (cfg.crate_id < 0 || cfg.mpd_id < 0 || cfg.adc_ch < 0)
                continue;
            json noise_arr = json::array();
            for (int s = 0; s < N_STRIPS; ++s)
                noise_arr.push_back(std::round(cfg.pedestal[s].noise * 10.f) / 10.f);
            apvs.push_back({{"id", i}, {"noise", std::move(noise_arr)}});
        }
    }
    result["apvs"] = std::move(apvs);
    return result;
}

void AppState::setGemZsSigma(float v)
{
    if (v < 0.f) v = 0.f;
    gem_sys.SetZeroSupThreshold(v);
}

//=============================================================================
// Clearing
//=============================================================================

void AppState::clearHistograms()
{
    std::lock_guard<std::mutex> lk(data_mtx);
    for (auto &[k, h] : histograms)        h.clear();
    for (auto &[k, h] : pos_histograms)   h.clear();
    for (auto &[k, h] : height_histograms) h.clear();
    occupancy.clear();
    occupancy_tcut.clear();
    events_processed = 0;
    cluster_energy_hist.clear();
    nclusters_hist.clear();
    nblocks_hist.clear();
    for (auto &h : cluster_energy_hist_by_ncl) h.clear();
    for (auto &h : nblocks_hist_by_ncl)        h.clear();
    energy_angle_hist.clear();
    moller_xy_hist.clear();
    moller_events = 0;
    hycal_xy_hist.clear();
    hycal_xy_events = 0;
    for (auto &h : gem_dx_hist) h.clear();
    for (auto &h : gem_dy_hist) h.clear();
    for (auto &n : gem_match_hits) n = 0;
    gem_match_events = 0;
    cluster_events_processed = 0;
    for (auto &h : gem_occupancy) h.clear();
    gem_nclusters_hist.clear();
    gem_theta_hist.clear();
}

void AppState::clearLms()
{
    std::lock_guard<std::mutex> lk(lms_mtx);
    lms_history.clear();
    latest_lms_integral.clear();
    latest_alpha_integral.clear();
    lms_events = 0;
    lms_first_ts = 0;
    sync_unix = 0;
    sync_rel_sec = 0.;
    pending_sync_unix = 0;
    pending_sync_ti = 0;
}

