//=============================================================================
// Replay.cpp — EVIO to ROOT tree conversion
//=============================================================================

#include "Replay.h"
#include "DaqConfig.h"
#include "EventData_io.h"
#include "PulseTemplateStore.h"
#include "HyCalSystem.h"
#include "GemSystem.h"
#include "HyCalCluster.h"
#include "GemCluster.h"
#include "MatchingTools.h"
#include "ConfigSetup.h"
#include "InstallPaths.h"
#include "gain_factor.h"

#include <nlohmann/json.hpp>
#include <fstream>
#include <iostream>

using json = nlohmann::json;

#ifndef DATABASE_DIR
#define DATABASE_DIR "."
#endif

namespace analysis {

void Replay::LoadDaqMap(const std::string &json_path)
{
    std::ifstream f(json_path);
    if (!f.is_open()) {
        std::cerr << "Replay: cannot open DAQ map: " << json_path << "\n";
        return;
    }
    auto j = json::parse(f, nullptr, false, true);
    if (j.is_array()) {
        for (auto &entry : j) {
            std::string name = entry.value("name", "");
            int crate   = entry.value("crate", -1);
            int slot    = entry.value("slot", -1);
            int channel = entry.value("channel", -1);
            if (!name.empty() && crate >= 0)
                daq_map_[std::to_string(crate) + "_" + std::to_string(slot) +
                         "_" + std::to_string(channel)] = name;
        }
    }
    std::cerr << "Replay: loaded " << daq_map_.size() << " DAQ map entries\n";
}

void Replay::LoadModulesInfo(const std::string &json_path)
{
    std::ifstream f(json_path);
    if (!f.is_open()) {
        std::cerr << "Replay: cannot open modules info: " << json_path << "\n";
        return;
    }
    auto j = json::parse(f, nullptr, false, true);
    if (!j.is_array()) return;

    auto parse_t = [](const std::string &t) {
        if (t == "PbGlass") return prad2::MOD_PbGlass;
        if (t == "PbWO4")   return prad2::MOD_PbWO4;
        if (t == "SCINT")   return prad2::MOD_SCINT;
        if (t == "LMS")     return prad2::MOD_LMS;
        return prad2::MOD_UNKNOWN;
    };

    for (auto &entry : j) {
        std::string name = entry.value("n", "");
        std::string t    = entry.value("t", "");
        if (name.empty()) continue;
        module_types_[name] = parse_t(t);
    }
    std::cerr << "Replay: loaded " << module_types_.size()
              << " module-type entries from " << json_path << "\n";
}

std::string Replay::moduleName(int roc, int slot, int ch) const
{
    auto it = daq_map_.find(std::to_string(roc) + "_" + std::to_string(slot) +
                            "_" + std::to_string(ch));
    return (it != daq_map_.end()) ? it->second : "";
}

prad2::ModuleType Replay::moduleType(int roc, int slot, int ch) const
{
    auto name = moduleName(roc, slot, ch);
    if (name.empty()) return prad2::MOD_UNKNOWN;
    auto it = module_types_.find(name);
    return (it != module_types_.end()) ? it->second : prad2::MOD_UNKNOWN;
}

int Replay::moduleID(int roc, int slot, int ch) const
{
    auto name = moduleName(roc, slot, ch);
    if (name.empty()) return -1;
    auto t = moduleType(roc, slot, ch);
    // Globally-unique ID encoding — see RawEventData docs.  The numeric
    // ranges are deliberately disjoint so HyCalSystem::module_by_id(...)
    // returns nullptr for SCINT/LMS, letting existing HyCal consumers
    // skip them via their existing nullptr / is_hycal() checks.
    switch (t) {
        case prad2::MOD_PbGlass:
            // "G<n>" → n.  std::stoi tolerates trailing junk; in practice the
            // map contains pure "G123" entries.
            try { return std::stoi(name.substr(1)); } catch (...) { return -1; }

        case prad2::MOD_PbWO4:
            try { return std::stoi(name.substr(1)) + 1000; } catch (...) { return -1; }

        case prad2::MOD_SCINT:
            // "V1".."V4" → 3001..3004
            if (name.size() >= 2 && name[0] == 'V')
                try { return 3000 + std::stoi(name.substr(1)); } catch (...) {}
            return -1;

        case prad2::MOD_LMS:
            // "LMSPin"=3100, "LMS1".."LMS3"=3101..3103
            if (name == "LMSPin") return 3100;
            if (name.rfind("LMS", 0) == 0 && name.size() == 4) {
                char d = name[3];
                if (d >= '1' && d <= '9') return 3100 + (d - '0');
            }
            return -1;

        case prad2::MOD_UNKNOWN:
        default:
            return -1;
    }
}

void Replay::clearEvent(EventVars &ev)
{
    ev.event_num = 0;
    ev.trigger_type = 0;
    ev.trigger_bits = 0;
    ev.timestamp = 0;
    ev.nch = 0;
    ev.gem_nch = 0;
    ev.ssp_raw.clear();
    std::fill(std::begin(ev.npeaks), std::end(ev.npeaks), 0);
    std::fill(&ev.peak_height[0][0],   &ev.peak_height[0][0]   + prad2::kMaxChannels * fdec::MAX_PEAKS, 0.f);
    std::fill(&ev.peak_time[0][0],     &ev.peak_time[0][0]     + prad2::kMaxChannels * fdec::MAX_PEAKS, 0.f);
    std::fill(&ev.peak_integral[0][0], &ev.peak_integral[0][0] + prad2::kMaxChannels * fdec::MAX_PEAKS, 0.f);
}

void Replay::clearReconEvent(EventVars_Recon &ev)
{
    ev.event_num = 0;
    ev.trigger_type = 0;
    ev.trigger_bits = 0;
    ev.timestamp = 0;
    ev.total_energy = 0.f;
    ev.n_clusters = 0;
    ev.n_gem_hits = 0;
    ev.matchNum = 0;
    std::fill(std::begin(ev.matchFlag), std::end(ev.matchFlag), 0);
    ev.veto_nch = 0;
    ev.lms_nch = 0;
    ev.ssp_raw.clear();
    std::fill(std::begin(ev.veto_npeaks), std::end(ev.veto_npeaks), 0);
    std::fill(&ev.veto_peak_time[0][0],     &ev.veto_peak_time[0][0]     + 4 * fdec::MAX_PEAKS, 0.f);
    std::fill(&ev.veto_peak_height[0][0],   &ev.veto_peak_height[0][0]   + 4 * fdec::MAX_PEAKS, 0.f);
    std::fill(&ev.veto_peak_integral[0][0], &ev.veto_peak_integral[0][0] + 4 * fdec::MAX_PEAKS, 0.f);
    std::fill(std::begin(ev.lms_npeaks), std::end(ev.lms_npeaks), 0);
    std::fill(&ev.lms_peak_time[0][0],     &ev.lms_peak_time[0][0]     + 4 * fdec::MAX_PEAKS, 0.f);
    std::fill(&ev.lms_peak_height[0][0],   &ev.lms_peak_height[0][0]   + 4 * fdec::MAX_PEAKS, 0.f);
    std::fill(&ev.lms_peak_integral[0][0], &ev.lms_peak_integral[0][0] + 4 * fdec::MAX_PEAKS, 0.f);
}

void Replay::setupBranches(TTree *tree, EventVars &ev, bool write_peaks)
{
    prad2::SetRawWriteBranches(tree, ev, write_peaks);
}

void Replay::setupReconBranches(TTree *tree, EventVars_Recon &ev)
{
    prad2::SetReconWriteBranches(tree, ev);
}

bool Replay::Process(const std::string &input_evio, const std::string &output_root, RunConfig &gRunConfig,
                     const std::string &db_dir,
                     int max_events, bool write_peaks , const std::string &daq_config_file)
{
    // build ROC tag → crate index mapping from DAQ config JSON
    std::unordered_map<int, int> roc_to_crate;
    if (!daq_config_file.empty()) {
        std::cout << "Loading DAQ config from " << daq_config_file << "\n";
        std::ifstream dcf(daq_config_file);
        if (dcf.is_open()) {
            auto dcj = nlohmann::json::parse(dcf, nullptr, false, true);
            if (dcj.contains("roc_tags") && dcj["roc_tags"].is_array()) {
                for (auto &entry : dcj["roc_tags"]) {
                    int tag   = std::stoi(entry.at("tag").get<std::string>(), nullptr, 16);
                    int crate = entry.at("crate").get<int>();
                    roc_to_crate[tag] = crate;
                }
            }
        }
    }
    else {
        std::cerr << "No DAQ config file provided, ROC tag to crate mapping will be unavailable.\n";
    }

    evc::EvChannel ch;
    ch.SetConfig(daq_cfg_);

    if (ch.OpenAuto(input_evio) != evc::status::success) {
        std::cerr << "Replay: cannot open " << input_evio << "\n";
        return false;
    }

    TFile *outfile = TFile::Open(output_root.c_str(), "RECREATE");
    if (!outfile || !outfile->IsOpen()) {
        std::cerr << "Replay: cannot create " << output_root << "\n";
        return false;
    }

    TTree *tree = new TTree("events", "PRad2 replay data");
    //EventVars ev;
    auto ev = std::make_unique<EventVars>();
    setupBranches(tree, *ev, write_peaks);

    auto event = std::make_unique<fdec::EventData>();
    auto ssp_evt = std::make_unique<ssp::SspEventData>();
    fdec::WaveAnalyzer ana(daq_cfg_.wave_cfg);
    // NNLS pile-up template store (loaded only if config asks for it;
    // failure is non-fatal — the analyzer falls back to local-maxima
    // peak heights).
    fdec::PulseTemplateStore template_store;
    if (daq_cfg_.wave_cfg.nnls_deconv.enabled
        && !daq_cfg_.wave_cfg.nnls_deconv.template_file.empty()) {
        template_store.LoadFromFile(
            db_dir + "/" + daq_cfg_.wave_cfg.nnls_deconv.template_file,
            daq_cfg_.wave_cfg);
    }
    ana.SetTemplateStore(&template_store);
    fdec::WaveResult wres;
    // Firmware-mode emulator (FADC250 Modes 1/2/3).  Configured from the
    // optional "fadc250_waveform.firmware" block in daq_config.json — defaults
    // are safe for DAQ signal studies but should be overridden to match the
    // actual run's TET/NSB/NSA/MAX_PULSES if comparing to firmware output.
    fdec::Fadc250FwAnalyzer fw_ana(daq_cfg_.fadc250_fw);
    fdec::DaqWaveResult dwres;
    int total = 0;

    int run_num = get_run_int(input_evio);
    auto gain_correction = prad2::ComputeGainCorrection(db_dir + "/" + gRunConfig.gain_data_dir, run_num, gRunConfig.gain_ref_run);

    while (ch.Read() == evc::status::success) {
        if (!ch.Scan()) continue;
        if (ch.GetEventType() != evc::EventType::Physics) continue;

        // Snapshot raw 0xE10C SSP trigger bank for this read group (one bank
        // per CODA event, shared by all sub-events from this Read()).
        std::vector<uint32_t> ssp_raw_snapshot;
        if (auto *n_e10c = ch.FindFirstByTag(0xE10C)) {
            const uint32_t *p = ch.GetData(*n_e10c);
            ssp_raw_snapshot.assign(p, p + n_e10c->data_words);
        }

        for (int ie = 0; ie < ch.GetNEvents(); ++ie) {
            event->clear();
            ssp_evt->clear();
            if (!ch.DecodeEvent(ie, *event, ssp_evt.get())) continue;
            if (max_events > 0 && total >= max_events) break;

            clearEvent(*ev);
            ev->event_num    = event->info.event_number;
            ev->trigger_type = event->info.trigger_type;
            ev->trigger_bits      = event->info.trigger_bits;
            ev->timestamp    = event->info.timestamp;
            ev->ssp_raw      = ssp_raw_snapshot;

            // Decode FADC250 data — single pass over all channels (HyCal +
            // Veto + LMS).  Type dispatch comes from hycal_modules.json,
            // not module-name prefix; module_type[nch] records the category.
            int nch = 0;
            for (int r = 0; r < event->nrocs; ++r) {
                auto &roc = event->rocs[r];
                if (!roc.present) continue;
                auto cit = roc_to_crate.find(roc.tag);
                int crate = (cit == roc_to_crate.end()) ? (int)roc.tag : cit->second;
                for (int s = 0; s < fdec::MAX_SLOTS; ++s) {
                    if (!roc.slots[s].present) continue;
                    for (int c = 0; c < 16; ++c) {
                        if (!(roc.slots[s].channel_mask & (1ull << c))) continue;
                        auto &cd = roc.slots[s].channels[c];
                        if (cd.nsamples <= 0 || nch >= prad2::kMaxChannels) continue;

                        int  mod_id   = moduleID(crate, s, c);
                        auto mod_type = moduleType(crate, s, c);
                        // Drop channels with no DAQ-map / module-info entry —
                        // we have no way to interpret them downstream.
                        if (mod_id < 0) continue;

                        ev->module_id[nch]   = static_cast<uint16_t>(mod_id);
                        ev->module_type[nch] = static_cast<uint8_t>(mod_type);
                        ev->nsamples[nch]    = static_cast<uint8_t>(cd.nsamples);
                        for (int i = 0; i < cd.nsamples && i < fdec::MAX_SAMPLES; ++i)
                            ev->samples[nch][i] = cd.samples[i];

                        // Gain correction is HyCal-only (PbGlass / PbWO4) —
                        // SCINT / LMS get unity factor.  Comes from a lookup
                        // table, no analyzer needed.
                        if (mod_type == prad2::MOD_PbWO4) {
                            ev->gain_factor[nch] = gain_correction.w[mod_id - 1000].avg;
                        } else if (mod_type == prad2::MOD_PbGlass) {
                            ev->gain_factor[nch] = gain_correction.g[mod_id].avg;
                        } else {
                            ev->gain_factor[nch] = 1.0f;
                        }

                        if (write_peaks) {
                            // Soft analyzer drives both peaks AND the
                            // pedestal estimate that the firmware analyzer
                            // consumes — only run it when its output is
                            // being written.
                            ana.SetChannelKey(roc.tag, s, c);
                            ana.Analyze(cd.samples, cd.nsamples, wres);
                            ev->ped_mean[nch]    = wres.ped.mean;
                            ev->ped_rms[nch]     = wres.ped.rms;
                            ev->ped_nused[nch]   = wres.ped.nused;
                            ev->ped_quality[nch] = wres.ped.quality;
                            ev->ped_slope[nch]   = wres.ped.slope;
                            ev->npeaks[nch]   = static_cast<uint8_t>(wres.npeaks);
                            for (int p = 0; p < wres.npeaks && p < fdec::MAX_PEAKS; p++) {
                                ev->peak_height[nch][p]   = wres.peaks[p].height;
                                ev->peak_time[nch][p]     = wres.peaks[p].time;
                                ev->peak_integral[nch][p] = wres.peaks[p].integral;
                                ev->peak_quality[nch][p]  = wres.peaks[p].quality;
                            }
                            fw_ana.Analyze(cd.samples, cd.nsamples, wres.ped.mean, dwres);
                            ev->daq_npeaks[nch] = static_cast<uint8_t>(dwres.npeaks);
                            for (int p = 0; p < dwres.npeaks && p < fdec::MAX_PEAKS; ++p) {
                                const auto &dp = dwres.peaks[p];
                                ev->daq_peak_vp[nch][p]       = dp.vpeak;
                                ev->daq_peak_integral[nch][p] = dp.integral;
                                ev->daq_peak_time[nch][p]     = dp.time_ns;
                                ev->daq_peak_cross[nch][p]    = dp.cross_sample;
                                ev->daq_peak_pos[nch][p]      = dp.peak_sample;
                                ev->daq_peak_coarse[nch][p]   = dp.coarse;
                                ev->daq_peak_fine[nch][p]     = dp.fine;
                                ev->daq_peak_quality[nch][p]  = dp.quality;
                            }
                        }
                        nch++;
                    }
                }
            }
            ev->nch = nch;

            // decode GEM SSP data
            int gem_ch = 0;
            for (int m = 0; m < ssp_evt->nmpds; ++m) {
                auto &mpd = ssp_evt->mpds[m];
                if (!mpd.present) continue;
                for (int a = 0; a < ssp::MAX_APVS_PER_MPD; ++a) {
                    auto &apv = mpd.apvs[a];
                    if (!apv.present) continue;
                    int idx = -1; // find APV index in GemSystem if needed
                    for (int s = 0; s < ssp::APV_STRIP_SIZE; ++s) {
                        if (!apv.hasStrip(s)) continue;
                        if (gem_ch >= prad2::kMaxGemStrips) continue;
                        
                        ev->mpd_crate[gem_ch] = mpd.crate_id;
                        ev->mpd_fiber[gem_ch] = mpd.mpd_id;
                        ev->apv[gem_ch]       = a;
                        ev->strip[gem_ch]     = s;
                        for (int t = 0; t < ssp::SSP_TIME_SAMPLES; t++)
                            ev->ssp_samples[gem_ch][t] = apv.strips[s][t];

                        gem_ch++;
                    }
                }
            }
            ev->gem_nch = gem_ch; // total channels = HyCal + GEM
            tree->Fill();
            total++;

            if (total % 1000 == 0)
                std::cerr << "\rReplay: " << total << " events processed" << std::flush;
        }
        if (max_events > 0 && total >= max_events) break;
    }

    std::cerr << "\rReplay: " << total << " events written to " << output_root << "\n";
    tree->Write();
    delete outfile;
    return true;
}

bool Replay::ProcessWithRecon(const std::string &input_evio, const std::string &output_root, RunConfig &gRunConfig,
                                const std::string &db_dir,
                                const std::string &daq_config_file, const std::string &gem_ped_file,
                                const float zerosup_override, bool prad1)
{
    // Similar to Process(), but with HyCal reconstruction and GEM hit reconstruction
    // before filling the ROOT tree.
    // The main differences are:
    // - After decoding, we run the HyCal clusterer to reconstruct clusters and hits.
    // - We also run the GemSystem reconstruction to get GEM hits.
    // - We fill a different TTree with reconstructed quantities instead of raw data.

    // build ROC tag → crate index mapping from DAQ config JSON
    std::unordered_map<int, int> roc_to_crate;
    if (!daq_config_file.empty()) {
        std::cout << "Loading DAQ config from " << daq_config_file << "\n";
        std::ifstream dcf(daq_config_file);
        if (dcf.is_open()) {
            auto dcj = nlohmann::json::parse(dcf, nullptr, false, true);
            if (dcj.contains("roc_tags") && dcj["roc_tags"].is_array()) {
                for (auto &entry : dcj["roc_tags"]) {
                    int tag   = std::stoi(entry.at("tag").get<std::string>(), nullptr, 16);
                    int crate = entry.at("crate").get<int>();
                    roc_to_crate[tag] = crate;
                }
            }
        }
    }
    else {
        std::cerr << "No DAQ config file provided, ROC tag to crate mapping will be unavailable.\n";
    }

    // Setup HyCal system and clusterer
    fdec::HyCalSystem hycal;
     std::string daq_map_file = db_dir + "/hycal_daq_map.json";
    if(prad1 == true)
        daq_map_file = db_dir + "/prad1/prad_hycal_daq_map.json";
    hycal.Init(db_dir + "/hycal_modules.json", daq_map_file);
    
    if(prad1 == true) evc::load_pedestals(db_dir + "/prad1/adc1881m_pedestals.json", daq_cfg_);

    std::string calib_file = db_dir + "/" + gRunConfig.energy_calib_file;
    int nmatched = hycal.LoadCalibration(calib_file);
    if (nmatched >= 0)
        std::cerr << "Calibration: " << calib_file << " (" << nmatched << " modules)\n";

    fdec::HyCalCluster clusterer(hycal);
    fdec::ClusterConfig cl_cfg;
    clusterer.SetConfig(cl_cfg);

    MatchingTools matching;

    // Initialize GEM system and clusterer
    std::unique_ptr<gem::GemSystem> gem_sys;
    std::unique_ptr<gem::GemCluster> gem_clusterer;
if(!prad1){
    gem_sys = std::make_unique<gem::GemSystem>();
    std::string gem_map_file = db_dir + "/gem_daq_map.json";
    gem_sys->Init(gem_map_file);
    std::cerr << "GEM map  : " << gem_map_file
                << " (" << gem_sys->GetNDetectors() << " detectors)\n";

    if (!gem_ped_file.empty()) {
        gem_sys->LoadPedestals(gem_ped_file);
        std::cerr << "GEM peds : " << gem_ped_file << "\n";
    }
    else {
        gem_sys->LoadPedestals(db_dir + "/" + gRunConfig.gem_pedestal_file);
            std::cerr << "GEM peds : " << db_dir + "/" + gRunConfig.gem_pedestal_file << "\n";

    }

    if (zerosup_override >= 0.f) {
        gem_sys->SetZeroSupThreshold(zerosup_override);
        std::cerr << "Zero-sup : " << zerosup_override << " sigma (override)\n";
    }
    
    gem_clusterer = std::make_unique<gem::GemCluster>();
}
    //open EVIO file and output ROOT file
    evc::EvChannel ch;
    ch.SetConfig(daq_cfg_);

    if (ch.OpenAuto(input_evio) != evc::status::success) {
        std::cerr << "Replay: cannot open " << input_evio << "\n";
        return false;
    }

    TFile *outfile = TFile::Open(output_root.c_str(), "RECREATE");
    if (!outfile || !outfile->IsOpen()) {
        std::cerr << "Replay: cannot create " << output_root << "\n";
        return false;
    }

    // create TTree and branches for reconstructed data
    TTree *tree = new TTree("recon", "PRad2 replay reconstruction");
    auto ev = std::make_unique<EventVars_Recon>();
    setupReconBranches(tree, *ev);

    //initialize tools for event decoder and cluster reconstruction
    auto event = std::make_unique<fdec::EventData>();
    auto ssp_evt = std::make_unique<ssp::SspEventData>();
    fdec::WaveAnalyzer ana(daq_cfg_.wave_cfg);
    fdec::PulseTemplateStore template_store;
    if (daq_cfg_.wave_cfg.nnls_deconv.enabled
        && !daq_cfg_.wave_cfg.nnls_deconv.template_file.empty()) {
        template_store.LoadFromFile(
            db_dir + "/" + daq_cfg_.wave_cfg.nnls_deconv.template_file,
            daq_cfg_.wave_cfg);
    }
    ana.SetTemplateStore(&template_store);
    fdec::WaveResult wres;

    int total = 0;

    int run_num = get_run_int(input_evio);
    auto gain_correction = prad2::ComputeGainCorrection(db_dir + "/" + gRunConfig.gain_data_dir, run_num, gRunConfig.gain_ref_run);

    while (ch.Read() == evc::status::success) {
        if (!ch.Scan()) continue;
        if (ch.GetEventType() != evc::EventType::Physics) continue;

        // Snapshot raw 0xE10C SSP trigger bank for this read group.
        std::vector<uint32_t> ssp_raw_snapshot;
        if (auto *n_e10c = ch.FindFirstByTag(0xE10C)) {
            const uint32_t *p = ch.GetData(*n_e10c);
            ssp_raw_snapshot.assign(p, p + n_e10c->data_words);
        }

        for (int ie = 0; ie < ch.GetNEvents(); ++ie) {
            event->clear();
            ssp_evt->clear();
            clusterer.Clear();
            if (!ch.DecodeEvent(ie, *event, ssp_evt.get())) continue;

            clearReconEvent(*ev);
            ev->event_num    = event->info.event_number;
            ev->trigger_type = event->info.trigger_type;
            ev->trigger_bits = event->info.trigger_bits;
            ev->timestamp    = event->info.timestamp;
            ev->ssp_raw      = ssp_raw_snapshot;

            // TODO: use config-driven trigger filter (monitor_config.json "physics" section
            // accept_trigger_bits/reject_trigger_bits) instead of hardcoded bit check.
            // Currently drops all non-SSP_RawSum events, including LMS.
            static constexpr uint32_t TBIT_sum = (1u << 8);
            static constexpr uint32_t TBIT_lms = (1u << 24);
            static constexpr uint32_t TBIT_alpha = (1u << 25);
            bool is_sum = (ev->trigger_bits & TBIT_sum) != 0;
            bool is_lms = (ev->trigger_bits & TBIT_lms) != 0;
            bool is_alpha = (ev->trigger_bits & TBIT_alpha) != 0;
            if (!is_sum && !is_lms && !is_alpha) continue;

            // decode FADC250 and reconstruct HyCal data
            int veto_nch = 0;
            int lms_nch = 0;
            int nch = 0;
            for (int r = 0; r < event->nrocs; ++r) {
                auto &roc = event->rocs[r];
                if (!roc.present) continue;
                auto cit = roc_to_crate.find(roc.tag);
                if (cit == roc_to_crate.end()) continue;
                int crate = cit->second;
                for (int s = 0; s < fdec::MAX_SLOTS; ++s) {
                    if (!roc.slots[s].present) continue;
                    for (int c = 0; c < 64; ++c) { //should be 16, a bigger number to adapt PRad1 data
                        if (!(roc.slots[s].channel_mask & (1ull << c))) continue;
                        auto &cd = roc.slots[s].channels[c];
                        if (cd.nsamples <= 0) continue;

                        std::string mod_name = moduleName(crate, s, c);
                        if(mod_name.empty()) continue;

                        if(is_lms || is_alpha) {
                            if(mod_name[0] == 'L'){
                                if(mod_name.length() != 4) continue;
                                if(lms_nch >= 4) { lms_nch++; continue; } // guard against overflow
                                if(mod_name[3] == 'P') ev->lms_id[lms_nch] = 0;
                                else ev->lms_id[lms_nch] = mod_name[3] - '0';
                                ana.SetChannelKey(roc.tag, s, c);
                                ana.Analyze(cd.samples, cd.nsamples, wres);
                                ev->lms_npeaks[lms_nch] = wres.npeaks;
                                if(wres.npeaks <= 0) continue;
                                for (int p = 0; p < wres.npeaks && p < fdec::MAX_PEAKS; ++p) {
                                    ev->lms_peak_height[lms_nch][p] = wres.peaks[p].height;
                                    ev->lms_peak_integral[lms_nch][p] = wres.peaks[p].integral;
                                    ev->lms_peak_time[lms_nch][p] = wres.peaks[p].time;
                                }
                                lms_nch++;
                            }
                            else continue;
                        }

                        if(is_sum && !is_lms) {
                            if(mod_name[0] == 'V'){
                                if(mod_name.length() != 2) continue;
                                if(veto_nch >= 4) { veto_nch++; continue; } // guard against overflow
                                ev->veto_id[veto_nch] = mod_name[1] - '0';
                                ana.SetChannelKey(roc.tag, s, c);
                                ana.Analyze(cd.samples, cd.nsamples, wres);
                                ev->veto_npeaks[veto_nch] = wres.npeaks;
                                if(wres.npeaks <= 0) continue;
                                for (int p = 0; p < wres.npeaks && p < fdec::MAX_PEAKS; ++p) {
                                    ev->veto_peak_height[veto_nch][p] = wres.peaks[p].height;
                                    ev->veto_peak_integral[veto_nch][p] = wres.peaks[p].integral;
                                    ev->veto_peak_time[veto_nch][p] = wres.peaks[p].time;
                                }
                                veto_nch++;
                            }
                            else{
                                const auto *mod = hycal.module_by_daq(crate, s, c);
                                if (!mod || !mod->is_hycal()) continue;
                                float adc = 0.f, time = 0.f;
                                if(prad1 == true) 
                                    adc = cd.samples[0] * 0.543; //0.543 for prad1 run1308,correct to 1.1GeV
                                else{
                                    ana.SetChannelKey(roc.tag, s, c);
                                    ana.Analyze(cd.samples, cd.nsamples, wres);
                                    if (wres.npeaks <= 0) continue;
                                    int bestIdx = -1;
                                    float bestHeight = -1.f;
                                    for(int p = 0; p < wres.npeaks && p < fdec::MAX_PEAKS; ++p){
                                        if(wres.peaks[p].time > gRunConfig.hc_time_win_lo &&
                                        wres.peaks[p].time < gRunConfig.hc_time_win_hi) {
                                            if(wres.peaks[p].height > bestHeight) {
                                                bestHeight = wres.peaks[p].height;
                                                bestIdx = p;
                                            }
                                        }
                                    }
                                    if (bestIdx < 0) continue;
                                    adc = wres.peaks[bestIdx].integral;
                                    time = wres.peaks[bestIdx].time;
                                }
                                //gain correction for HyCal modules
                                if(mod->id > 1000) adc *= gain_correction.w[mod->id-1000].avg;
                                else adc *= gain_correction.g[mod->id].avg;

                                float energy = static_cast<float>(mod->energize(adc));
                                clusterer.AddHit(mod->index, energy, time);
                                ev->total_energy += energy;
                                nch++;
                            }
                        }
                    }
                }
            }
            ev->veto_nch = veto_nch;
            ev->lms_nch = lms_nch;
            if(nch > 500) continue; // too many hits, likely noise, skip the event

            clusterer.FormClusters();
            std::vector<fdec::ClusterHit> hits;
            clusterer.ReconstructHits(hits);
            //HyCal event reconstrued, fill root tree and histograms
            ev->n_clusters = std::min((int)hits.size(), prad2::kMaxClusters);
            for (int i = 0; i < ev->n_clusters; ++i) {
                ev->cl_nblocks[i] = hits[i].nblocks;
                ev->cl_time[i]    = hits[i].time;
                //transform the cluster positions to the lab coordinate
                HCHit local_hit = {hits[i].x, hits[i].y, fdec::shower_depth(hits[i].center_id, hits[i].energy), 
                    hits[i].energy, static_cast<uint16_t>(hits[i].center_id), hits[i].flag};
                RotateDetData(local_hit, gRunConfig);
                TransformDetData(local_hit, gRunConfig);
                GetProjection(local_hit, gRunConfig.hycal_z);
                ev->cl_x[i] = local_hit.x;
                ev->cl_y[i] = local_hit.y;
                ev->cl_z[i] = local_hit.z;
                ev->cl_energy[i] = local_hit.energy;
                ev->cl_center[i] = local_hit.center_id;
                ev->cl_flag[i] = local_hit.flag;
            }

            //decode GEM data and reconstruct GEM hits
        if(!prad1){
            if (gem_sys) {
                gem_sys->Clear();
                gem_sys->ProcessEvent(*ssp_evt);
                if (gem_clusterer)
                    gem_sys->Reconstruct(*gem_clusterer);
            }
            else {
                ev->n_gem_hits = 0;
                std::cerr << "Warning: GEM system not initialized, skipping GEM reconstruction\n";
            }
            auto &all_hits = gem_sys->GetAllHits();
            ev->n_gem_hits = std::min((int)all_hits.size(), prad2::kMaxGemHits);
            for (int i = 0; i < ev->n_gem_hits; i++) {
                auto &h = all_hits[i];
                ev->det_id[i] = h.det_id;
                ev->gem_x_charge[i] = h.x_charge;
                ev->gem_y_charge[i] = h.y_charge;
                ev->gem_x_peak[i] = h.x_peak;
                ev->gem_y_peak[i] = h.y_peak;
                ev->gem_x_size[i] = h.x_size;
                ev->gem_y_size[i] = h.y_size;
                ev->gem_x_mTbin[i] = h.x_max_timebin;
                ev->gem_y_mTbin[i] = h.y_max_timebin;
                //transform the GEM hit positions to the lab coordinate
                GEMHit local_hit = {h.x, h.y, 0.f, static_cast<uint8_t>(h.det_id)};
                RotateDetData(local_hit, gRunConfig);
                TransformDetData(local_hit, gRunConfig);
                ev->gem_x[i] = local_hit.x;
                ev->gem_y[i] = local_hit.y;
                ev->gem_z[i] = local_hit.z;
            }

            // Perform matching between HyCal clusters and GEM hits
            //store all the hits on HyCal and GEMs in this event
            std::vector<HCHit> hc_hits;
            std::vector<GEMHit> gem_hits[4]; // separate vector for each GEM
            for (int i = 0; i < ev->n_clusters; ++i)
                hc_hits.push_back({ev->cl_x[i], ev->cl_y[i], ev->cl_z[i], ev->cl_energy[i], ev->cl_center[i], ev->cl_flag[i]});
            for (int i = 0; i < ev->n_gem_hits; ++i)
                gem_hits[ev->det_id[i]].push_back(GEMHit{ev->gem_x[i], ev->gem_y[i], ev->gem_z[i], ev->det_id[i]});
            
            // already transform to the coordinates

            matching.SetMatchRange(gRunConfig.matching_radius); // matching radius in mm, 15mm default
            matching.SetSquareSelection(gRunConfig.matching_use_square); // square/circular cut
            std::vector<MatchHit> matched_hits = matching.Match(hc_hits, gem_hits[0], gem_hits[1], gem_hits[2], gem_hits[3]);
            std::vector<MatchHit_perChamber> matched_hits_chamber = matching.MatchPerChamber(hc_hits, gem_hits[0], gem_hits[1], gem_hits[2], gem_hits[3]); 
            
            for(int i = 0; i < matched_hits_chamber.size(); i++){
                auto &m = matched_hits_chamber[i];
                int cl_idx = m.hycal_idx;
                if( cl_idx != i) std::cerr << "Warning: cluster index mismatch in matched_hits_chamber: " << cl_idx << " vs " << i << "\n";
                for(int j = 0; j < 4; j++){
                    ev->matchGEMx[i][j] = m.gem_hits[j][0];
                    ev->matchGEMy[i][j] = m.gem_hits[j][1];
                    ev->matchGEMz[i][j] = m.gem_hits[j][2];
                }
                ev->matchFlag[i] = 0;
                ev->matchFlag[i] = m.mflag;
            }

            ev->matchNum = std::min((int)matched_hits.size(), prad2::kMaxClusters);
            for (int i = 0; i < ev->matchNum; i++){
                // save the matched GEM hit (must 2 matchings) info in mHit_ arrays for quick check
                ev->mHit_E[i] = matched_hits[i].hycal_hit.energy;
                ev->mHit_x[i] = matched_hits[i].hycal_hit.x;
                ev->mHit_y[i] = matched_hits[i].hycal_hit.y;
                ev->mHit_z[i] = matched_hits[i].hycal_hit.z;
                for(int j = 0; j < 2; j++) {
                    ev->mHit_gx[i][j] =  matched_hits[i].gem[j].x;
                    ev->mHit_gy[i][j] =  matched_hits[i].gem[j].y;
                    ev->mHit_gz[i][j] =  matched_hits[i].gem[j].z;
                    ev->mHit_gid[i][j] = matched_hits[i].gem[j].det_id; // placeholder for GEM hit ID if needed
                }
            }

        } //end of if(PRad1)
            tree->Fill();
            total++;
            if (total % 1000 == 0)
                std::cerr << "\rReplay: " << total << " events processed" << std::flush;
        }
    }
    std::cerr << "\rReplay: " << total << " events reconstructed -> " << output_root << "\n";
    tree->Write();
    delete outfile;

    return true;
}
} // namespace analysis
