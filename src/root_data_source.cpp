// =========================================================================
// root_data_source.cpp — ROOT file data source implementations
// =========================================================================

#ifdef WITH_ROOT

#include "root_data_source.h"

#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>

#include <iostream>
#include <algorithm>

// =========================================================================
// Factory: detect tree type
// =========================================================================

std::unique_ptr<DataSource> createRootDataSource(
    const std::string &path,
    const std::unordered_map<int, uint32_t> &crate_to_roc,
    const fdec::HyCalSystem *hycal)
{
    std::unique_ptr<TFile> f(TFile::Open(path.c_str(), "READ"));
    if (!f || f->IsZombie()) return nullptr;

    if (f->Get<TTree>("events"))
        return std::make_unique<RootRawDataSource>(crate_to_roc, hycal);
    if (f->Get<TTree>("recon"))
        return std::make_unique<RootReconDataSource>();

    return nullptr;
}

// =========================================================================
// RootRawDataSource
// =========================================================================

std::string RootRawDataSource::open(const std::string &path)
{
    close();
    file_.reset(TFile::Open(path.c_str(), "READ"));
    if (!file_ || file_->IsZombie()) {
        file_.reset();
        return "cannot open ROOT file";
    }
    tree_ = file_->Get<TTree>("events");
    if (!tree_) { close(); return "no 'events' tree in ROOT file"; }

    n_entries_ = static_cast<int>(tree_->GetEntries());

    // event-level branches (shared RawEventData struct)
    tree_->SetBranchAddress("event_num",    &ev_.event_num);
    if (tree_->GetBranch("trigger_type"))
        tree_->SetBranchAddress("trigger_type", &ev_.trigger_type);
    if (tree_->GetBranch("trigger_bits"))
        tree_->SetBranchAddress("trigger_bits", &ev_.trigger_bits);
    tree_->SetBranchAddress("timestamp",     &ev_.timestamp);

    // HyCal per-channel data (indexed by module_id, not crate/slot/channel)
    tree_->SetBranchAddress("hycal.nch",        &ev_.nch);
    tree_->SetBranchAddress("hycal.module_id",  ev_.module_id);
    tree_->SetBranchAddress("hycal.nsamples",   ev_.nsamples);
    tree_->SetBranchAddress("hycal.samples",    ev_.samples);
    tree_->SetBranchAddress("hycal.ped_mean",   ev_.ped_mean);
    tree_->SetBranchAddress("hycal.ped_rms",    ev_.ped_rms);
    tree_->SetBranchAddress("hycal.integral",   ev_.integral);

    has_peaks_ = (tree_->GetBranch("hycal.npeaks") != nullptr);
    if (has_peaks_) {
        tree_->SetBranchAddress("hycal.npeaks",        ev_.npeaks);
        tree_->SetBranchAddress("hycal.peak_height",   ev_.peak_height);
        tree_->SetBranchAddress("hycal.peak_time",     ev_.peak_time);
        tree_->SetBranchAddress("hycal.peak_integral", ev_.peak_integral);
    }

    has_gem_ = (tree_->GetBranch("gem.nch") != nullptr);
    if (has_gem_) {
        tree_->SetBranchAddress("gem.nch",         &ev_.gem_nch);
        tree_->SetBranchAddress("gem.mpd_crate",   ev_.mpd_crate);
        tree_->SetBranchAddress("gem.mpd_fiber",   ev_.mpd_fiber);
        tree_->SetBranchAddress("gem.apv",         ev_.apv);
        tree_->SetBranchAddress("gem.strip",       ev_.strip);
        tree_->SetBranchAddress("gem.ssp_samples", ev_.ssp_samples);
    }

    if (!hycal_) {
        std::cerr << "ROOT raw: warning — no HyCalSystem provided; "
                     "module_id → DAQ mapping unavailable, HyCal channels will be skipped\n";
    }

    std::cerr << "ROOT raw: " << n_entries_ << " events"
              << (has_peaks_ ? ", peaks" : "")
              << (has_gem_   ? ", GEM"   : "") << "\n";
    return "";
}

void RootRawDataSource::close()
{
    tree_ = nullptr;
    file_.reset();
    n_entries_ = 0;
    has_peaks_ = false;
    has_gem_ = false;
}

DataSourceCaps RootRawDataSource::capabilities() const
{
    return {
        true,       // has_waveforms
        has_peaks_, // has_peaks
        true,       // has_pedestals
        true,       // has_clusters (computed)
        has_gem_,   // has_gem_raw
        has_gem_,   // has_gem_hits (computed)
        false,      // has_epics
        false,      // has_sync
        "root_raw"
    };
}

void RootRawDataSource::fillEventData(fdec::EventData &evt) const
{
    evt.clear();
    evt.info.event_number = ev_.event_num;
    evt.info.trigger_type = ev_.trigger_type;
    evt.info.trigger_bits = ev_.trigger_bits;
    evt.info.timestamp = static_cast<uint64_t>(ev_.timestamp);

    // Without the HyCalSystem we cannot reverse module_id → (crate, slot, ch).
    // Downstream code indexes by ROC/slot/channel, so skip HyCal channels.
    if (!hycal_) return;

    for (int i = 0; i < ev_.nch && i < prad2::kMaxChannels; ++i) {
        const fdec::Module *mod = hycal_->module_by_id(ev_.module_id[i]);
        if (!mod) continue;

        int crate = mod->daq.crate;
        int sl    = mod->daq.slot;
        int ch    = mod->daq.channel;
        if (sl < 0 || sl >= fdec::MAX_SLOTS || ch < 0 || ch >= fdec::MAX_CHANNELS) continue;

        // translate crate ID → ROC tag
        uint32_t roc_tag = static_cast<uint32_t>(crate);
        auto it = crate_to_roc_.find(crate);
        if (it != crate_to_roc_.end()) roc_tag = it->second;

        int roc_idx = -1;
        for (int r = 0; r < evt.nrocs; ++r) {
            if (evt.rocs[r].tag == roc_tag) { roc_idx = r; break; }
        }
        if (roc_idx < 0) {
            if (evt.nrocs >= fdec::MAX_ROCS) continue;
            roc_idx = evt.nrocs++;
            evt.rocs[roc_idx].present = true;
            evt.rocs[roc_idx].tag = roc_tag;
        }

        auto &slot = evt.rocs[roc_idx].slots[sl];
        slot.present = true;
        slot.channel_mask |= (1ull << ch);
        auto &cd = slot.channels[ch];
        cd.nsamples = std::min((int)ev_.nsamples[i], fdec::MAX_SAMPLES);
        for (int s = 0; s < cd.nsamples; ++s)
            cd.samples[s] = ev_.samples[i][s];
    }
}

std::string RootRawDataSource::decodeEvent(int index, fdec::EventData &evt,
                                            ssp::SspEventData *ssp)
{
    if (index < 0 || index >= n_entries_) return "event out of range";
    std::lock_guard<std::mutex> lk(mtx_);
    tree_->GetEntry(index);
    fillEventData(evt);
    if (ssp) ssp->clear();
    return "";
}

void RootRawDataSource::iterateAll(EventCallback ev_cb, ReconCallback /*recon_cb*/,
                                    ControlCallback /*ctrl_cb*/, EpicsCallback /*epics_cb*/,
                                    DscCallback /*dsc_cb*/, int /*dsc_bank_tag*/)
{
    if (!tree_ || !ev_cb) return;

    std::lock_guard<std::mutex> lk(mtx_);  // block concurrent decodeEvent calls
    auto event_ptr = std::make_unique<fdec::EventData>();
    for (int i = 0; i < n_entries_; ++i) {
        tree_->GetEntry(i);
        fillEventData(*event_ptr);
        ev_cb(i, *event_ptr, nullptr);
    }
}

// =========================================================================
// RootReconDataSource
// =========================================================================

std::string RootReconDataSource::open(const std::string &path)
{
    close();
    file_.reset(TFile::Open(path.c_str(), "READ"));
    if (!file_ || file_->IsZombie()) {
        file_.reset();
        return "cannot open ROOT file";
    }
    tree_ = file_->Get<TTree>("recon");
    if (!tree_) { close(); return "no 'recon' tree in ROOT file"; }

    n_entries_ = static_cast<int>(tree_->GetEntries());

    tree_->SetBranchAddress("event_num",    &ev_.event_num);
    if (tree_->GetBranch("trigger_type"))
        tree_->SetBranchAddress("trigger_type", &ev_.trigger_type);
    if (tree_->GetBranch("trigger_bits"))
        tree_->SetBranchAddress("trigger_bits", &ev_.trigger_bits);
    tree_->SetBranchAddress("timestamp",    &ev_.timestamp);

    // HyCal clusters
    tree_->SetBranchAddress("n_clusters",   &ev_.n_clusters);
    tree_->SetBranchAddress("cl_x",         ev_.cl_x);
    tree_->SetBranchAddress("cl_y",         ev_.cl_y);
    tree_->SetBranchAddress("cl_energy",    ev_.cl_energy);
    tree_->SetBranchAddress("cl_nblocks",   ev_.cl_nblocks);
    tree_->SetBranchAddress("cl_center",    ev_.cl_center);

    // GEM hits
    tree_->SetBranchAddress("n_gem_hits",   &ev_.n_gem_hits);
    tree_->SetBranchAddress("det_id",       ev_.det_id);
    tree_->SetBranchAddress("gem_x",        ev_.gem_x);
    tree_->SetBranchAddress("gem_y",        ev_.gem_y);
    tree_->SetBranchAddress("gem_x_charge", ev_.gem_x_charge);
    tree_->SetBranchAddress("gem_y_charge", ev_.gem_y_charge);
    tree_->SetBranchAddress("gem_x_peak",   ev_.gem_x_peak);
    tree_->SetBranchAddress("gem_y_peak",   ev_.gem_y_peak);
    tree_->SetBranchAddress("gem_x_size",   ev_.gem_x_size);
    tree_->SetBranchAddress("gem_y_size",   ev_.gem_y_size);

    std::cerr << "ROOT recon: " << n_entries_ << " events\n";
    return "";
}

void RootReconDataSource::close()
{
    tree_ = nullptr;
    file_.reset();
    n_entries_ = 0;
}

DataSourceCaps RootReconDataSource::capabilities() const
{
    return {
        false,        // has_waveforms
        false,        // has_peaks
        false,        // has_pedestals
        true,         // has_clusters (pre-computed)
        false,        // has_gem_raw
        true,         // has_gem_hits (pre-computed)
        false,        // has_epics
        false,        // has_sync
        "root_recon"
    };
}

std::string RootReconDataSource::decodeEvent(int index, fdec::EventData &evt,
                                              ssp::SspEventData *ssp)
{
    if (index < 0 || index >= n_entries_) return "event out of range";
    std::lock_guard<std::mutex> lk(mtx_);
    tree_->GetEntry(index);
    evt.clear();
    evt.info.event_number = ev_.event_num;
    evt.info.trigger_type = ev_.trigger_type;
    evt.info.trigger_bits = ev_.trigger_bits;
    evt.info.timestamp = static_cast<uint64_t>(ev_.timestamp);
    if (ssp) ssp->clear();
    return "";
}

void RootReconDataSource::fillRecon(ReconEventData &recon) const
{
    recon.event_num = ev_.event_num;
    recon.trigger_type = ev_.trigger_type;
    recon.trigger_bits = ev_.trigger_bits;
    recon.timestamp = static_cast<uint64_t>(ev_.timestamp);
    recon.clusters.clear();
    for (int i = 0; i < ev_.n_clusters && i < prad2::kMaxClusters; ++i)
        recon.clusters.push_back({ev_.cl_x[i], ev_.cl_y[i], ev_.cl_energy[i],
                                   static_cast<int>(ev_.cl_nblocks[i]),
                                   static_cast<int>(ev_.cl_center[i])});
    recon.gem_hits.clear();
    for (int i = 0; i < ev_.n_gem_hits && i < prad2::kMaxGemHits; ++i)
        recon.gem_hits.push_back({static_cast<int>(ev_.det_id[i]),
                                   ev_.gem_x[i], ev_.gem_y[i],
                                   ev_.gem_x_charge[i], ev_.gem_y_charge[i],
                                   ev_.gem_x_peak[i], ev_.gem_y_peak[i],
                                   static_cast<int>(ev_.gem_x_size[i]),
                                   static_cast<int>(ev_.gem_y_size[i])});
}

bool RootReconDataSource::decodeReconEvent(int index, ReconEventData &recon)
{
    if (index < 0 || index >= n_entries_) return false;
    std::lock_guard<std::mutex> lk(mtx_);
    tree_->GetEntry(index);
    fillRecon(recon);
    return true;
}

void RootReconDataSource::iterateAll(EventCallback /*ev_cb*/, ReconCallback recon_cb,
                                      ControlCallback /*ctrl_cb*/, EpicsCallback /*epics_cb*/,
                                      DscCallback /*dsc_cb*/, int /*dsc_bank_tag*/)
{
    if (!tree_ || !recon_cb) return;

    std::lock_guard<std::mutex> lk(mtx_);  // block concurrent decodeReconEvent calls
    ReconEventData recon;
    for (int i = 0; i < n_entries_; ++i) {
        tree_->GetEntry(i);
        fillRecon(recon);
        recon_cb(i, recon);
    }
}

#endif // WITH_ROOT
