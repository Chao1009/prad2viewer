#pragma once
//=============================================================================
// PhysicsTools.h — physics analysis tools for PRad2
//
// Provides kinematic calculations, energy loss corrections, and
// per-module energy histogram management with ROOT.
// Depends on prad2det (HyCalSystem) and ROOT (TH1F/TH2F).
//=============================================================================

#include "HyCalSystem.h"
#include <TH1F.h>
#include <TH2F.h>
#include <array>
#include <string>
#include <vector>
#include <memory>

namespace analysis {

// two simple data structures used in physics analysis
struct GEMHit {
    float x = 0.f;
    float y = 0.f;
    float z = 0.f;
    uint8_t det_id = 5; // 0-3 for GEM1-GEM4
};

struct HCHit {
    float x = 0.f;
    float y = 0.f;
    float z = 0.f;
    float energy = 0.f;
    uint16_t center_id = 0; // index of central block
    uint32_t flag = -1;
};

//transfrom detector coordinates to target and beam center coordinates
// only used for offline analysis
void TransformDetData(std::vector<HCHit> &hc_hits,  float beamX, float beamY, float ZfromTarget);
void TransformDetData(std::vector<GEMHit> &gem_hits, float beamX, float beamY, float ZfromTarget);

class PhysicsTools
{
public:
    explicit PhysicsTools(fdec::HyCalSystem &hycal);
    ~PhysicsTools();

    // --- per-module cluster energy histograms --------------------------------
    void FillModuleEnergy(int module_id, float energy);
    TH1F *GetModuleEnergyHist(int module_id) const;

    // --- 2D energy vs module index -------------------------------------------
    void FillEnergyVsModule(int module_id, float energy);
    TH2F *GetEnergyVsModuleHist() const { return h2_energy_module_.get(); }

    // --- energy vs scattering angle -------------------------------------------
    void FillEnergyVsTheta(float theta_deg, float energy);
    TH2F *GetEnergyVsThetaHist() const { return h2_energy_theta_.get(); }

    // physics event yield histograms (caller owns the returned histogram)
    std::unique_ptr<TH1F> GetEpYieldHist(TH2F *energy_theta, float Ebeam);
    std::unique_ptr<TH1F> GetEeYieldHist(TH2F *energy_theta, float Ebeam);
    std::unique_ptr<TH1F> GetYieldRatioHist(TH1F *ep_hist, TH1F *ee_hist);

    // --- Moller event Hist ------------------------------------------------
    void Fill2armMollerPosHist(float x, float y);
    TH2F *Get2armMollerPosHist() const { return h2_moller_pos_.get(); }

    // --- peak / resolution analysis ------------------------------------------
    // Returns {peak, resolution} from Gaussian fit. Resolution = sigma/mean.
    std::array<float, 2> FitPeakResolution(int module_id) const;
    void Resolution2Database(int run_id);

    // --- kinematics ----------------------------------------------------------
    // Expected energy for elastic e-p or e-e scattering.
    //   theta: scattering angle in degrees
    //   Ebeam: beam energy in MeV
    //   type:  "ep" or "ee"
    static float ExpectedEnergy(float theta_deg, float Ebeam, const std::string &type);

    // Energy loss correction for electron passing through target + windows.
    //   theta: scattering angle in degrees
    //   E:     measured energy in MeV
    static float EnergyLoss(float theta_deg, float E);

    static float GetShowerDepth(int primex_id, const float &E);

    //physics analysis helpers
    //data structure for storing reconstructed Moller events used for analysis
    struct DataPoint
    {
        float x;
        float y;
        float z;
        float E;

        DataPoint() {};
        DataPoint(float xi, float yi, float zi, float Ei) : x(xi), y(yi), z(zi), E(Ei) {};
    };
    typedef std::pair<DataPoint, DataPoint> MollerEvent;
    typedef std::vector<MollerEvent> MollerData;

    // Get the center of the Moller distribution in x-y space
    // enter two moller events, find the intersection of 2 lines
    // output the x-y coordinates of the center for every 2 moller events
    std::array<float, 2> GetMollerCenter( MollerEvent &event1, MollerEvent &event2);

    float GetMollerZdistance(MollerEvent &event, float Ebeam);

    //Get azimuthal angle difference(should be around 180 degrees) of the Moller event
    float GetMollerPhiDiff(MollerEvent &event1);

    float GetPhiAngle(float x, float y);

    void FillMollerPhiDiff(float phi_diff) { if (moller_phi_diff_) moller_phi_diff_->Fill(phi_diff); }
    void FillMollerXY(float x, float y) { if (moller_x_) moller_x_->Fill(x); if (moller_y_) moller_y_->Fill(y); }
    void FillMollerZ(float z) { if (moller_z_) moller_z_->Fill(z); }

    TH1F *GetMollerPhiDiffHist() const { return moller_phi_diff_.get(); };
    TH1F *GetMollerXHist() const { return moller_x_.get(); };
    TH1F *GetMollerYHist() const { return moller_y_.get(); };
    TH1F *GetMollerZHist() const { return moller_z_.get(); };

    //fill and get gain monitoring replay histograms
    void Fill_lmsCH_lmsHeight(int lms_id, float height)
        { if (lms_id >= 0 && lms_id < 4 && h_lmsCH_lmsHeight_[lms_id]) h_lmsCH_lmsHeight_[lms_id]->Fill(height); }
    void Fill_lmsCH_lmsIntegral(int lms_id, float integral)
        { if (lms_id >= 0 && lms_id < 4 && h_lmsCH_lmsIntegral_[lms_id]) h_lmsCH_lmsIntegral_[lms_id]->Fill(integral); }
    void Fill_lmsCH_alphaHeight(int lms_id, float height)
        { if (lms_id >= 0 && lms_id < 4 && h_lmsCH_alphaHeight_[lms_id]) h_lmsCH_alphaHeight_[lms_id]->Fill(height); }
    void Fill_lmsCH_alphaIntegral(int lms_id, float integral)
        { if (lms_id >= 0 && lms_id < 4 && h_lmsCH_alphaIntegral_[lms_id]) h_lmsCH_alphaIntegral_[lms_id]->Fill(integral); }
    void Fill_modCH_lmsHeight(int module_id, float height)
        { int module_index = hycal_.id_to_index(module_id); 
          if (module_index >= 0 && module_index < (int)h_modCH_lmsHeight_.size()) h_modCH_lmsHeight_[module_index]->Fill(height); }
    void Fill_modCH_lmsIntegral(int module_id, float integral)
        { int module_index = hycal_.id_to_index(module_id);
          if (module_index >= 0 && module_index < (int)h_modCH_lmsIntegral_.size()) h_modCH_lmsIntegral_[module_index]->Fill(integral); }

    TH1F *Get_lmsCH_lmsHeightHist(int lms_id) const 
        { return (lms_id >= 0 && lms_id < 4) ? h_lmsCH_lmsHeight_[lms_id].get() : nullptr; };
    TH1F *Get_lmsCH_lmsIntegralHist(int lms_id) const 
        { return (lms_id >= 0 && lms_id < 4) ? h_lmsCH_lmsIntegral_[lms_id].get() : nullptr; };
    TH1F *Get_lmsCH_alphaHeightHist(int lms_id) const 
        { return (lms_id >= 0 && lms_id < 4) ? h_lmsCH_alphaHeight_[lms_id].get() : nullptr; };
    TH1F *Get_lmsCH_alphaIntegralHist(int lms_id) const 
        { return (lms_id >= 0 && lms_id < 4) ? h_lmsCH_alphaIntegral_[lms_id].get() : nullptr; };
    TH1F *Get_modCH_lmsHeightHist(int module_id) const
        { 
            int module_index = hycal_.id_to_index(module_id);
            return (module_index >= 0 && module_index < (int)h_modCH_lmsHeight_.size()) ? h_modCH_lmsHeight_[module_index].get() : nullptr; 
        };
    TH1F *Get_modCH_lmsIntegralHist(int module_id) const
        { 
            int module_index = hycal_.id_to_index(module_id);
            return (module_index >= 0 && module_index < (int)h_modCH_lmsIntegral_.size()) ? h_modCH_lmsIntegral_[module_index].get() : nullptr; 
        };

private:
    fdec::HyCalSystem &hycal_;
    std::vector<std::unique_ptr<TH1F>> module_hists_;  // one per module
    std::unique_ptr<TH2F> h2_energy_module_;
    std::unique_ptr<TH2F> h2_energy_theta_;
    std::unique_ptr<TH2F> h2_moller_pos_;
    std::unique_ptr<TH1F> moller_phi_diff_;
    std::unique_ptr<TH1F> moller_x_;
    std::unique_ptr<TH1F> moller_y_;
    std::unique_ptr<TH1F> moller_z_;

    //histograms for gain monitoring replay
    std::unique_ptr<TH1F> h_lmsCH_lmsHeight_[4];
    std::unique_ptr<TH1F> h_lmsCH_lmsIntegral_[4];
    std::unique_ptr<TH1F> h_lmsCH_alphaHeight_[4];
    std::unique_ptr<TH1F> h_lmsCH_alphaIntegral_[4];
    std::vector<std::unique_ptr<TH1F>> h_modCH_lmsHeight_; //per module
    std::vector<std::unique_ptr<TH1F>> h_modCH_lmsIntegral_;
};

} // namespace analysis
