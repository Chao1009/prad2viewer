//=============================================================================
// PhysicsTools.cpp — physics analysis tools
//=============================================================================

#include "PhysicsTools.h"
#include <TF1.h>
#include <TMath.h>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>

#ifndef DATABASE_DIR
#define DATABASE_DIR "."
#endif

namespace analysis {

// Physical constants
static constexpr float M_PROTON  = 938.272f;   // MeV
static constexpr float M_ELECTRON = 0.511f;    // MeV
static constexpr float DEG2RAD = 3.14159265f / 180.f;

PhysicsTools::PhysicsTools(fdec::HyCalSystem &hycal)
    : hycal_(hycal)
{
    int nmod = hycal_.module_count();
    module_hists_.resize(nmod);
    for (int i = 0; i < nmod; ++i) {
        auto &mod = hycal_.module(i);
        std::string name = "h_" + mod.name;
        std::string title = mod.name + " cluster energy;Energy (MeV);Counts";
        module_hists_[i] = std::make_unique<TH1F>(name.c_str(), title.c_str(), 300, 0, 3000);
    }
    h2_energy_module_ = std::make_unique<TH2F>(
        "h2_energy_module", "Energy vs Module;Module Index;Energy (MeV)",
        nmod, 0, nmod, 2000, 0, 4000);
    h2_energy_theta_ = std::make_unique<TH2F>(
        "h2_energy_theta", "Energy vs Theta;Theta (deg);Energy (MeV)",
        80, 0, 8, 2000, 0, 4000);

    h2_moller_pos_ = std::make_unique<TH2F>(
        "h2_moller_pos", "Moller 2-arm Hit Position;X (mm);Y (mm)",
        200, -500, 500, 200, -500, 500);

    moller_phi_diff_ = std::make_unique<TH1F>(
        "h_moller_phi_diff", "Moller Phi Difference;Phi_{e1} - Phi_{e2} (deg);Counts",
        40, -20, 20);
    moller_x_ = std::make_unique<TH1F>(
        "h_moller_x", "Moller Center X Position (HyCal);X (mm);Counts",
        100, -10, 10);
    moller_y_ = std::make_unique<TH1F>(
        "h_moller_y", "Moller Center Y Position (HyCal);Y (mm);Counts",
        100, -10, 10);
    moller_z_ = std::make_unique<TH1F>(
        "h_moller_z", "Moller Z Position (HyCal);Z (mm);Counts",
        1000, 5000, 8000);
}

PhysicsTools::~PhysicsTools() = default;

// transfrom detector coordinates to target and beam center coordinates
// only used for offline analysis
void TransformDetData(std::vector<HCHit> &hc_hits, float beamX, float beamY, float ZfromTarget)
{
    // Transform HyCal hit coordinates from detector frame to target frame
    for (auto &hc_hit : hc_hits) {
        hc_hit.x -= beamX;
        hc_hit.y -= beamY;
        hc_hit.z += ZfromTarget;
    }
}

void TransformDetData(std::vector<GEMHit> &gem_hits, float beamX, float beamY, float ZfromTarget)
{
    // Transform GEM hit coordinates from detector frame to target frame
    for (auto &gem_hit : gem_hits) {
        gem_hit.x -= beamX;
        gem_hit.y -= beamY;
        gem_hit.z += ZfromTarget;
    }
}

void PhysicsTools::FillModuleEnergy(int module_index, float energy)
{
    if (module_index >= 0 && module_index < (int)module_hists_.size())
        module_hists_[module_index]->Fill(energy);
}

TH1F *PhysicsTools::GetModuleEnergyHist(int module_index) const
{
    if (module_index >= 0 && module_index < (int)module_hists_.size())
        return module_hists_[module_index].get();
    return nullptr;
}

void PhysicsTools::FillEnergyVsModule(int module_index, float energy)
{
    if (h2_energy_module_)
        h2_energy_module_->Fill(module_index, energy);
}

void PhysicsTools::FillEnergyVsTheta(float theta_deg, float energy)
{
    if (h2_energy_theta_)
        h2_energy_theta_->Fill(theta_deg, energy);
}

std::unique_ptr<TH1F> PhysicsTools::GetEpYieldHist(TH2F *energy_theta, float Ebeam)
{
    if (!energy_theta) return nullptr;

    auto h_ep = std::make_unique<TH1F>("ep_yield", "Elastic e-p Yield;Scattering Angle (deg);Counts", 80, 0, 8);
    h_ep->SetDirectory(nullptr);
    for (int i = 1; i <= energy_theta->GetNbinsX(); i++) {
        for (int j = 1; j <= energy_theta->GetNbinsY(); j++) {
            float theta = energy_theta->GetXaxis()->GetBinCenter(i);
            float E = energy_theta->GetYaxis()->GetBinCenter(j);
            float E_expected = ExpectedEnergy(theta, Ebeam, "ep");
            if (std::abs(E - E_expected) < E_expected*0.026f/std::sqrt(E_expected/1000.f)) {
                float count = energy_theta->GetBinContent(i, j);
                h_ep->Fill(theta, count);
            }
        }
    }
    return h_ep;
}

std::unique_ptr<TH1F> PhysicsTools::GetEeYieldHist(TH2F *energy_theta, float Ebeam)
{
    if (!energy_theta) return nullptr;

    auto h_ee = std::make_unique<TH1F>("ee_yield", "Elastic e-e Yield;Scattering Angle (deg);Counts", 80, 0, 8);
    h_ee->SetDirectory(nullptr);
    for (int i = 1; i <= energy_theta->GetNbinsX(); i++) {
        for (int j = 1; j <= energy_theta->GetNbinsY(); j++) {
            float theta = energy_theta->GetXaxis()->GetBinCenter(i);
            float E = energy_theta->GetYaxis()->GetBinCenter(j);
            float E_expected = ExpectedEnergy(theta, Ebeam, "ee");
            if (std::abs(E - E_expected) < E_expected*0.026f/std::sqrt(E_expected/1000.f)) {
                float count = energy_theta->GetBinContent(i, j);
                h_ee->Fill(theta, count);
            }
        }
    }
    return h_ee;
}

std::unique_ptr<TH1F> PhysicsTools::GetYieldRatioHist(TH1F *ep_hist, TH1F *ee_hist)
{
    if (!ep_hist || !ee_hist) return nullptr;

    auto h_ratio = std::make_unique<TH1F>("yield_ratio", "Yield Ratio (e-p / e-e);Scattering Angle (deg);Ratio", 80, 0, 8);
    h_ratio->SetDirectory(nullptr);
    for (int i = 1; i <= ep_hist->GetNbinsX(); i++) {
        float theta = ep_hist->GetXaxis()->GetBinCenter(i);
        float ep_count = ep_hist->GetBinContent(i);
        float ee_count = ee_hist->GetBinContent(i);
        if (ee_count > 0) {
            h_ratio->Fill(theta, ep_count / ee_count);
        }
    }
    return h_ratio;
}

void PhysicsTools::Fill2armMollerPosHist(float x, float y)
{
    if (h2_moller_pos_)
        h2_moller_pos_->Fill(x, y);
}

std::array<float, 2> PhysicsTools::FitPeakResolution(int module_index) const
{
    if (module_index < 0 || module_index >= (int)module_hists_.size())
        return {0.f, 0.f};

    TH1F *h = module_hists_[module_index].get();
    if (!h || h->GetEntries() < 50) return {0.f, 0.f};

    // find peak bin, fit Gaussian around it
    int maxbin = h->GetMaximumBin();
    float peak = h->GetBinCenter(maxbin);
    float rms  = h->GetRMS();

    TF1 gaus("gfit", "gaus", peak - 2 * rms, peak + 2 * rms);
    h->Fit(&gaus, "QNR");

    float mean  = gaus.GetParameter(1);
    float sigma = gaus.GetParameter(2);
    float resolution = (mean > 0) ? sigma / mean : 0.f;
    return {mean, resolution};
}

void PhysicsTools::Resolution2Database(int run_id)
{
    std::string db_dir = DATABASE_DIR;
    if (const char *env = std::getenv("PRAD2_DATABASE_DIR"))  db_dir = env;
    std::string filename = db_dir + Form("/recon/run_%d.dat", run_id);

    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "Failed to open file for writing: " << filename << std::endl;
        return;
    }

    int module_count = hycal_.module_count();
    for (int m = 0; m < module_count; m++) {
        auto [peak, sigma] = FitPeakResolution(m);
        if (peak > 0 && sigma > 0) {
            std::string name = hycal_.module(m).name;
            out << name << " " << peak << " " << sigma << "\n";
        }
    }
}

float PhysicsTools::ExpectedEnergy(float theta_deg, float Ebeam, const std::string &type)
{
    float theta = theta_deg * DEG2RAD;
    float cos_t = std::cos(theta);
    float sin_t = std::sin(theta);

    if (type == "ep") {
        // elastic e-p: E' = E * M / (M + E*(1 - cos_t))
        // where M = proton mass
        float expectE = Ebeam * M_PROTON / (M_PROTON + Ebeam * (1.f - cos_t));
        float eloss = EnergyLoss(theta_deg, expectE);
        return expectE - eloss;
    }
    if (type == "ee") {
        // Moller scattering: E' = E * cos^2(theta) / (1 + (E/m)(sin^2(theta)))
        // simplified from CM frame kinematics
        float gamma = Ebeam / M_ELECTRON;
        float num = (gamma + 1.f) * cos_t * cos_t;
        float den = (gamma + 1.f) - (gamma - 1.f) * cos_t * cos_t;
        if (den <= 0) return 0.f;
        float expectE = M_ELECTRON * num / den;
        float eloss = EnergyLoss(theta_deg, expectE);
        return expectE - eloss;
    }
    return 0.f;
}

float PhysicsTools::EnergyLoss(float theta_deg, float E)
{
    // simplified energy loss through target materials
    // path lengths scale as 1/cos(theta) for small angles
    float theta = theta_deg * DEG2RAD;
    float sec = (std::cos(theta) > 0.01f) ? (1.f / std::cos(theta)) : 100.f;

    // material thicknesses (mm) and dE/dx (MeV/mm) — approximate values
    // aluminum window: 0.025 mm, dE/dx ~ 1.6 MeV/mm
    // GEM foils: ~0.05 mm effective, dE/dx ~ 2.0 MeV/mm
    // kapton window: ~0.05 mm, dE/dx ~ 1.8 MeV/mm
    float eloss = 0.f;
    eloss += 0.025f * 1.6f * sec;  // Al window
    eloss += 0.050f * 2.0f * sec;  // GEM
    eloss += 0.050f * 1.8f * sec;  // kapton cover

    return eloss;  // total energy loss in MeV
}

// get shower depth, unit is in MeV
float PhysicsTools::GetShowerDepth(int primex_id, const float &E)
{
    if(E > 0.) {
        // here all the values are hard coded, because these are all physical
        // values corresponding to the material, so no need to change
        // it returns the maximum shower depth that
        // t = X0*(ln(E0/Ec) - Cf),
        // where X0 is radiation length, Ec is critical energy, Cf = -0.5 for
        // electron induced shower and 0.5 for photon
        // units are in mm and MeV
        if(primex_id >= fdec::PWO_ID0) //module_type PbWO4
            return 8.6*(log(E/1.1) - 0.5);

        // -101.2 is the surface difference between Lead Glass and Lead Tungstate modules
        if(primex_id < fdec::PWO_ID0) //module_type PbGlass
            return 26.7*(log(E/2.84) - 0.5);
    }

    return 0.;
}

std::array<float, 2> PhysicsTools::GetMollerCenter(MollerEvent &event1, MollerEvent &event2)
{
    float x1[2], y1[2];
    float x2[2], y2[2];

    x1[0] = event1.first.x; y1[0] = event1.first.y;
    x1[1] = event1.second.x; y1[1] = event1.second.y;
    x2[0] = event2.first.x; y2[0] = event2.first.y;
    x2[1] = event2.second.x; y2[1] = event2.second.y;

    //two lines: y = ax + b, y = cx + d
    float dx1 = x1[0] - x1[1];
    float dx2 = x2[0] - x2[1];
    if (std::abs(dx1) < 1e-6f || std::abs(dx2) < 1e-6f)
        return {0.f, 0.f};  // vertical line — degenerate

    float a = (y1[0] - y1[1]) / dx1;
    float b = y1[0] - a * x1[0];
    float c = (y2[0] - y2[1]) / dx2;
    float d = y2[0] - c * x2[0];

    if (std::abs(a - c) < 1e-6f)
        return {0.f, 0.f};  // parallel lines — no intersection

    float x_cross = (d - b) / (a - c);
    float y_cross = a * x_cross + b;

    return {x_cross, y_cross};

}

float PhysicsTools::GetMollerZdistance(MollerEvent &event, float Ebeam)
{
    float R1 = sqrt(event.first.x*event.first.x + event.first.y*event.first.y);
    float R2 = sqrt(event.second.x*event.second.x + event.second.y*event.second.y);
    float z = sqrt( (Ebeam + M_ELECTRON) * R1 * R2 / (2.*M_ELECTRON) );
    return z;
}

float PhysicsTools::GetMollerPhiDiff(MollerEvent &event1)
{
    // Calculate the azimuthal angle difference (phi) for a Moller event
    float x1 = event1.first.x, y1 = event1.first.y;
    float x2 = event1.second.x, y2 = event1.second.y;
    float phi1 = GetPhiAngle(x1, y1);
    float phi2 = GetPhiAngle(x2, y2);
    float phi_diff = fabs(phi1 - phi2) - 180.f; // Expecting back-to-back, so difference should be around 180 degrees
    return phi_diff;
}

float PhysicsTools::GetPhiAngle(float x, float y)
{
    // atan2 handles all quadrants and x==0 correctly
    float phi = std::atan2(y, x) * 180.f / static_cast<float>(TMath::Pi());
    if (phi < 0) phi += 360.f;
    return phi;
}

} // namespace analysis
