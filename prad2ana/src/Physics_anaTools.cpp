
#include "Physics_anaTools.h"

#include <cmath>
#include <array>

#include "TMath.h"
#include "TF1.h"

physics_tools::physics_tools(HyCalSystem &hycal_system) : hycal(hycal_system)
{
    for (int i = 0; i < 1156; ++i) {
        CrystalClusterEHist[i] = new TH1F(Form("CrystalClusterE_W%d", i+1),
                                          Form("Crystal W%d Cluster E;E (MeV);Counts", i+1),
                                          2000, 0, 4000);
    }
    ClusterE2DHist = new TH2F("ClusterE2D", "Cluster Energy vs Theta;Theta (degrees);E (MeV)", 70, 0, 7, 2000, 0, 4000);
    for (int i = 0; i < 900; ++i) {
        LGClusterEHist[i] = new TH1F(Form("LGClusterE_G%d", i+1),
                                     Form("LG G%d Cluster E;E (MeV);Counts", i+1),
                                     2000, 0, 4000);
    }
}

physics_tools::~physics_tools()
{
    for (int i = 0; i < 1156; ++i) {
        delete CrystalClusterEHist[i];
    }
    delete ClusterE2DHist;
    for (int i = 0; i < 900; ++i)
        delete LGClusterEHist[i];
}

void physics_tools::FillModuleClusterEHist(int module_id, float energy)
{
    if (module_id >= PWO_ID0) {
        int idx = module_id - PWO_ID0 - 1;
        if (idx >= 0 && idx < 1156 && CrystalClusterEHist[idx])
            CrystalClusterEHist[idx]->Fill(energy);
    } else {
        int idx = module_id - 1;
        if (idx >= 0 && idx < 900 && LGClusterEHist[idx])
            LGClusterEHist[idx]->Fill(energy);
    }
}

void physics_tools::FillClusterE2DHist(int module_id, float energy)
{
    auto mod = hycal.module_by_id(module_id);
    if (mod && ClusterE2DHist) {
        double theta = atan(sqrt(mod->x*mod->x + mod->y*mod->y)/6225.);
        theta *= 180. / TMath::Pi(); // convert to degrees
        ClusterE2DHist->Fill(float(theta), energy);
    }
}

TH1F* physics_tools::GetModuleClusterEHist(int module_id)
{
    if (module_id >= PWO_ID0) {
        int idx = module_id - PWO_ID0 - 1;
        if (idx >= 0 && idx < 1156) return CrystalClusterEHist[idx];
    } else {
        int idx = module_id - 1;
        if (idx >= 0 && idx < 900) return LGClusterEHist[idx];
    }
    return nullptr;
}

std::array<float, 2> physics_tools::GetPeakAndResolution(int module_id)
{
    TH1F* hist = GetModuleClusterEHist(module_id);
    if (!hist) return {-1.f, -1.f};

    int n_bins = hist->GetNbinsX();
    int max_bin = hist->GetMaximumBin();
    float peak = hist->GetBinCenter(max_bin);
    double fit_min = peak - 2. * 0.025 / sqrt(peak/1000.) * peak;
    double fit_max = peak + 2. * 0.025 / sqrt(peak/1000.) * peak;

    TF1 *fitFunc = new TF1("fitFunc", "gaus", fit_min, fit_max);
    hist->Fit(fitFunc, "RQ");

    peak = fitFunc->GetParameter(1);
    float sigma = fitFunc->GetParameter(2);
    delete fitFunc;

    std::array<float, 2> result = {peak, sigma};
    return result;
}

float physics_tools::GetElossIonElectron(float &theta, float& E)
{
    // Calculates energy loss dE/dx in MeV/mm due to ionization for relativistic electrons/positrons.
    //
    // For formula used, see:
    // Track fitting with energy loss
    // Stampfer, Regler and Fruehwirth
    // Computer Physics Communication 79 (1994), 157-164
    //
    // ZoverA:     atomic number / atomic mass of passed material
    // density:    density of material in g/mm^3
    // I : mean excitation energy in MeV

    //only the Al thin window now, need to add GEM, GEM frame, and the cover between HyCal and GEM to be exact

    float ZoverA[3]   = { 13./27., 10.6/21.8 , 0.49919};
    float density[3]  = {2.699, 0.1117433, 1.205e-3};                      // g/cm^2
    float I[3]        = {166*1.e-6, 106.6e-6, 85.7e-6};                       // MeV
    float de       = 5.0989 * 1.e-25;                 // 4*pi*re*me*c^2 in MeV * mm^2 with re = classical electron radius
    float avogadro = TMath::Na();                     // Avogadro constant in 1/mol
    float me       = 0.5109989181;   // electron mass in MeV/c^2
    float gamma    = E / me;                          // Relativistic gamma-factor.

    // Formula is slightly different for electrons and positrons.
    float gammaFac = 3.;
    float corr     = 1.95;

    float eDep = 0;
    float length[3] = { 0.2, 1.5, 50 };
    for (int i=0; i<3; i++) {
      length[i] /= cos(theta);
      float dedx = 0.5 * de * avogadro * density[i] * ZoverA[i] * (2 * TMath::Log(2*me/I[i]) + gammaFac * TMath::Log(gamma) - corr);
      eDep += dedx*length[i];
    }

    return eDep;
}

//new for expected energy calculation of e-p and e-e on each module
//unit : MeV
float physics_tools::GetExpectedEnergy(int primex_id, float Ebeam, std::string type)
{   
    const Module *m = hycal.module_by_id(primex_id);
    if (!m) return -1;
    float crystal_z = 5817.; // distance from target to crystal front face in mm
    float theta = atan(sqrt(m->x*m->x + m->y*m->y)/crystal_z);
    float expectE = 0.;
    if(type == "e-p")
        expectE = Ebeam*938.272046 / ( Ebeam*(1.-cos(theta)) + 938.272046 );
    else if(type == "e-e")
        expectE = Ebeam / ( 1. + (Ebeam/0.511)*(1.-cos(theta)) );
    else
        return -1;
    float eLoss = GetElossIonElectron(theta, expectE);
    //eloss->Fill(expectE, eLoss);
    return expectE - eLoss;
}