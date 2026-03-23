#pragma once
//=============================================================================
// Physics_anaTools.h — physics analysis tools for PRad2
//=============================================================================

#include "HyCalSystem.h"

#include "TH1F.h"
#include "TH2F.h"

using namespace fdec;

class physics_tools
{
    public:
        physics_tools(HyCalSystem &hycal_system);
        ~physics_tools();
        
        float GetElossIonElectron(float &theta, float& E);

        //physics histograms
        void FillModuleClusterEHist(int module_id, float energy);
        void FillClusterE2DHist(int module_id, float energy);
        TH1F* GetModuleClusterEHist(int module_id);
        TH2F* GetClusterE2DHist() { return ClusterE2DHist; };

        //about the reconstruction resolution monitoring
        std::array<float, 2> GetPeakAndResolution(int module_id);
        void Rseolution2database(int module_id, int run_id, float peak, float resolution);
        void ResolutionHistory(int module_id);

        //physics analysis tools
        float GetExpectedEnergy(int primex_id, float Ebeam, std::string type);


    private:
        HyCalSystem &hycal;
        
        TH1F* CrystalClusterEHist[1156];
        TH1F* LGClusterEHist[900];
        TH2F* ClusterE2DHist;
};