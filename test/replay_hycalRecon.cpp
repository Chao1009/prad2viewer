#include "HyCalSystem.h"
#include "HyCalCluster.h"
#include "load_daq_config.h"
#include "replay.h"

#include <TFile.h>
#include <TTree.h>
#include <iostream>

/*data structure for reconstructed hits:
struct ClusterHit {
int   center_id;    // PrimEx ID of center module
float x, y;         // reconstructed position (mm)
float energy;       // total cluster energy (MeV)
int   nblocks;      // number of modules in cluster
int   npos;         // number of modules used in position reconstruction
uint32_t flag;      // cluster flags
*/
//variables for root tree
const int MAX_HITS = 20;
int reco_event_num, reco_trigger, Nhits, hit_cid[MAX_HITS];
Long64_t reco_timestamp;
float hitE[MAX_HITS], hitX[MAX_HITS], hitY[MAX_HITS];
int hit_nblocks[MAX_HITS], hit_npos[MAX_HITS];
uint32_t hit_flag[MAX_HITS];

void setupTTreeRecon(TTree *tree){
    tree->Branch("event_num", &reco_event_num, "event_num/I");
    tree->Branch("trigger", &reco_trigger, "trigger/I");
    tree->Branch("timestamp", &reco_timestamp, "timestamp/L");
    tree->Branch("Nhits", &Nhits, "Nhits/I");
    tree->Branch("cid", &hit_cid[0], "cid[Nhits]/I");
    tree->Branch("E", &hitE[0], "E[Nhits]/F");
    tree->Branch("X", &hitX[0], "X[Nhits]/F");
    tree->Branch("Y", &hitY[0], "Y[Nhits]/F");
    tree->Branch("nblocks", &hit_nblocks[0], "nblocks[Nhits]/I");
    tree->Branch("npos", &hit_npos[0], "npos[Nhits]/I");
    tree->Branch("flag", &hit_flag[0], "flag[Nhits]/i");
}

static void usage(const char *prog)
{
    std::cerr << "Usage:\n"
              << "  " << prog << " <input_root_file> <output_root_file>\n"
              << "  Options:\n";
}

int main(int argc, char *argv[]){

    if (argc < 3) {
        usage(argv[0]);
        return 1;
    }
    std::string input_root_file = argv[1];
    std::string output_root_file = argv[2];

    //load reconstruction config
    std::string db_dir = std::string(DATABASE_DIR);
    std::string calib_file = std::string(db_dir + "/prad1/prad_calibration.json");
    fdec::ClusterConfig cluster_cfg;
    
    HyCalSystem sys;
    if(!sys.Init(db_dir + "/hycal_modules.json", db_dir + "/daq_map.json")) {
        std::cerr << "Failed to initialize HyCalSystem with database files.\n";
        return 1;
    }
    HyCalCluster clusterer(sys);
    clusterer.SetConfig(cluster_cfg);

    int nmatched = sys.LoadCalibration(calib_file);
    if (nmatched >= 0)
        std::cerr << "Calibration: " << calib_file << " (" << nmatched << " modules)\n";

    //read raw data from root, event loop, fill clusterer, form clusters, print results
    TFile *file_in = TFile::Open(input_root_file.c_str(), "READ");
     if (!file_in || file_in->IsZombie()) {
        std::cerr << "Failed to open input ROOT file: " << input_root_file << "\n";
        return 1;
    }
    TTree *tree = file_in->Get<TTree>("T");  // ROOT 6 style
    if (!tree) {
        std::cerr << "Failed to get TTree from input ROOT file: " << input_root_file << "\n";
        return 1;
    }
    replay rep;
    rep.setupTTreeRead(tree, true); // true: read peaks, false: read raw samples

    // Prepare output ROOT file and tree for reconstructed hits
    TFile *file_out = new TFile(output_root_file.c_str(), "RECREATE");
    TTree *tree_out = new TTree("T", "T");
    setupTTreeRecon(tree_out);

    int nentries = tree->GetEntries();
    for (int i = 0; i < nentries; ++i) {
        tree->GetEntry(i);
        if(i % 1000 == 0) {
            std::cout << "Processing event " << reco_event_num << " / total " << nentries << "\r";
            std::cout.flush();
        }
        clusterer.Clear();
        for (int ch = 0; ch < Nchs; ch++) {
            if (module_id[ch] <= 0) continue; // skip unknown modules
            const fdec::Module *mod = sys.module_by_id(module_id[ch]);
            if (!mod)   continue;
            int module_index = mod->index;
            //float ch_energy = integral[ch]; // or use peak info if needed
            float ch_energy = (float)mod->energize(integral[ch]);
            clusterer.AddHit(module_index, ch_energy);
        }
        clusterer.FormClusters();

        std::vector<fdec::ClusterHit> reco_hits;
        clusterer.ReconstructHits(reco_hits);

        Nhits = (int)reco_hits.size();
        if (Nhits > MAX_HITS) {
            std::cerr << "Warning: number of reconstructed hits (" << Nhits << ") exceeds MAX_HITS (" << MAX_HITS << "). Truncating.\n";
            Nhits = MAX_HITS;
        }
        for(int h = 0; h < Nhits; h++){
            hit_cid[h] = reco_hits[h].center_id;
            hitE[h] = reco_hits[h].energy;
            hitX[h] = reco_hits[h].x;
            hitY[h] = reco_hits[h].y;
            hit_nblocks[h] = reco_hits[h].nblocks;
            hit_npos[h] = reco_hits[h].npos;
            hit_flag[h] = reco_hits[h].flag;
        }
        reco_event_num = event_num;
        reco_trigger = ftrigger;
        reco_timestamp = ftimestamp;

        tree_out->Fill();
    }
    file_out->cd();
    file_out->Write();
    file_out->Close();
    file_in->cd();
    file_in->Close();
    std::cout << "Processing event " << reco_event_num << " / total " << nentries << "\n";
    std::cout << "Finished processing " << nentries << " events. Output saved to " << output_root_file << "\n";
    return 0;
}