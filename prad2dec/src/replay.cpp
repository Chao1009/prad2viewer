
//=============================================================================
// replay.cpp — convert raw evio data to ROOT TTree for analysis
#include "replay.h"

using namespace fdec;

DaqMap replay::loadDaqMap(const std::string &json_path)
{
    DaqMap map;
    std::ifstream f(json_path);
    if (!f) {
        std::cerr << "Failed to open DAQ map file: " << json_path << "\n";
        return map;
    }
    auto json = nlohmann::json::parse(f);
    for(auto &entry : json){
        int crate = entry["crate"].get<int>();
        int slot = entry["slot"].get<int>();
        int channel = entry["channel"].get<int>();
        std::string name = entry["name"].get<std::string>();
        std::string key = std::to_string(crate) + "_" + std::to_string(slot) + "_" + std::to_string(channel);
        map[key] = name;
    }
    return map;
}

std::string replay::GetModuleName(int roc, int slot, int channel)
{
    std::string key = std::to_string(roc) + "_" + std::to_string(slot) + "_" + std::to_string(channel);
    auto it = daq_map.find(key);
    if (it != daq_map.end()) return it->second;
    else return "Unknown";
}

void replay::clear(){
    ftrigger = 0;
    ftimestamp = 0;
    Nchs = 0;
    memset(crate_id, 0, sizeof(crate_id));
    memset(slot_id, 0, sizeof(slot_id));
    memset(channel_id, 0, sizeof(channel_id));
    memset(module_id, 0, sizeof(module_id));
    memset(sample_num, 0, sizeof(sample_num));
    memset(samples, 0, sizeof(samples));
    memset(ped, 0, sizeof(ped));
    memset(ped_err, 0, sizeof(ped_err));
    memset(integral, 0, sizeof(integral));
    memset(Npeaks, 0, sizeof(Npeaks));
    memset(peak, 0, sizeof(peak));
    memset(peak_time, 0, sizeof(peak_time));
    memset(peak_inte, 0, sizeof(peak_inte));
}

void replay::setupTTree(TTree *tree, bool mode_peaks) {
    tree->Branch("event_num", &event_num, "event_num/I");
    tree->Branch("trigger", &ftrigger, "trigger/I");
    tree->Branch("timestamp", &ftimestamp, "timestamp/L");
    tree->Branch("Nchs", &Nchs, "Nchs/I");
    tree->Branch("crate_id", &crate_id[0], "crate_id[Nchs]/I");
    tree->Branch("slot_id", &slot_id[0], "slot_id[Nchs]/I");
    tree->Branch("channel_id", &channel_id[0], "channel_id[Nchs]/I");
    tree->Branch("module_id", &module_id[0], "module_id[Nchs]/I");
    if(!mode_peaks) {
        tree->Branch("sample_num", &sample_num[0], "sample_num[Nchs]/I");
        tree->Branch("samples", &samples[0][0], Form("samples[Nchs][%d]/I", fdec::MAX_SAMPLES));
    } 
    else {
        tree->Branch("ped", &ped[0], "ped[Nchs]/F");
        tree->Branch("ped_err", &ped_err[0], "ped_err[Nchs]/F");
        tree->Branch("integral", &integral[0], "integral[Nchs]/F");
        tree->Branch("Npeaks", &Npeaks[0], "Npeaks[Nchs]/I");
        tree->Branch("peak_height", &peak[0][0], Form("peak_height[Nchs][%d]/F", MAX_PEAKS));
        tree->Branch("peak_time", &peak_time[0][0], Form("peak_time[Nchs][%d]/F", MAX_PEAKS));
        tree->Branch("peak_integral", &peak_inte[0][0], Form("peak_integral[Nchs][%d]/F", MAX_PEAKS));
    }
}

void replay::setupTTreeRead(TTree *tree, bool mode_peaks) {
    tree->SetBranchAddress("event_num", &event_num);
    tree->SetBranchAddress("trigger", &ftrigger);
    tree->SetBranchAddress("timestamp", &ftimestamp);
    tree->SetBranchAddress("Nchs", &Nchs);
    tree->SetBranchAddress("crate_id", &crate_id[0]);
    tree->SetBranchAddress("slot_id", &slot_id[0]);
    tree->SetBranchAddress("channel_id", &channel_id[0]);
    tree->SetBranchAddress("module_id", &module_id[0]);
    if(!mode_peaks) {
        tree->SetBranchAddress("sample_num", &sample_num[0]);
        tree->SetBranchAddress("samples", &samples[0][0]);
    } 
    else {
        tree->SetBranchAddress("ped", &ped[0]);
        tree->SetBranchAddress("ped_err", &ped_err[0]);
        tree->SetBranchAddress("integral", &integral[0]);
        tree->SetBranchAddress("Npeaks", &Npeaks[0]);
        tree->SetBranchAddress("peak_height", &peak[0][0]);
        tree->SetBranchAddress("peak_time", &peak_time[0][0]);
        tree->SetBranchAddress("peak_integral", &peak_inte[0][0]);
    }
}

bool replay::rawdata2root(const std::string &input_path, const std::string &output_path, const int max_events, const bool mode_peaks)
{
    evc::EvChannel ch;
    if (ch.Open(input_path) != evc::status::success) {
        std::cerr << "Failed to open: " << input_path << "\n"; return false;
    }

    TFile *file_out = new TFile((output_path + ".root").c_str(), "RECREATE");
    TTree *tree = new TTree("T", "T");
    setupTTree(tree, mode_peaks);

    daq_map = loadDaqMap("/home/liyuan/PRad2_script/prad2evviewer/database/daq_map.json");

    EventData event;
    event_num = 0;

    while (ch.Read() == evc::status::success) {
        if (!ch.Scan()) continue;

        int nevt = ch.GetNEvents();
        for (int i = 0; i < nevt; ++i) {
            clear();

            // decode this sub-event
            if (ch.DecodeEvent(i, event)) {
                ftrigger = event.info.trigger_bits;
                ftimestamp = event.info.timestamp;
                for (int r = 0; r < event.nrocs; ++r) {
                    if (!event.rocs[r].present) continue;
                    auto &roc = event.rocs[r];

                    int roc_id = -1;
                    if(roc.tag == 0x80) roc_id = 0;
                    if(roc.tag == 0x82) roc_id = 1;
                    if(roc.tag == 0x84) roc_id = 2;
                    if(roc.tag == 0x86) roc_id = 3;
                    if(roc.tag == 0x88) roc_id = 4;
                    if(roc.tag == 0x8a) roc_id = 5;
                    if(roc.tag == 0x8c) roc_id = 6;

                    for (int s = 0; s < MAX_SLOTS; ++s) {
                        auto &slot = roc.slots[s];
                        if (!slot.present) continue;

                        for (int c = 0; c < FADC_CHANNELS; ++c){
                            if (!(slot.channel_mask & (1u << c)))
                                continue;
                            if(Nchs >= max_chs){
                                std::cerr << "Warning: reached max channels per event (" << max_chs << "), some data may be truncated.\n";
                                break;
                            }

                            auto &cd = slot.channels[c];

                            if(!mode_peaks){
                                sample_num[Nchs] = cd.nsamples;
                                for(int j = 0; j < cd.nsamples; j++){
                                    samples[Nchs][j] = cd.samples[j];
                                }
                            }else{
                                fdec::WaveAnalyzer ana;
                                fdec::WaveResult wave_res;
                                ana.Analyze(cd.samples, cd.nsamples, wave_res);

                                ped[Nchs] = wave_res.ped.mean;
                                ped_err[Nchs] = wave_res.ped.rms;
                                integral[Nchs] = GetIntegral(cd, ped[Nchs]);
                                Npeaks[Nchs] = wave_res.npeaks;
                                for (int p = 0; p < wave_res.npeaks; p++) {
                                    if (p >= MAX_PEAKS) break; // safety check for max peaks
                                    auto &pk = wave_res.peaks[p];
                                    peak[Nchs][p] = pk.height;
                                    peak_time[Nchs][p] = pk.time;
                                    peak_inte[Nchs][p] = pk.integral;
                                }
                            }
                            crate_id[Nchs] = roc_id;
                            slot_id[Nchs] = s;
                            channel_id[Nchs] = c;

                            std::string module_name = GetModuleName(roc_id, s, c);
                            if(module_name[0] == 'G') 
                                module_id[Nchs] = std::stoi(module_name.substr(1));
                            else if(module_name[0] == 'W') 
                                module_id[Nchs] = std::stoi(module_name.substr(1)) + 1000;
                            else 
                                module_id[Nchs] = -1; // Unknown module

                            Nchs++;
                        }
                    }
                }
                event_num++;
                tree->Fill();
                if (max_events > 0 && event_num >= max_events) goto done;
                if (event_num % PROGRESS_COUNT == 0) {
                        std::cout << "------[Decoded " << event_num << " events ]---"
                            << "\r" << std::flush;
                    }
            }
        }
    }
    done:
    ch.Close();
    file_out->cd();
    file_out->Write();
    file_out->Close();
    std::cout << "------[Decoded " << event_num << " events ]------\n";
    return true;
}



float replay::GetIntegral(const fdec::ChannelData &cd, const float pedestal) {
    float integral = 0;
    for (int i = 0; i < cd.nsamples; i++) {
        integral += float(cd.samples[i]) - pedestal;
    }
    return integral;
}