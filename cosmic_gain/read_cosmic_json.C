

void read_cosmic_json(){

    // ---- Load cosmic_modules_run555.json ----
    // New format: "WX": [{"peak_height_mean":V,"peak_height_sigma":V,"peak_height_diff":V,
    //                      "peak_integral_mean":V,"peak_integral_sigma":V,"peak_integral_diff":V,"count":N}]
    // Arrays indexed by module number - 1: W1 -> [0], W1156 -> [1155]
    const int NMODULES = 1156;
    const int LG_num = 76;
    int LG_module_id[LG_num] = 
    {156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174,
    186, 216, 246, 276, 306, 336, 366, 396, 426, 456, 486, 516, 546, 576, 606, 636, 666, 696, 726,
    175, 205, 235, 265, 295, 325, 355, 385, 415, 445, 475, 505, 535, 565, 595, 625, 655, 685, 715,
    727, 728, 729, 730, 731, 732, 733, 734, 735, 736, 737, 738, 739, 740, 741, 742, 743, 744, 745
    };

    auto find_lg_idx = [&](int gid) -> int {
        for (int j = 0; j < LG_num; j++) {
            if (LG_module_id[j] == gid) return NMODULES + j;
        }
        return -1;
    };

    const int NTOTAL = NMODULES + LG_num;
    double ph_mean[NTOTAL]   = {};   // peak_height_mean
    double ph_sigma[NTOTAL]  = {};   // peak_height_sigma
    double ph_diff[NTOTAL]   = {};   // peak_height_diff
    double pi_mean[NTOTAL]   = {};   // peak_integral_mean
    double pi_sigma[NTOTAL]  = {};   // peak_integral_sigma
    double pi_diff[NTOTAL]   = {};   // peak_integral_diff
    int    count[NTOTAL]     = {};   // count

    std::string json_file = "cosmic_modules_run570.json";
    std::ifstream fin(json_file);
    if (!fin.is_open()) {
        std::cerr << "Cannot open json file" << std::endl;
        return;
    }
    std::string line;
    while (std::getline(fin, line)) {
        int    mod = 0;
        double v_ph_mean = 0, v_ph_sigma = 0, v_ph_diff = 0;
        double v_pi_mean = 0, v_pi_sigma = 0, v_pi_diff = 0;
        int    cnt = 0;
        // parse: "WX": [{"run": N, "peak_height_mean": V, ...}]
        int run_dummy = 0;
        if (sscanf(line.c_str(),
                    " \"W%d\": [{\"run\": %d, \"peak_height_mean\": %lf, \"peak_height_sigma\": %lf, \"peak_height_diff\": %lf,"
                    " \"peak_integral_mean\": %lf, \"peak_integral_sigma\": %lf, \"peak_integral_diff\": %lf,"
                    " \"count\": %d",
                    &mod, &run_dummy, &v_ph_mean, &v_ph_sigma, &v_ph_diff,
                    &v_pi_mean, &v_pi_sigma, &v_pi_diff, &cnt) == 9) {
            int idx = mod - 1;  // W1 -> 0, W1156 -> 1155
            if (idx >= 0 && idx < NMODULES) {
                ph_mean[idx]  = v_ph_mean;
                ph_sigma[idx] = v_ph_sigma;
                ph_diff[idx]  = v_ph_diff;
                pi_mean[idx]  = v_pi_mean;
                pi_sigma[idx] = v_pi_sigma;
                pi_diff[idx]  = v_pi_diff;
                count[idx]    = cnt;
            }
        } else if (sscanf(line.c_str(),
                    " \"G%d\": [{\"run\": %d, \"peak_height_mean\": %lf, \"peak_height_sigma\": %lf, \"peak_height_diff\": %lf,"
                    " \"peak_integral_mean\": %lf, \"peak_integral_sigma\": %lf, \"peak_integral_diff\": %lf,"
                    " \"count\": %d",
                    &mod, &run_dummy, &v_ph_mean, &v_ph_sigma, &v_ph_diff,
                    &v_pi_mean, &v_pi_sigma, &v_pi_diff, &cnt) == 9) {
            int idx = find_lg_idx(mod);
            if (idx >= 0) {
                ph_mean[idx]  = v_ph_mean;
                ph_sigma[idx] = v_ph_sigma;
                ph_diff[idx]  = v_ph_diff;
                pi_mean[idx]  = v_pi_mean;
                pi_sigma[idx] = v_pi_sigma;
                pi_diff[idx]  = v_pi_diff;
                count[idx]    = cnt;
            }
        }
    }

    TH1D *h_ph_mean = new TH1D("h_ph_mean", ("Peak Height Mean per Module "+json_file+"; Peak Height ADC; Counts").c_str(), 50, 0, 100);
    TH1D *h_pi_mean = new TH1D("h_pi_mean", ("Peak Integral Mean per Module "+json_file+"; Integral ADC; Counts").c_str(), 50, 0, 500);
    for (int i = 0; i < NTOTAL; i++) {
        if (count[i] > 0) {
            h_ph_mean->Fill(ph_mean[i]);
            h_pi_mean->Fill(pi_mean[i]);
        }
    }

    TCanvas *c_peak = new TCanvas("c_peak", "Cosmic Peak per Module", 1400, 500);
    c_peak->Divide(2, 1);

    c_peak->cd(1);
    h_ph_mean->SetLineColor(kBlue);
    h_ph_mean->SetLineWidth(2);
    h_ph_mean->Draw();
    TLine *line_ph = new TLine(35, 0, 35, h_ph_mean->GetMaximum());
    line_ph->SetLineColor(kRed);
    line_ph->SetLineWidth(2);
    line_ph->SetLineStyle(2);
    line_ph->Draw();

    c_peak->cd(2);
    h_pi_mean->SetLineColor(kRed);
    h_pi_mean->SetLineWidth(2);
    h_pi_mean->Draw();
    TLine *line_pi = new TLine(250, 0, 250, h_pi_mean->GetMaximum());
    line_pi->SetLineColor(kBlue);
    line_pi->SetLineWidth(2);
    line_pi->SetLineStyle(2);
    line_pi->Draw();

}