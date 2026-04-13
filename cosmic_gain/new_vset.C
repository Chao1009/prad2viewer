void new_vset(){

    const int NMODULES = 1156;

    const int LG_num = 76;
    int LG_module_id[LG_num] = 
    {156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174,
    186, 216, 246, 276, 306, 336, 366, 396, 426, 456, 486, 516, 546, 576, 606, 636, 666, 696, 726,
    175, 205, 235, 265, 295, 325, 355, 385, 415, 445, 475, 505, 535, 565, 595, 625, 655, 685, 715,
    727, 728, 729, 730, 731, 732, 733, 734, 735, 736, 737, 738, 739, 740, 741, 742, 743, 744, 745
    };

    // vset3 from vset_iter3_new.json, mean from cosmic_modules_run561.json
    double vset[NMODULES+LG_num] = {};
    double mean[NMODULES+LG_num]  = {};
    bool   has_vset[NMODULES+LG_num] = {};
    bool   has_run[NMODULES+LG_num] = {};

    // Helper: find index for G module id, returns NMODULES+j or -1
    auto find_lg_idx = [&](int gid) -> int {
        for (int j = 0; j < LG_num; j++) {
            if (LG_module_id[j] == gid) return NMODULES + j;
        }
        return -1;
    };

    std::string vset_file = "vset_iter1.json";
    std::string cosmic_file = "cosmic_modules_run575.json";
    std::string output_file = "vset_iter2.json";

    // ---- Read vset_iter3_new.json: V0Set per channel ----
    {
        std::ifstream fin(vset_file);
        if (!fin.is_open()) { std::cerr << "Cannot open " << vset_file << std::endl; return; }
        std::string line;
        int current_mod = -1;
        bool current_is_g = false;
        while (std::getline(fin, line)) {
            int mod = 0;
            if (line.find("\"name\"") != std::string::npos) {
                if (sscanf(line.c_str(), " \"name\": \"W%d\"", &mod) == 1) {
                    current_mod = mod;
                    current_is_g = false;
                } else if (sscanf(line.c_str(), " \"name\": \"G%d\"", &mod) == 1) {
                    current_mod = mod;
                    current_is_g = true;
                }
            }
            int idx = -1;
            if (!current_is_g && current_mod > 0 && current_mod <= NMODULES) {
                idx = current_mod - 1;
            } else if (current_is_g) {
                idx = find_lg_idx(current_mod);
            }
            if (idx >= 0) {
                double val;
                if (line.find("\"V0Set\"") != std::string::npos) {
                    if (sscanf(line.c_str(), " \"V0Set\": %lf", &val) == 1) {
                        vset[idx] = val;
                        has_vset[idx] = true;
                    }
                }
            }
        }
    }

    // ---- Read cosmic_modules_run561.json: peak_height_mean ----
    {
        std::ifstream fin(cosmic_file);
        if (!fin.is_open()) { std::cerr << "Cannot open " << cosmic_file << std::endl; return; }
        std::string line;
        while (std::getline(fin, line)) {
            int mod = 0;
            int run_dummy = 0;
            double v_ph_mean = 0, v_ph_sigma = 0, v_ph_diff = 0;
            double v_pi_mean = 0, v_pi_sigma = 0, v_pi_diff = 0;
            int cnt = 0;
            if (sscanf(line.c_str(),
                       " \"W%d\": [{\"run\": %d, \"peak_height_mean\": %lf, \"peak_height_sigma\": %lf, \"peak_height_diff\": %lf,"
                       " \"peak_integral_mean\": %lf, \"peak_integral_sigma\": %lf, \"peak_integral_diff\": %lf,"
                       " \"count\": %d",
                       &mod, &run_dummy, &v_ph_mean, &v_ph_sigma, &v_ph_diff,
                       &v_pi_mean, &v_pi_sigma, &v_pi_diff, &cnt) == 9) {
                int idx = mod - 1;
                if (idx >= 0 && idx < NMODULES) {
                    mean[idx] = v_ph_mean;
                    has_run[idx] = true;
                }
            } else if (sscanf(line.c_str(),
                       " \"G%d\": [{\"run\": %d, \"peak_height_mean\": %lf, \"peak_height_sigma\": %lf, \"peak_height_diff\": %lf,"
                       " \"peak_integral_mean\": %lf, \"peak_integral_sigma\": %lf, \"peak_integral_diff\": %lf,"
                       " \"count\": %d",
                       &mod, &run_dummy, &v_ph_mean, &v_ph_sigma, &v_ph_diff,
                       &v_pi_mean, &v_pi_sigma, &v_pi_diff, &cnt) == 9) {
                int idx = find_lg_idx(mod);
                if (idx >= 0) {
                    mean[idx] = v_ph_mean;
                    has_run[idx] = true;
                }
            }
        }
    }

    // ---- Compute vset4 ----
    double vsetnew[NMODULES+LG_num] = {};
    bool   valid[NMODULES+LG_num] = {};

    int n_valid = 0, n_skip = 0;
    int n_increase[3] = {0}, n_decrease[3] = {0}, n_unchanged = 0;

    for (int i = 0; i < NMODULES+LG_num; i++) {
        if (!has_vset[i] || !has_run[i] || i == 1018 || i == 1019) { n_skip++; continue; }

        vsetnew[i] = vset[i];

        if (mean[i] > 45.0) {
            vsetnew[i] = vset[i] - 20.0;
            n_decrease[2]++;
        }else if (mean[i] > 40.0) {
            vsetnew[i] = vset[i] - 10.0;
            n_decrease[1]++;
        }else if (mean[i] > 37.0) {
            vsetnew[i] = vset[i] - 5.0;
            n_decrease[0]++;
        }

        if(mean[i] < 25.0) {
            vsetnew[i] = vset[i] + 20.0;
            n_increase[2]++;
        }else if(mean[i] < 30.0) {
            vsetnew[i] = vset[i] + 10.0;
            n_increase[1]++;
        }else if(mean[i] < 33.0) {
            vsetnew[i] = vset[i] + 5.0;
            n_increase[0]++;
        }
        
        if(mean[i] >= 33.0 && mean[i] <= 37.0) {
            vsetnew[i] = vset[i];
            n_unchanged++;
        }
        if (vsetnew[i] > 1270.0 && i < NMODULES) vsetnew[i] = 1270.0;
        if (vsetnew[i] > 1800 && i >= NMODULES) vsetnew[i] = 1800.0;
        valid[i] = true;
        n_valid++;
    }

    printf("\n=== Summary: %d valid, %d skipped ===\n", n_valid, n_skip);
    printf("Voltage increased (+20V, height<25):  %d\n", n_increase[2]);
    printf("Voltage increased (+10V, height<30):  %d\n", n_increase[1]);
    printf("Voltage increased (+5V, height<33):  %d\n", n_increase[0]);
    printf("Voltage decreased (-20V, height>45):  %d\n", n_decrease[2]);
    printf("Voltage decreased (-10V, height>40):  %d\n", n_decrease[1]);
    printf("Voltage decreased (-5V, height>37):  %d\n", n_decrease[0]);
    printf("Voltage unchanged (33~37):            %d\n\n", n_unchanged);

    printf("%-6s %10s %10s %10s %8s\n", "Ch", "vset", "height", "vsetnew", "action");
    printf("----------------------------------------------------------\n");
    for (int i = 0; i < NMODULES + LG_num; i++) {
        if (!valid[i]) continue;
        const char* action = "keep";
        double dv = vsetnew[i] - vset[i];
        if (dv > 15.0) action = "+20V";
        else if (dv < -15.0) action = "-20V";
        else if (dv > 7.0) action = "+10V";
        else if (dv < -7.0) action = "-10V";
        else if (dv > 2.0) action = "+5V";
        else if (dv < -2.0) action = "-5V";
        if (i < NMODULES)
            printf("W%-5d %10.1f %10.2f %10.1f %8s\n",
                   i+1, vset[i], mean[i], vsetnew[i], action);
        else
            printf("G%-5d %10.1f %10.2f %10.1f %8s\n",
                   LG_module_id[i - NMODULES], vset[i], mean[i], vsetnew[i], action);
    }

    // ---- Write vset_iter4.json: based on vset_iter3_new.json, V0Set replaced by vsetnew ----
    {
        std::ifstream fin(vset_file);
        std::ofstream fout(output_file);
        if (!fin.is_open() || !fout.is_open()) {
            std::cerr << "Cannot open files for writing" << std::endl;
            return;
        }
        std::string line;
        int current_mod = -1;
        bool current_is_g = false;
        while (std::getline(fin, line)) {
            int mod = 0;
            if (line.find("\"name\"") != std::string::npos) {
                if (sscanf(line.c_str(), " \"name\": \"W%d\"", &mod) == 1) {
                    current_mod = mod;
                    current_is_g = false;
                } else if (sscanf(line.c_str(), " \"name\": \"G%d\"", &mod) == 1) {
                    current_mod = mod;
                    current_is_g = true;
                }
            }
            int idx = -1;
            if (!current_is_g && current_mod > 0 && current_mod <= NMODULES) {
                idx = current_mod - 1;
            } else if (current_is_g) {
                idx = find_lg_idx(current_mod);
            }
            if (line.find("\"V0Set\"") != std::string::npos && idx >= 0) {
                if (valid[idx]) {
                    bool has_comma = (line.find(",") != std::string::npos &&
                                      line.rfind(",") > line.find("V0Set"));
                    char buf[256];
                    if (has_comma)
                        snprintf(buf, sizeof(buf), "        \"V0Set\": %.1f,", vsetnew[idx]);
                    else
                        snprintf(buf, sizeof(buf), "        \"V0Set\": %.1f", vsetnew[idx]);
                    fout << buf << "\n";
                    continue;
                }
            }
            fout << line << "\n";
        }
        printf("\n%s written.\n", output_file.c_str());
    }
}
