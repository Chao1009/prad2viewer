#pragma once
//=============================================================================
// viewer_utils.h — shared utilities for evc_viewer and evc_monitor
//=============================================================================

#include "Fadc250Data.h"
#include "WaveAnalyzer.h"

#include <fstream>
#include <string>
#include <vector>
#include <cmath>

// --- Detector coordinate transform ------------------------------------------
// Reusable for HyCal, GEMs, or any planar detector.
// Transforms detector-plane coordinates (x, y, 0) to lab frame.
// Convention: rotation applied first (Rx * Ry * Rz), then translation.
struct DetectorTransform {
    float x=0, y=0, z=0;               // detector origin in lab frame (mm)
    float rx=0, ry=0, rz=0;            // tilting angles (degrees)

    // Transform a point from detector plane to lab frame.
    void toLab(float dx, float dy, float &lx, float &ly, float &lz) const {
        // convert angles to radians
        const float DEG = 3.14159265f / 180.f;
        float cx=std::cos(rx*DEG), sx=std::sin(rx*DEG);
        float cy=std::cos(ry*DEG), sy=std::sin(ry*DEG);
        float cz=std::cos(rz*DEG), sz=std::sin(rz*DEG);
        // R = Rx * Ry * Rz applied to (dx, dy, 0)
        float px =  cy*cz*dx - cy*sz*dy;
        float py = (sx*sy*cz + cx*sz)*dx + (-sx*sy*sz + cx*cz)*dy;
        float pz = (-cx*sy*cz + sx*sz)*dx + (cx*sy*sz + sx*cz)*dy;
        lx = px + x;
        ly = py + y;
        lz = pz + z;
    }
};

// --- TI timestamp conversion ------------------------------------------------
// TI clock runs at 250 MHz → 4 ns per tick
static constexpr double TI_TICK_SEC = 4e-9;

// --- file I/O helpers -------------------------------------------------------
inline std::string readFile(const std::string &path) {
    std::ifstream f(path);
    if (!f) return "";
    return {std::istreambuf_iterator<char>(f), {}};
}

inline std::string findFile(const std::string &name, const std::string &base) {
    { std::ifstream f(name); if (f.good()) return name; }
    std::string p = base + "/" + name;
    { std::ifstream f(p); if (f.good()) return p; }
    return "";
}

inline std::string contentType(const std::string &path) {
    if (path.size() >= 5 && path.substr(path.size()-5) == ".html") return "text/html; charset=utf-8";
    if (path.size() >= 4 && path.substr(path.size()-4) == ".css")  return "text/css; charset=utf-8";
    if (path.size() >= 3 && path.substr(path.size()-3) == ".js")   return "application/javascript; charset=utf-8";
    return "application/octet-stream";
}

// --- Histogram (used by both viewer and monitor) ----------------------------
struct Histogram {
    int underflow = 0, overflow = 0;
    std::vector<int> bins;
    void init(int n) { bins.assign(n, 0); underflow = overflow = 0; }
    void fill(float v, float bmin, float bstep) {
        if (v < bmin) { ++underflow; return; }
        int b = (int)((v - bmin) / bstep);
        if (b >= (int)bins.size()) { ++overflow; return; }
        ++bins[b];
    }
    void clear() { std::fill(bins.begin(), bins.end(), 0); underflow = overflow = 0; }
};

struct Histogram2D {
    int nx = 0, ny = 0;
    std::vector<int> bins;  // row-major: bins[iy*nx + ix]
    void init(int nx_, int ny_) { nx = nx_; ny = ny_; bins.assign(nx * ny, 0); }
    void fill(float vx, float vy, float xmin, float xstep, float ymin, float ystep) {
        int ix = (int)((vx - xmin) / xstep);
        int iy = (int)((vy - ymin) / ystep);
        if (ix < 0 || ix >= nx || iy < 0 || iy >= ny) return;
        bins[iy * nx + ix]++;
    }
    void clear() { std::fill(bins.begin(), bins.end(), 0); }
};

// --- Histogram config -------------------------------------------------------
struct HistConfig {
    float time_min  = 170;
    float time_max  = 190;
    float bin_min   = 0;
    float bin_max   = 20000;
    float bin_step  = 100;
    float threshold = 3.0;
    float pos_min   = 0;
    float pos_max   = 400;
    float pos_step  = 4;
    float min_peak_ratio = 0.3f;
};

// --- LMS entry (shared between viewer FileData and monitor globals) ---------
struct LmsEntry {
    double time_sec;    // seconds since first LMS event (from TI timestamp)
    float  integral;    // peak integral within timing cut (or raw ADC for ADC1881M)
};

// --- Peak extraction helper -------------------------------------------------
// Find best peak integral within time window. Returns -1 if no peak found.
inline float bestPeakInWindow(const fdec::WaveResult &wres,
                               float threshold, float time_min, float time_max)
{
    float best = -1;
    for (int p = 0; p < wres.npeaks; ++p) {
        auto &pk = wres.peaks[p];
        if (pk.height < threshold) continue;
        if (pk.time >= time_min && pk.time <= time_max)
            if (pk.integral > best) best = pk.integral;
    }
    return best;
}
