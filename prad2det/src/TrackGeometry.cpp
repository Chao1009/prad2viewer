// TrackGeometry.cpp — straight-line tracking primitives
//=============================================================================
// Direct port of the anonymous-namespace helpers that used to live at
// the top of src/app_state.cpp (Line3D / seedLine / fitWeightedLine /
// projectLineToLocal).  Kept here so prad2det owns the math and any
// caller — viewer, replay tools, pybind users — gets the same numbers.
//=============================================================================

#include "TrackGeometry.h"

#include <cmath>

namespace prad2::trk {

Line3D seedLine(float x1, float y1, float z1,
                float x2, float y2, float z2)
{
    Line3D L{};
    float dz = z2 - z1;
    if (std::abs(dz) < 1e-6f) { L.ax = x1; L.ay = y1; return L; }
    L.bx = (x2 - x1) / dz;  L.ax = x1 - L.bx * z1;
    L.by = (y2 - y1) / dz;  L.ay = y1 - L.by * z1;
    return L;
}

bool fitWeightedLine(int N,
                     const float *z, const float *x, const float *y,
                     const float *wx, const float *wy,
                     Line3D &out)
{
    if (N < 2) return false;
    if (wy == nullptr) wy = wx;

    double Swx=0, Szx=0, Szzx=0, Sx=0, Sxz=0;
    for (int i = 0; i < N; ++i) {
        double wi = wx[i];
        Swx  += wi;
        Szx  += wi * z[i];
        Szzx += wi * z[i] * z[i];
        Sx   += wi * x[i];
        Sxz  += wi * x[i] * z[i];
    }
    double Dx = Swx * Szzx - Szx * Szx;
    if (std::abs(Dx) < 1e-9) return false;
    double bx = (Swx * Sxz - Szx * Sx) / Dx;
    double ax = (Sx - bx * Szx) / Swx;

    double Swy=0, Szy=0, Szzy=0, Sy=0, Syz=0;
    for (int i = 0; i < N; ++i) {
        double wi = wy[i];
        Swy  += wi;
        Szy  += wi * z[i];
        Szzy += wi * z[i] * z[i];
        Sy   += wi * y[i];
        Syz  += wi * y[i] * z[i];
    }
    double Dy = Swy * Szzy - Szy * Szy;
    if (std::abs(Dy) < 1e-9) return false;
    double by = (Swy * Syz - Szy * Sy) / Dy;
    double ay = (Sy - by * Szy) / Swy;

    out.ax = (float)ax; out.bx = (float)bx;
    out.ay = (float)ay; out.by = (float)by;
    int dof = 2 * N - 4;
    if (dof > 0) {
        double chi2 = 0;
        for (int i = 0; i < N; ++i) {
            double dxp = (ax + bx * z[i]) - x[i];
            double dyp = (ay + by * z[i]) - y[i];
            chi2 += wx[i] * dxp * dxp + wy[i] * dyp * dyp;
        }
        out.chi2_per_dof = (float)(chi2 / dof);
    } else {
        out.chi2_per_dof = 0.f;
    }
    return true;
}

void projectLineToLocal(const DetectorTransform &xform, const Line3D &L,
                        float &px, float &py)
{
    float ax1 = L.ax,                 ay1 = L.ay,                 z1 = 0.f;
    float ax2 = L.ax + L.bx * 1000.f, ay2 = L.ay + L.by * 1000.f, z2 = 1000.f;
    float l1x, l1y, l1z, l2x, l2y, l2z;
    xform.labToLocal(ax1, ay1, z1, l1x, l1y, l1z);
    xform.labToLocal(ax2, ay2, z2, l2x, l2y, l2z);
    float dz = l2z - l1z;
    if (std::abs(dz) < 1e-6f) { px = l1x; py = l1y; return; }
    float s = -l1z / dz;
    px = l1x + s * (l2x - l1x);
    py = l1y + s * (l2y - l1y);
}

} // namespace prad2::trk
