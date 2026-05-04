#pragma once
//=============================================================================
// TrackGeometry.h — straight-line tracking primitives
//
// Lab-frame two-point seed line + 4-parameter weighted-LSQ line fit + a
// helper to project the resulting line onto a planar detector's local
// frame.  No magnetic field, no curvature corrections — PRad-II tracks
// are straight lines because the spectrometer has no field.
//
// Lifted out of src/app_state.cpp so AppState, the analysis tools, and
// any future track matcher can share one implementation.
//=============================================================================

#include "DetectorTransform.h"

namespace prad2::trk {

// 4-parameter line: x(z) = ax + bx·z, y(z) = ay + by·z, plus the χ²/dof
// from the most recent fit (0 for seed lines that never went through
// fitWeightedLine).
struct Line3D {
    float ax = 0.f, bx = 0.f;
    float ay = 0.f, by = 0.f;
    float chi2_per_dof = 0.f;
};

// Two-point seed line in lab frame.  When |z2 - z1| is degenerate the
// returned line has bx = by = 0 and (ax, ay) = (x1, y1) — caller decides
// whether that's acceptable.
Line3D seedLine(float x1, float y1, float z1,
                float x2, float y2, float z2);

// Independent weighted LSQ fits in (z, x) and (z, y).  `wy = nullptr`
// reuses `wx` for both axes (the common case where σ_x = σ_y); pass
// distinct arrays for anisotropic per-point uncertainties (target point
// in target-in-fit mode, where σ_target_z couples differently into σ_x
// and σ_y via the slope).  dof = 2N - 4; for N == 2 χ²/dof is set to 0.
//
// Returns false if N < 2 or the normal-equations determinant is
// degenerate; `out` is unspecified in that case.
bool fitWeightedLine(int N,
                     const float *z, const float *x, const float *y,
                     const float *wx, const float *wy,
                     Line3D &out);

inline bool fitWeightedLine(int N,
                            const float *z, const float *x, const float *y,
                            const float *w, Line3D &out)
{
    return fitWeightedLine(N, z, x, y, w, nullptr, out);
}

// Project a lab-frame line onto a detector's local plane (z_local = 0)
// using the labToLocal-then-1D-interpolate trick.  Returns the local
// (x, y) where the line crosses z_local = 0.
void projectLineToLocal(const DetectorTransform &xform, const Line3D &L,
                        float &px, float &py);

} // namespace prad2::trk
