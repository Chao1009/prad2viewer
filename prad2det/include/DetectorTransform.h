#pragma once
//=============================================================================
// DetectorTransform.h — planar detector coordinate transform
//
// Transforms detector-plane coordinates (x, y) to lab frame via
// Euler rotation (Rx * Ry * Rz) then translation.
// Reusable for HyCal, GEMs, or any planar detector.
//
// The rotation matrix is lazily computed on first use and cached.
// Call prepare() explicitly to force precomputation, or just use
// toLab()/rotate() — they auto-prepare if needed.
//=============================================================================

#include <cmath>

struct DetectorTransform {
    float x=0, y=0, z=0;               // detector origin in lab frame (mm)
    float rx=0, ry=0, rz=0;            // tilting angles (degrees)

    // Precomputed rotation matrix elements (full 3x3, row-major: rIJ).
    struct Matrix {
        float r00=1, r01=0, r02=0;
        float r10=0, r11=1, r12=0;
        float r20=0, r21=0, r22=1;
        float tx=0, ty=0, tz=0;
    };

    // Mark the cached rotation matrix as stale so the next prepare() /
    // toLab() / rotate() / matrix() call rebuilds it.  Use this after
    // any direct mutation of x/y/z/rx/ry/rz — those raw field writes
    // don't know about the cache.  set() and the Python property setters
    // call this for you.
    void invalidate() { prepared_ = false; }

    // Set pose (translation in mm + tilts in degrees) and rebuild the
    // cached rotation matrix in one shot.  Equivalent to writing the six
    // fields plus invalidate() + prepare(); preferred over field-at-a-
    // time mutation because it can't leave the cache half-built.
    void set(float x_, float y_, float z_,
             float rx_, float ry_, float rz_) {
        x = x_; y = y_; z = z_;
        rx = rx_; ry = ry_; rz = rz_;
        invalidate();
        prepare();
    }

    // Force matrix precomputation (idempotent — first call only).  Call
    // invalidate() first if you mutated any field after the previous
    // prepare(), or just use set() to do both at once.
    void prepare() const {
        if (prepared_) return;
        const float DEG = 3.14159265f / 180.f;
        float cx=std::cos(rx*DEG), sx=std::sin(rx*DEG);
        float cy=std::cos(ry*DEG), sy=std::sin(ry*DEG);
        float cz=std::cos(rz*DEG), sz=std::sin(rz*DEG);
        // R = Rx * Ry * Rz (Euler XYZ, intrinsic).
        mat_.r00 =  cy*cz;              mat_.r01 = -cy*sz;              mat_.r02 =  sy;
        mat_.r10 =  sx*sy*cz + cx*sz;   mat_.r11 = -sx*sy*sz + cx*cz;   mat_.r12 = -sx*cy;
        mat_.r20 = -cx*sy*cz + sx*sz;   mat_.r21 =  cx*sy*sz + sx*cz;   mat_.r22 =  cx*cy;
        mat_.tx = x;  mat_.ty = y;  mat_.tz = z;
        prepared_ = true;
    }

    // Transform a point from the detector plane (z_local = 0) to lab frame.
    void toLab(float dx, float dy, float &lx, float &ly, float &lz) const {
        prepare();
        lx = mat_.r00*dx + mat_.r01*dy + mat_.tx;
        ly = mat_.r10*dx + mat_.r11*dy + mat_.ty;
        lz = mat_.r20*dx + mat_.r21*dy + mat_.tz;
    }

    // Transform a 3D local point to lab frame (e.g. HyCal cluster with
    // shower depth as z_local).  Equivalent to lab = R * [dx,dy,dz] + [tx,ty,tz].
    void toLab(float dx, float dy, float dz,
               float &lx, float &ly, float &lz) const {
        prepare();
        lx = mat_.r00*dx + mat_.r01*dy + mat_.r02*dz + mat_.tx;
        ly = mat_.r10*dx + mat_.r11*dy + mat_.r12*dz + mat_.ty;
        lz = mat_.r20*dx + mat_.r21*dy + mat_.r22*dz + mat_.tz;
    }

    // Inverse: lab → detector-local. R is orthonormal so R^{-1} = R^T.
    void labToLocal(float lx, float ly, float lz,
                    float &dx, float &dy, float &dz) const {
        prepare();
        float ux = lx - mat_.tx, uy = ly - mat_.ty, uz = lz - mat_.tz;
        dx = mat_.r00*ux + mat_.r10*uy + mat_.r20*uz;
        dy = mat_.r01*ux + mat_.r11*uy + mat_.r21*uz;
        dz = mat_.r02*ux + mat_.r12*uy + mat_.r22*uz;
    }

    // Rotation only (no translation). For drawing in detector-local space.
    void rotate(float dx, float dy, float &ox, float &oy) const {
        prepare();
        ox = mat_.r00*dx + mat_.r01*dy;
        oy = mat_.r10*dx + mat_.r11*dy;
    }

    // The detector's surface normal in the lab frame (third column of R).
    void normal(float &nx, float &ny, float &nz) const {
        prepare();
        nx = mat_.r02;  ny = mat_.r12;  nz = mat_.r22;
    }

    // Access the cached matrix directly (auto-prepares).
    const Matrix& matrix() const { prepare(); return mat_; }

private:
    mutable Matrix mat_;
    mutable bool   prepared_ = false;
};
