#pragma once
//=============================================================================
// HyCalSystem.h — HyCal detector geometry, module lookup, and neighbor system
//
// Initialized once from JSON config files. Immutable after Init().
// Provides O(1) module lookup by name, ID, or DAQ address, plus pre-computed
// neighbor lists with quantized distances for clustering.
//
// No per-event state lives here — callers provide flat energy arrays.
//=============================================================================

#include <cstdint>
#include <string>
#include <vector>
#include <unordered_map>
#include <array>
#include <cmath>

namespace fdec
{

// --- capacity limits --------------------------------------------------------
static constexpr int HYCAL_MAX_MODULES  = 1800;  // 1728 real + margin
static constexpr int MAX_NEIGHBORS      = 12;    // max neighbors per module

// --- module types -----------------------------------------------------------
enum class ModuleType : int8_t { PbGlass = 0, PbWO4 = 1, LMS = 2, Unknown = -1 };

// crystal module IDs start from 1000 (PrimEx convention)
static constexpr int PWO_ID0 = 1000;

// --- layout flags (bit positions) -------------------------------------------
enum LayoutFlag : uint32_t {
    kPbGlass     = 0,
    kPbWO4       = 1,
    kTransition  = 2,   // boundary between PbWO4 and PbGlass
    kInnerBound  = 3,   // inner boundary of HyCal (beam hole)
    kOuterBound  = 4,   // outer boundary of HyCal
    kDeadModule  = 5,
    kDeadNeighbor = 6,
    // cluster flags (used on ModuleCluster::flag, not module layout)
    kSplit       = 7,   // cluster was split from a multi-maximum group
    kLeakCorr    = 8,   // leakage correction applied
};

inline void     set_bit(uint32_t &f, uint32_t b) { f |= (1u << b); }
inline void     clear_bit(uint32_t &f, uint32_t b) { f &= ~(1u << b); }
inline bool     test_bit(uint32_t f, uint32_t b) { return (f >> b) & 1u; }

// --- sector types -----------------------------------------------------------
enum class Sector : int { Center = 0, Top = 1, Right = 2, Bottom = 3, Left = 4, Max = 5 };

struct SectorInfo {
    int         id      = -1;
    ModuleType  mtype   = ModuleType::Unknown;
    double      msize_x = 0.;
    double      msize_y = 0.;
    // boundary rectangle: [x_min, y_min] to [x_max, y_max]
    // stored as 4 corner points (BL, TL, TR, BR) for intersection tests
    struct Point { double x, y; };
    std::array<Point, 4> boundpts{};

    void set_boundary(double x1, double y1, double x2, double y2)
    {
        // counter-clockwise: BL, TL, TR, BR  (matches old code: [3]=TL,[0]=BL,[1]=TR,[2]=BR)
        // old code order: (x1,y2), (x1,y1), (x2,y1), (x2,y2)
        boundpts[0] = {x1, y2};  // top-left
        boundpts[1] = {x1, y1};  // bottom-left
        boundpts[2] = {x2, y1};  // bottom-right
        boundpts[3] = {x2, y2};  // top-right
    }

    void get_boundary(double &x1, double &y1, double &x2, double &y2) const
    {
        x1 = boundpts[1].x;  y1 = boundpts[1].y;
        x2 = boundpts[3].x;  y2 = boundpts[3].y;
    }
};

// --- neighbor info ----------------------------------------------------------
struct NeighborInfo {
    int   index;        // index into Module array
    float dx, dy;       // quantized distances
    float dist;         // sqrt(dx*dx + dy*dy)
};

// --- DAQ address ------------------------------------------------------------
struct DaqAddr {
    int crate   = -1;
    int slot    = -1;
    int channel = -1;
};

// --- module -----------------------------------------------------------------
struct Module {
    // identity
    std::string name;
    int         id          = -1;       // PrimEx ID (G: 1-576, W: 1001-2152)
    int         index       = -1;       // contiguous array index

    // geometry
    ModuleType  type        = ModuleType::Unknown;
    double      x           = 0.;
    double      y           = 0.;
    double      size_x      = 0.;
    double      size_y      = 0.;

    // layout
    uint32_t    flag        = 0;
    int         sector      = -1;       // Sector enum as int
    int         row         = 0;
    int         column      = 0;

    // DAQ mapping
    DaqAddr     daq;

    // calibration (per-module)
    double      cal_factor      = 0.;   // MeV per ADC count (or per integral unit)
    double      cal_base_energy = 0.;   // calibration beam energy (MeV)
    double      cal_non_linear  = 0.;   // non-linear correction coefficient

    // convert ADC value (pedestal-subtracted) to energy in MeV
    // includes non-linear correction: E + nl * (E - E_cal) / 1000
    double energize(double adc) const {
        if (adc < 0.) return 0.;
        double E = cal_factor * adc;
        return E + cal_non_linear * (E - cal_base_energy) / 1000.;
    }

    // pre-computed neighbors (filled by InitLayout)
    int         neighbor_count = 0;
    NeighborInfo neighbors[MAX_NEIGHBORS];

    // helpers
    bool is_pwo4()   const { return type == ModuleType::PbWO4; }
    bool is_glass()  const { return type == ModuleType::PbGlass; }
    bool is_hycal()  const { return is_pwo4() || is_glass(); }

    bool is_neighbor(int other_index, bool include_corners = true) const
    {
        for (int i = 0; i < neighbor_count; ++i) {
            if (neighbors[i].index == other_index) {
                if (include_corners)
                    return (std::abs(neighbors[i].dx) < 1.01f &&
                            std::abs(neighbors[i].dy) < 1.01f);
                else
                    return neighbors[i].dist < 1.2f;
            }
        }
        return false;
    }
};

// --- HyCalSystem ------------------------------------------------------------
class HyCalSystem
{
public:
    HyCalSystem() = default;

    // Initialize from JSON config files.
    // modules_path: hycal_modules.json (geometry)
    // daq_path:     daq_map.json (crate/slot/channel mapping)
    // Returns false on error.
    bool Init(const std::string &modules_path, const std::string &daq_path);

    // Load per-module calibration from JSON file.
    // Format: [{"name":"W735","factor":0.37,"base_energy":2138.67,"non_linear":0.006}, ...]
    // Returns number of modules matched, or -1 on error.
    int LoadCalibration(const std::string &calib_path);
    double GetCalibConstant(const int primex_id) const;
    void SetCalibConstant(const int primex_id, const double factor);

    // --- module access (all O(1)) -------------------------------------------
    int             module_count()                     const { return n_modules_; }
    const Module   &module(int index)                  const { return modules_[index]; }
    Module         &module(int index)                        { return modules_[index]; }
    const Module   *module_by_name(const std::string &name) const;
    const Module   *module_by_id(int primex_id)             const;
    const Module   *module_by_daq(int crate, int slot, int ch) const;

    // --- sector info --------------------------------------------------------
    const SectorInfo &sector_info(int s)               const { return sectors_[s]; }
    int  get_sector_id(double x, double y)             const;

    // --- quantized distance -------------------------------------------------
    void qdist(double x1, double y1, int s1,
               double x2, double y2, int s2,
               double &dx, double &dy)                 const;

    void qdist(const Module &m1, const Module &m2,
               double &dx, double &dy)                 const
    {
        qdist(m1.x, m1.y, m1.sector, m2.x, m2.y, m2.sector, dx, dy);
    }

    // --- static helpers -----------------------------------------------------
    static int          name_to_id(const std::string &name);
    static std::string  id_to_name(int id);
    static ModuleType   parse_type(const std::string &t);

private:
    void  compute_sectors();
    void  assign_layout(Module &m) const;
    void  build_neighbors();

    // line-segment intersection (ported from cana::intersection)
    static int line_intersect(double x1, double y1, double x2, double y2,
                              double x3, double y3, double x4, double y4,
                              double &xc, double &yc);

    int                  n_modules_ = 0;
    std::vector<Module>  modules_;

    // lookup maps
    std::unordered_map<std::string, int> name_map_;     // name → index
    std::unordered_map<int, int>         id_map_;       // primex_id → index
    std::unordered_map<uint64_t, int>    daq_map_;      // packed daq addr → index

    // sector info
    std::array<SectorInfo, static_cast<int>(Sector::Max)> sectors_;

    // pack DAQ address into a single key
    static uint64_t pack_daq(int crate, int slot, int ch)
    {
        return (static_cast<uint64_t>(crate) << 32) |
               (static_cast<uint64_t>(slot)  << 16) |
               static_cast<uint64_t>(ch);
    }
};

} // namespace fdec
