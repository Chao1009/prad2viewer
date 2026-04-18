// bind_det.cpp — pybind11 bindings for prad2det (prad2py.det submodule).
//
// Phase 2a: GEM reconstruction.
//   prad2py.det.GemSystem            (detector hierarchy + processing)
//   prad2py.det.GemCluster           (clustering algorithm)
//   prad2py.det.ClusterConfig        (tuning knobs for GemCluster)
//   prad2py.det.StripHit / StripCluster / GEMHit   (per-event outputs)
//   prad2py.det.ApvConfig / PlaneConfig / DetectorConfig / ApvPedestal
//
// Phases 2b (HyCal) and 2c (DetectorTransform + EpicsStore) will plug in
// alongside these — each gets its own `py::class_<...>` block below.
//
// Usage sketch (matches the C++ driver in test/gem_dump.cpp):
//
//     from prad2py import dec, det
//     cfg = dec.load_daq_config()
//     ch  = dec.EvChannel(); ch.set_config(cfg); ch.open(path)
//
//     gsys = det.GemSystem()
//     gsys.init("database/gem_map.json")
//     gsys.load_pedestals("database/gem_ped.json")    # optional
//     gcl = det.GemCluster()
//
//     while ch.read() == dec.Status.success:
//         if not ch.scan() or ch.get_event_type() != dec.EventType.Physics:
//             continue
//         for i in range(ch.get_n_events()):
//             ch.select_event(i)
//             ssp = ch.gem()                  # SspEventData
//             gsys.clear()
//             gsys.process_event(ssp)
//             gsys.reconstruct(gcl)
//             for h in gsys.get_all_hits():   # list[GEMHit]
//                 print(h.det_id, h.x, h.y, h.x_charge, h.y_charge)

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "GemSystem.h"
#include "GemCluster.h"
#include "HyCalSystem.h"
#include "HyCalCluster.h"
#include "DetectorTransform.h"
#include "EpicsStore.h"
#include "SspData.h"

#include <cstring>
#include <cstdio>
#include <string>

namespace py = pybind11;

// -------------------------------------------------------------------------
// GEM bindings
// -------------------------------------------------------------------------
static void bind_gem(py::module_ &m)
{
    // --- stateless strip-mapping (shared with GemSystem::buildStripMap) -----
    // Exposed at module level so Python layout / diagnostic scripts
    // (scripts/gem_strip_map.py, scripts/gem_layout.py) hit the same C++
    // implementation the reconstruction uses — no more drifting duplicates.

    m.def("map_strip", &gem::MapStrip,
          py::arg("ch"), py::arg("plane_index"), py::arg("orient"),
          py::arg("pin_rotate")     = 0,
          py::arg("shared_pos")     = -1,
          py::arg("hybrid_board")   = true,
          py::arg("apv_channels")   = 128,
          py::arg("readout_center") = 32,
          "Map one APV25 channel index to the plane-wide strip number.  "
          "Pure function — same 6-step pipeline used by "
          "GemSystem::buildStripMap.");

    m.def("map_apv_strips", &gem::MapApvStrips,
          py::arg("plane_index"), py::arg("orient"),
          py::arg("pin_rotate")     = 0,
          py::arg("shared_pos")     = -1,
          py::arg("hybrid_board")   = true,
          py::arg("apv_channels")   = 128,
          py::arg("readout_center") = 32,
          "Compute plane-wide strip numbers for every channel of one APV.  "
          "Returns a Python list of length `apv_channels`.");

    // --- configuration leaves ------------------------------------------------

    py::class_<gem::ApvPedestal>(m, "ApvPedestal",
        "Per-strip pedestal offset + noise RMS.")
        .def(py::init<>())
        .def_readwrite("offset", &gem::ApvPedestal::offset)
        .def_readwrite("noise",  &gem::ApvPedestal::noise);

    py::class_<gem::PlaneConfig>(m, "PlaneConfig",
        "One plane (X or Y) of a GEM detector.")
        .def_readwrite("type",   &gem::PlaneConfig::type)
        .def_readwrite("size",   &gem::PlaneConfig::size)
        .def_readwrite("n_apvs", &gem::PlaneConfig::n_apvs)
        .def_readwrite("pitch",  &gem::PlaneConfig::pitch);

    py::class_<gem::DetectorConfig>(m, "DetectorConfig",
        "One GEM detector — name, id, type, and its two planes.")
        .def_readwrite("name",   &gem::DetectorConfig::name)
        .def_readwrite("id",     &gem::DetectorConfig::id)
        .def_readwrite("type",   &gem::DetectorConfig::type)
        .def_property_readonly("plane_x",
            [](const gem::DetectorConfig &d) -> const gem::PlaneConfig& { return d.planes[0]; },
            py::return_value_policy::reference_internal)
        .def_property_readonly("plane_y",
            [](const gem::DetectorConfig &d) -> const gem::PlaneConfig& { return d.planes[1]; },
            py::return_value_policy::reference_internal);

    py::class_<gem::ApvConfig>(m, "ApvConfig",
        "Per-APV configuration: DAQ address + detector mapping + strip-mapping "
        "params + per-strip pedestals.  One entry per APV in the gem_map.json.")
        .def_readwrite("crate_id",     &gem::ApvConfig::crate_id)
        .def_readwrite("mpd_id",       &gem::ApvConfig::mpd_id)
        .def_readwrite("adc_ch",       &gem::ApvConfig::adc_ch)
        .def_readwrite("det_id",       &gem::ApvConfig::det_id)
        .def_readwrite("plane_type",   &gem::ApvConfig::plane_type)
        .def_readwrite("orient",       &gem::ApvConfig::orient)
        .def_readwrite("plane_index",  &gem::ApvConfig::plane_index)
        .def_readwrite("det_pos",      &gem::ApvConfig::det_pos)
        .def_readwrite("pin_rotate",   &gem::ApvConfig::pin_rotate)
        .def_readwrite("shared_pos",   &gem::ApvConfig::shared_pos)
        .def_readwrite("hybrid_board", &gem::ApvConfig::hybrid_board)
        .def_readwrite("match",        &gem::ApvConfig::match)
        .def_readwrite("cm_range_min", &gem::ApvConfig::cm_range_min)
        .def_readwrite("cm_range_max", &gem::ApvConfig::cm_range_max)
        .def("pedestal",
            [](const gem::ApvConfig &c, int ch) -> const gem::ApvPedestal& {
                if (ch < 0 || ch >= 128)
                    throw py::index_error("APV channel out of range [0,128)");
                return c.pedestal[ch];
            },
            py::arg("ch"),
            py::return_value_policy::reference_internal,
            "Pedestal (offset + noise) for one of 128 APV channels.");

    // --- per-event output types ---------------------------------------------

    py::class_<gem::StripHit>(m, "StripHit",
        "One strip's zero-suppressed pulse — one entry per channel with charge "
        "above the ZS threshold.  `ts_adc` holds the 6 pedestal/CM-corrected "
        "time samples.")
        .def_readonly("strip",       &gem::StripHit::strip)
        .def_readonly("charge",      &gem::StripHit::charge)
        .def_readonly("max_timebin", &gem::StripHit::max_timebin)
        .def_readonly("position",    &gem::StripHit::position)
        .def_readonly("cross_talk",  &gem::StripHit::cross_talk)
        .def_property_readonly("ts_adc",
            [](const gem::StripHit &h) {
                // Copy into a fresh numpy array — safer than viewing into
                // the C++ vector, which could dangle if the StripHit is
                // moved / the owning GemSystem is cleared.
                auto arr = py::array_t<float>(
                    static_cast<py::ssize_t>(h.ts_adc.size()));
                if (!h.ts_adc.empty())
                    std::memcpy(arr.mutable_data(),
                                h.ts_adc.data(),
                                h.ts_adc.size() * sizeof(float));
                return arr;
            },
            "Time-sample ADC values after pedestal + common-mode correction "
            "(numpy float32, one entry per SSP time sample).");

    py::class_<gem::StripCluster>(m, "StripCluster",
        "One 1-D plane cluster: a group of consecutive StripHits with "
        "charge-weighted position.")
        .def_readonly("position",     &gem::StripCluster::position)
        .def_readonly("peak_charge",  &gem::StripCluster::peak_charge)
        .def_readonly("total_charge", &gem::StripCluster::total_charge)
        .def_readonly("max_timebin",  &gem::StripCluster::max_timebin)
        .def_readonly("cross_talk",   &gem::StripCluster::cross_talk)
        .def_readonly("hits",         &gem::StripCluster::hits);

    py::class_<gem::GEMHit>(m, "GEMHit",
        "2-D reconstructed GEM hit (X cluster × Y cluster match).")
        .def_readonly("x",             &gem::GEMHit::x)
        .def_readonly("y",             &gem::GEMHit::y)
        .def_readonly("z",             &gem::GEMHit::z)
        .def_readonly("det_id",        &gem::GEMHit::det_id)
        .def_readonly("x_charge",      &gem::GEMHit::x_charge)
        .def_readonly("y_charge",      &gem::GEMHit::y_charge)
        .def_readonly("x_peak",        &gem::GEMHit::x_peak)
        .def_readonly("y_peak",        &gem::GEMHit::y_peak)
        .def_readonly("x_max_timebin", &gem::GEMHit::x_max_timebin)
        .def_readonly("y_max_timebin", &gem::GEMHit::y_max_timebin)
        .def_readonly("x_size",        &gem::GEMHit::x_size)
        .def_readonly("y_size",        &gem::GEMHit::y_size)
        .def_readonly("sig_pos",       &gem::GEMHit::sig_pos)
        .def("__repr__", [](const gem::GEMHit &h) {
            char buf[160];
            std::snprintf(buf, sizeof(buf),
                "<GEMHit det=%d x=%.3f y=%.3f z=%.3f qx=%.1f qy=%.1f>",
                h.det_id, h.x, h.y, h.z, h.x_charge, h.y_charge);
            return std::string(buf);
        });

    // --- GemCluster (configurable clustering algorithm) ---------------------

    py::class_<gem::ClusterConfig>(m, "ClusterConfig",
        "Tuning knobs for GemCluster.  Defaults reproduce the mpd_gem_view_ssp "
        "reconstruction chain.")
        .def(py::init<>())
        .def_readwrite("min_cluster_hits",   &gem::ClusterConfig::min_cluster_hits)
        .def_readwrite("max_cluster_hits",   &gem::ClusterConfig::max_cluster_hits)
        .def_readwrite("consecutive_thres",  &gem::ClusterConfig::consecutive_thres)
        .def_readwrite("split_thres",        &gem::ClusterConfig::split_thres)
        .def_readwrite("cross_talk_width",   &gem::ClusterConfig::cross_talk_width)
        .def_readwrite("position_res",       &gem::ClusterConfig::position_res)
        .def_readwrite("charac_dists",       &gem::ClusterConfig::charac_dists)
        .def_readwrite("match_mode",         &gem::ClusterConfig::match_mode)
        .def_readwrite("match_adc_asymmetry",&gem::ClusterConfig::match_adc_asymmetry)
        .def_readwrite("match_time_diff",    &gem::ClusterConfig::match_time_diff)
        .def_readwrite("ts_period",          &gem::ClusterConfig::ts_period);

    py::class_<gem::GemCluster>(m, "GemCluster",
        "Strip clustering + X/Y cluster matching.  Used by GemSystem during "
        "Reconstruct(); can also be driven directly for custom analysis.")
        .def(py::init<>())
        .def("set_config", &gem::GemCluster::SetConfig, py::arg("cfg"))
        .def("get_config", &gem::GemCluster::GetConfig,
             py::return_value_policy::reference_internal)
        .def("form_clusters",
            [](const gem::GemCluster &self, std::vector<gem::StripHit> &hits) {
                std::vector<gem::StripCluster> out;
                self.FormClusters(hits, out);
                return out;
            },
            py::arg("hits"),
            "Group the supplied StripHit list into StripClusters.  Modifies "
            "`hits` in place (sorted by strip).")
        .def("cartesian_reconstruct",
            [](const gem::GemCluster &self,
               const std::vector<gem::StripCluster> &x_clusters,
               const std::vector<gem::StripCluster> &y_clusters,
               int det_id, float resolution) {
                std::vector<gem::GEMHit> out;
                self.CartesianReconstruct(x_clusters, y_clusters, out,
                                          det_id, resolution);
                return out;
            },
            py::arg("x_clusters"), py::arg("y_clusters"),
            py::arg("det_id"), py::arg("resolution"),
            "Match X clusters with Y clusters on one detector to produce "
            "2-D GEMHits.");

    // --- GemSystem (the main entry point) -----------------------------------

    py::class_<gem::GemSystem>(m, "GemSystem",
        "PRad-II GEM detector system: loads the gem_map.json / gem_ped.json, "
        "processes SspEventData (pedestal subtraction, common-mode correction, "
        "zero suppression, strip mapping), and hands off to GemCluster for "
        "2-D reconstruction.")
        .def(py::init<>())

        // initialization
        .def("init",
            [](gem::GemSystem &self, const std::string &path) {
                py::gil_scoped_release rel;
                self.Init(path);
            },
            py::arg("map_file"),
            "Load the detector hierarchy and APV mapping from a JSON file "
            "(typically database/gem_map.json).")
        .def("load_pedestals",
            [](gem::GemSystem &self, const std::string &path) {
                py::gil_scoped_release rel;
                self.LoadPedestals(path);
            },
            py::arg("ped_file"),
            "Load per-strip pedestal mean/RMS.  Required before ProcessEvent "
            "for real zero suppression — defaults keep all strips silent.")
        .def("load_common_mode_range",
            [](gem::GemSystem &self, const std::string &path) {
                py::gil_scoped_release rel;
                self.LoadCommonModeRange(path);
            },
            py::arg("cm_file"),
            "Optional per-APV common-mode suppression window file.")

        // per-event processing
        .def("clear", &gem::GemSystem::Clear,
            "Reset per-event working buffers.  Call before every "
            "process_event().")
        .def("process_event",
            [](gem::GemSystem &self, const ssp::SspEventData &evt) {
                py::gil_scoped_release rel;
                self.ProcessEvent(evt);
            },
            py::arg("ssp_evt"),
            "Run pedestal + common-mode + zero-suppression over every APV in "
            "the given SspEventData.  Results feed reconstruct().")
        .def("reconstruct",
            [](gem::GemSystem &self, gem::GemCluster &cl) {
                py::gil_scoped_release rel;
                self.Reconstruct(cl);
            },
            py::arg("clusterer"),
            "Cluster the per-plane strip hits and match X/Y clusters into "
            "2-D GEMHits via the supplied GemCluster instance.")

        // accessors
        .def("get_n_detectors", &gem::GemSystem::GetNDetectors)
        .def("get_detectors",   &gem::GemSystem::GetDetectors,
             py::return_value_policy::reference_internal)
        .def("get_plane_hits", &gem::GemSystem::GetPlaneHits,
             py::arg("det"), py::arg("plane"),
             py::return_value_policy::reference_internal,
             "Per-plane StripHit list (plane=0 for X, 1 for Y).")
        .def("get_plane_clusters", &gem::GemSystem::GetPlaneClusters,
             py::arg("det"), py::arg("plane"),
             py::return_value_policy::reference_internal,
             "Per-plane StripCluster list (plane=0 for X, 1 for Y).  Populated "
             "after reconstruct().")
        .def("get_hits", &gem::GemSystem::GetHits, py::arg("det"),
             py::return_value_policy::reference_internal,
             "Per-detector reconstructed 2-D GEMHits.")
        .def("get_all_hits", &gem::GemSystem::GetAllHits,
             py::return_value_policy::reference_internal,
             "All GEMHits across every detector (flattened view).")

        // APV diagnostics
        .def("get_n_apvs", &gem::GemSystem::GetNApvs)
        .def("get_apv_config", &gem::GemSystem::GetApvConfig, py::arg("index"),
             py::return_value_policy::reference_internal)
        .def("find_apv_index", &gem::GemSystem::FindApvIndex,
             py::arg("crate"), py::arg("mpd"), py::arg("adc"),
             "O(1) DAQ→APV index lookup; -1 if no matching APV in the map.")
        .def("get_hole_x_offset", &gem::GemSystem::GetHoleXOffset,
             "Beam-hole X offset (mm) inferred from the `match` APVs in the map.")

        // post-process-event diagnostics
        .def("is_channel_hit", &gem::GemSystem::IsChannelHit,
             py::arg("apv_index"), py::arg("ch"),
             "True if this strip survived zero suppression in the last "
             "process_event().")
        .def("has_apv_zs_hits", &gem::GemSystem::HasApvZsHits,
             py::arg("apv_index"),
             "True if any channel in this APV survived zero suppression.")
        .def("get_processed_adc", &gem::GemSystem::GetProcessedAdc,
             py::arg("apv_index"), py::arg("ch"), py::arg("ts"),
             "Pedestal + common-mode-corrected ADC for (APV, channel, time "
             "sample); valid after process_event().")

        // threshold knobs
        .def_property("common_mode_threshold",
            &gem::GemSystem::GetCommonModeThreshold,
            &gem::GemSystem::SetCommonModeThreshold)
        .def_property("zero_sup_threshold",
            &gem::GemSystem::GetZeroSupThreshold,
            &gem::GemSystem::SetZeroSupThreshold)
        .def_property_readonly("cross_talk_threshold",
            &gem::GemSystem::GetCrossTalkThreshold);
}

// -------------------------------------------------------------------------
// HyCal bindings
// -------------------------------------------------------------------------
static void bind_hycal(py::module_ &m)
{
    // --- enums --------------------------------------------------------------

    py::enum_<fdec::ModuleType>(m, "ModuleType")
        .value("PbGlass", fdec::ModuleType::PbGlass)
        .value("PbWO4",   fdec::ModuleType::PbWO4)
        .value("LMS",     fdec::ModuleType::LMS)
        .value("Unknown", fdec::ModuleType::Unknown);

    py::enum_<fdec::Sector>(m, "Sector")
        .value("Center", fdec::Sector::Center)
        .value("Top",    fdec::Sector::Top)
        .value("Right",  fdec::Sector::Right)
        .value("Bottom", fdec::Sector::Bottom)
        .value("Left",   fdec::Sector::Left);

    // --- small leaf types ---------------------------------------------------

    py::class_<fdec::DaqAddr>(m, "DaqAddr",
        "HyCal module DAQ mapping: (crate, slot, channel).")
        .def_readwrite("crate",   &fdec::DaqAddr::crate)
        .def_readwrite("slot",    &fdec::DaqAddr::slot)
        .def_readwrite("channel", &fdec::DaqAddr::channel);

    py::class_<fdec::NeighborInfo>(m, "NeighborInfo",
        "One cross-sector neighbor of a HyCal module.")
        .def_readonly("index", &fdec::NeighborInfo::index)
        .def_readonly("dx",    &fdec::NeighborInfo::dx)
        .def_readonly("dy",    &fdec::NeighborInfo::dy)
        .def_readonly("dist",  &fdec::NeighborInfo::dist);

    py::class_<fdec::SectorInfo>(m, "SectorInfo",
        "HyCal sector (Center / Top / Right / Bottom / Left) geometry.")
        .def_readonly("id",      &fdec::SectorInfo::id)
        .def_readonly("mtype",   &fdec::SectorInfo::mtype)
        .def_readonly("msize_x", &fdec::SectorInfo::msize_x)
        .def_readonly("msize_y", &fdec::SectorInfo::msize_y)
        .def("boundary",
            [](const fdec::SectorInfo &s) {
                double x1, y1, x2, y2;
                s.get_boundary(x1, y1, x2, y2);
                return py::make_tuple(x1, y1, x2, y2);
            },
            "Return (x_min, y_min, x_max, y_max) of this sector's rectangle.");

    // --- Module -------------------------------------------------------------

    py::class_<fdec::Module>(m, "Module",
        "One HyCal module: identity, geometry, DAQ mapping, calibration, and "
        "pre-computed neighbor list.")
        .def_readonly("name",   &fdec::Module::name)
        .def_readonly("id",     &fdec::Module::id)
        .def_readonly("index",  &fdec::Module::index)
        .def_readonly("type",   &fdec::Module::type)
        .def_readonly("x",      &fdec::Module::x)
        .def_readonly("y",      &fdec::Module::y)
        .def_readonly("size_x", &fdec::Module::size_x)
        .def_readonly("size_y", &fdec::Module::size_y)
        .def_readonly("flag",   &fdec::Module::flag)
        .def_readonly("sector", &fdec::Module::sector)
        .def_readonly("row",    &fdec::Module::row)
        .def_readonly("column", &fdec::Module::column)
        .def_readonly("daq",    &fdec::Module::daq)
        .def_readonly("cal_factor",      &fdec::Module::cal_factor)
        .def_readonly("cal_base_energy", &fdec::Module::cal_base_energy)
        .def_readonly("cal_non_linear",  &fdec::Module::cal_non_linear)
        .def("energize", &fdec::Module::energize, py::arg("adc"),
             "Convert a pedestal-subtracted ADC value to MeV, including the "
             "non-linear correction term.")
        .def("is_pwo4",  &fdec::Module::is_pwo4)
        .def("is_glass", &fdec::Module::is_glass)
        .def("is_hycal", &fdec::Module::is_hycal)
        .def_property_readonly("neighbors",
            [](const fdec::Module &m) {
                py::list out;
                for (int i = 0; i < m.neighbor_count; ++i)
                    out.append(m.neighbors[i]);
                return out;
            },
            "Pre-computed cross-sector neighbor list (NeighborInfo[]).  "
            "Same-sector neighbors are resolved via the sector grid instead.")
        .def("__repr__", [](const fdec::Module &m) {
            char buf[96];
            std::snprintf(buf, sizeof(buf),
                "<Module %s id=%d idx=%d x=%.1f y=%.1f>",
                m.name.c_str(), m.id, m.index, m.x, m.y);
            return std::string(buf);
        });

    // --- HyCalSystem --------------------------------------------------------

    py::class_<fdec::HyCalSystem>(m, "HyCalSystem",
        "HyCal detector geometry + DAQ map + calibration.  Initialized once "
        "per job and then immutable — no per-event state lives here.")
        .def(py::init<>())
        .def("init",
            [](fdec::HyCalSystem &self, const std::string &modules_path,
               const std::string &daq_path) {
                py::gil_scoped_release rel;
                return self.Init(modules_path, daq_path);
            },
            py::arg("modules_path"), py::arg("daq_path"),
            "Load module geometry (hycal_modules.json) and DAQ map "
            "(daq_map.json).  Returns True on success.")
        .def("load_calibration",
            [](fdec::HyCalSystem &self, const std::string &path) {
                py::gil_scoped_release rel;
                return self.LoadCalibration(path);
            },
            py::arg("calib_path"),
            "Load per-module calibration constants from a JSON file.  "
            "Returns the number of modules matched, or -1 on error.")

        .def("module_count", &fdec::HyCalSystem::module_count)
        .def("module",
            static_cast<const fdec::Module &(fdec::HyCalSystem::*)(int) const>(
                &fdec::HyCalSystem::module),
            py::arg("index"),
            py::return_value_policy::reference_internal,
            "Module by array index (0 .. module_count()-1).")
        .def("module_by_name",
            [](const fdec::HyCalSystem &self, const std::string &name) -> py::object {
                const fdec::Module *m = self.module_by_name(name);
                if (!m) return py::none();
                return py::cast(m, py::return_value_policy::reference_internal);
            },
            py::arg("name"),
            "Module by name (e.g. 'W735'); None if not found.")
        .def("module_by_id",
            [](const fdec::HyCalSystem &self, int id) -> py::object {
                const fdec::Module *m = self.module_by_id(id);
                if (!m) return py::none();
                return py::cast(m, py::return_value_policy::reference_internal);
            },
            py::arg("primex_id"),
            "Module by PrimEx id (G: 1-576, W: 1001-2152); None if not found.")
        .def("module_by_daq",
            [](const fdec::HyCalSystem &self, int crate, int slot, int ch)
                -> py::object {
                const fdec::Module *m = self.module_by_daq(crate, slot, ch);
                if (!m) return py::none();
                return py::cast(m, py::return_value_policy::reference_internal);
            },
            py::arg("crate"), py::arg("slot"), py::arg("channel"),
            "Module by DAQ address; None if not found.")

        .def("get_calib_constant", &fdec::HyCalSystem::GetCalibConstant,
             py::arg("primex_id"))
        .def("set_calib_constant", &fdec::HyCalSystem::SetCalibConstant,
             py::arg("primex_id"), py::arg("factor"))
        .def("print_calib_constants", &fdec::HyCalSystem::PrintCalibConstants,
             py::arg("output_file"))

        .def("sector_info",
            &fdec::HyCalSystem::sector_info, py::arg("sector"),
            py::return_value_policy::reference_internal)
        .def("get_sector_id", &fdec::HyCalSystem::get_sector_id,
             py::arg("x"), py::arg("y"),
             "Return the sector id containing the (x, y) point, or -1 if none.")

        .def("qdist",
            [](const fdec::HyCalSystem &self,
               double x1, double y1, int s1,
               double x2, double y2, int s2) {
                double dx, dy;
                self.qdist(x1, y1, s1, x2, y2, s2, dx, dy);
                return py::make_tuple(dx, dy);
            },
            py::arg("x1"), py::arg("y1"), py::arg("s1"),
            py::arg("x2"), py::arg("y2"), py::arg("s2"),
            "Quantized (dx, dy) distance between two points across "
            "(potentially different) sectors.")
        .def("qdist_modules",
            [](const fdec::HyCalSystem &self,
               const fdec::Module &m1, const fdec::Module &m2) {
                double dx, dy;
                self.qdist(m1, m2, dx, dy);
                return py::make_tuple(dx, dy);
            },
            py::arg("m1"), py::arg("m2"),
            "Quantized (dx, dy) distance between two modules.")

        .def("for_each_neighbor",
            [](const fdec::HyCalSystem &self, int module_index,
               bool include_corners, const py::function &fn) {
                self.for_each_neighbor(module_index, include_corners,
                    [&](int ni) { fn(ni); });
            },
            py::arg("module_index"), py::arg("include_corners"),
            py::arg("fn"),
            "Walk every neighbor of `module_index`, same-sector (grid) + "
            "cross-sector (pre-computed list), calling fn(neighbor_index).")

        .def_static("name_to_id",  &fdec::HyCalSystem::name_to_id,
                    py::arg("name"))
        .def_static("id_to_name",  &fdec::HyCalSystem::id_to_name,
                    py::arg("id"))
        .def_static("parse_type",  &fdec::HyCalSystem::parse_type,
                    py::arg("type_string"));

    // --- clustering types ---------------------------------------------------

    py::class_<fdec::ClusterConfig>(m, "HyCalClusterConfig",
        "HyCal island-clustering tuning parameters.  Note: named "
        "`HyCalClusterConfig` in Python so it doesn't collide with the GEM "
        "`ClusterConfig`.")
        .def(py::init<>())
        .def_readwrite("min_module_energy",  &fdec::ClusterConfig::min_module_energy)
        .def_readwrite("min_center_energy",  &fdec::ClusterConfig::min_center_energy)
        .def_readwrite("min_cluster_energy", &fdec::ClusterConfig::min_cluster_energy)
        .def_readwrite("min_cluster_size",   &fdec::ClusterConfig::min_cluster_size)
        .def_readwrite("corner_conn",        &fdec::ClusterConfig::corner_conn)
        .def_readwrite("split_iter",         &fdec::ClusterConfig::split_iter)
        .def_readwrite("least_split",        &fdec::ClusterConfig::least_split)
        .def_readwrite("log_weight_thres",   &fdec::ClusterConfig::log_weight_thres);

    py::class_<fdec::ModuleHit>(m, "ModuleHit",
        "One HyCal module hit fed into the clusterer.")
        .def(py::init<>())
        .def_readwrite("index",  &fdec::ModuleHit::index)
        .def_readwrite("energy", &fdec::ModuleHit::energy);

    py::class_<fdec::ModuleCluster>(m, "ModuleCluster",
        "Module-level cluster: seed + constituent hits + total energy + flags.")
        .def_readonly("center", &fdec::ModuleCluster::center)
        .def_readonly("hits",   &fdec::ModuleCluster::hits)
        .def_readonly("energy", &fdec::ModuleCluster::energy)
        .def_readonly("flag",   &fdec::ModuleCluster::flag);

    py::class_<fdec::ClusterHit>(m, "ClusterHit",
        "Reconstructed HyCal cluster with (x, y) position in mm.")
        .def_readonly("center_id", &fdec::ClusterHit::center_id)
        .def_readonly("x",         &fdec::ClusterHit::x)
        .def_readonly("y",         &fdec::ClusterHit::y)
        .def_readonly("energy",    &fdec::ClusterHit::energy)
        .def_readonly("nblocks",   &fdec::ClusterHit::nblocks)
        .def_readonly("npos",      &fdec::ClusterHit::npos)
        .def_readonly("flag",      &fdec::ClusterHit::flag)
        .def("__repr__", [](const fdec::ClusterHit &c) {
            char buf[160];
            std::snprintf(buf, sizeof(buf),
                "<ClusterHit center_id=%d x=%.2f y=%.2f E=%.1f nblocks=%d>",
                c.center_id, c.x, c.y, c.energy, c.nblocks);
            return std::string(buf);
        });

    py::class_<fdec::HyCalCluster::RecoResult>(m, "HyCalRecoResult",
        "Paired ModuleCluster + reconstructed ClusterHit from "
        "HyCalCluster.reconstruct_matched().  `cluster` is a reference back "
        "into HyCalCluster's internal list — invalidated on clear().")
        .def_property_readonly("cluster",
            [](const fdec::HyCalCluster::RecoResult &r) -> const fdec::ModuleCluster& {
                return *r.cluster;
            },
            py::return_value_policy::reference_internal)
        .def_readonly("hit", &fdec::HyCalCluster::RecoResult::hit);

    // --- HyCalCluster -------------------------------------------------------

    py::class_<fdec::HyCalCluster>(m, "HyCalCluster",
        "Island clustering on top of a HyCalSystem: accepts (module_index, "
        "energy) hits per event, groups connected modules, splits multi-"
        "maximum groups, and produces ClusterHit positions via log-weighted "
        "center-of-gravity.")
        .def(py::init<const fdec::HyCalSystem &>(), py::arg("system"),
             py::keep_alive<1, 2>(),
             "Construct with a reference to an initialized HyCalSystem — "
             "the system must outlive this clusterer.")
        .def("set_config", &fdec::HyCalCluster::SetConfig, py::arg("cfg"))
        .def("get_config", &fdec::HyCalCluster::GetConfig,
             py::return_value_policy::reference_internal)
        .def("clear", &fdec::HyCalCluster::Clear,
             "Reset per-event state.  Call before add_hit()s for each event.")
        .def("add_hit", &fdec::HyCalCluster::AddHit,
             py::arg("module_index"), py::arg("energy"),
             "Add a hit for the module at `module_index` with the given "
             "calibrated energy (MeV).")
        .def("form_clusters",
            [](fdec::HyCalCluster &self) {
                py::gil_scoped_release rel;
                self.FormClusters();
            },
            "Run the island grouping + splitting algorithm over the hits "
            "accumulated via add_hit().")
        .def("reconstruct_hits",
            [](const fdec::HyCalCluster &self) {
                std::vector<fdec::ClusterHit> out;
                { py::gil_scoped_release rel; self.ReconstructHits(out); }
                return out;
            },
            "Return the list of ClusterHits (x, y, energy, nblocks, …) for "
            "every cluster that passed min_cluster_energy / min_cluster_size.")
        .def("reconstruct_matched",
            [](const fdec::HyCalCluster &self) {
                std::vector<fdec::HyCalCluster::RecoResult> out;
                { py::gil_scoped_release rel; self.ReconstructMatched(out); }
                return out;
            },
            "Return [(ModuleCluster, ClusterHit)] pairs for every cluster "
            "that passed thresholds — avoids the fragile parallel-iteration "
            "between get_clusters() and reconstruct_hits().")
        .def("get_clusters", &fdec::HyCalCluster::GetClusters,
             py::return_value_policy::reference_internal,
             "Low-level module-level clusters (ModuleCluster[]).  Only valid "
             "between form_clusters() and the next clear().");
}

// -------------------------------------------------------------------------
// Helper bindings (Phase 2c)
// -------------------------------------------------------------------------
static void bind_transform(py::module_ &m)
{
    py::class_<DetectorTransform::Matrix>(m, "TransformMatrix",
        "Pre-computed 2×3 rotation + translation matrix cached inside a "
        "DetectorTransform.")
        .def_readonly("r00", &DetectorTransform::Matrix::r00)
        .def_readonly("r01", &DetectorTransform::Matrix::r01)
        .def_readonly("r10", &DetectorTransform::Matrix::r10)
        .def_readonly("r11", &DetectorTransform::Matrix::r11)
        .def_readonly("r20", &DetectorTransform::Matrix::r20)
        .def_readonly("r21", &DetectorTransform::Matrix::r21)
        .def_readonly("tx",  &DetectorTransform::Matrix::tx)
        .def_readonly("ty",  &DetectorTransform::Matrix::ty)
        .def_readonly("tz",  &DetectorTransform::Matrix::tz);

    py::class_<DetectorTransform>(m, "DetectorTransform",
        "Planar detector pose (origin + tilts in degrees).  to_lab(x, y) "
        "returns a lab-frame 3-vector; rotate(x, y) skips the translation.")
        .def(py::init<>())
        .def_readwrite("x",  &DetectorTransform::x)
        .def_readwrite("y",  &DetectorTransform::y)
        .def_readwrite("z",  &DetectorTransform::z)
        .def_readwrite("rx", &DetectorTransform::rx)
        .def_readwrite("ry", &DetectorTransform::ry)
        .def_readwrite("rz", &DetectorTransform::rz)
        .def("prepare", &DetectorTransform::prepare,
             "Force matrix precomputation (idempotent).  Called implicitly "
             "by to_lab / rotate / matrix().")
        .def("to_lab",
            [](const DetectorTransform &self, float dx, float dy) {
                float lx, ly, lz;
                self.toLab(dx, dy, lx, ly, lz);
                return py::make_tuple(lx, ly, lz);
            },
            py::arg("dx"), py::arg("dy"),
            "Map a (dx, dy) point on the detector plane to lab-frame "
            "(x, y, z).")
        .def("rotate",
            [](const DetectorTransform &self, float dx, float dy) {
                float ox, oy;
                self.rotate(dx, dy, ox, oy);
                return py::make_tuple(ox, oy);
            },
            py::arg("dx"), py::arg("dy"),
            "Rotate a (dx, dy) vector by the detector tilts — no "
            "translation.  Handy for drawing in detector-local coords.")
        .def("matrix", &DetectorTransform::matrix,
             py::return_value_policy::reference_internal);
}

static void bind_epics(py::module_ &m)
{
    py::class_<fdec::EpicsStore::Snapshot>(m, "EpicsSnapshot",
        "One EPICS snapshot (event number + timestamp + channel values).")
        .def_readonly("event_number", &fdec::EpicsStore::Snapshot::event_number)
        .def_readonly("timestamp",    &fdec::EpicsStore::Snapshot::timestamp)
        .def_readonly("values",       &fdec::EpicsStore::Snapshot::values);

    py::class_<fdec::EpicsStore>(m, "EpicsStore",
        "Accumulator for EPICS slow-control snapshots.  Feed it raw EPICS "
        "text via feed() whenever an EPICS event arrives (channel discovery "
        "is automatic), then query get_value() for any subsequent physics "
        "event — returns the most recent snapshot at or before that event.")
        .def(py::init<>())
        .def("feed", &fdec::EpicsStore::Feed,
             py::arg("event_number"), py::arg("timestamp"), py::arg("text"),
             "Parse a raw EPICS payload and store a snapshot.  Text format: "
             "one `value  channel_name` line per channel.")
        .def("get_value",
            [](const fdec::EpicsStore &self, int32_t event_number,
               const std::string &channel) -> py::object {
                float v = 0.f;
                if (self.GetValue(event_number, channel, v))
                    return py::float_(v);
                return py::none();
            },
            py::arg("event_number"), py::arg("channel"),
            "Return the most recent value of `channel` at or before "
            "`event_number`, or None if unknown / not yet seen.")
        .def("find_snapshot",
            [](const fdec::EpicsStore &self, int32_t event_number) -> py::object {
                const auto *s = self.FindSnapshot(event_number);
                if (!s) return py::none();
                return py::cast(s, py::return_value_policy::reference_internal);
            },
            py::arg("event_number"),
            "Return the whole snapshot at or before `event_number`, or None.")
        .def("get_channel_count", &fdec::EpicsStore::GetChannelCount)
        .def("get_channel_id",    &fdec::EpicsStore::GetChannelId,
             py::arg("name"),
             "Return the integer id assigned to `name`, or -1 if unknown.")
        .def("get_channel_name",  &fdec::EpicsStore::GetChannelName,
             py::arg("id"), py::return_value_policy::reference_internal)
        .def("get_channel_names", &fdec::EpicsStore::GetChannelNames,
             py::return_value_policy::reference_internal)
        .def("get_snapshot_count", &fdec::EpicsStore::GetSnapshotCount)
        .def("get_snapshot",       &fdec::EpicsStore::GetSnapshot,
             py::arg("index"), py::return_value_policy::reference_internal)
        .def("trim", &fdec::EpicsStore::Trim, py::arg("max_count"),
             "Drop oldest snapshots so at most `max_count` remain.")
        .def("clear", &fdec::EpicsStore::Clear);
}

// -------------------------------------------------------------------------
// Submodule entry point — called from prad2py.cpp
// -------------------------------------------------------------------------
void register_det(py::module_ &m)
{
    auto det = m.def_submodule("det",
        "prad2det bindings — GEM + HyCal reconstruction and slow-control "
        "helpers.");

    bind_gem(det);       // 2a
    bind_hycal(det);     // 2b
    bind_transform(det); // 2c
    bind_epics(det);     // 2c
}
