// bind_dec.cpp — pybind11 bindings for prad2dec (prad2py.dec submodule)
//
// Exposes:
//   prad2py.dec.EventType            (enum)
//   prad2py.dec.Status               (enum)
//   prad2py.dec.DaqConfig            (config struct)
//   prad2py.dec.EventInfo            (per-event metadata)
//   prad2py.dec.ChannelData / SlotData / RocData / EventData  (fdec)
//   prad2py.dec.ApvAddress / ApvData / MpdData / SspEventData (ssp)
//   prad2py.dec.EcPeak / EcCluster / VtpBlock / VtpEventData  (vtp)
//   prad2py.dec.TdcHit / TdcEventData                         (tdc)
//   prad2py.dec.SyncInfo                                      (sync)
//   prad2py.dec.EvChannel            (evio reader)
//   prad2py.dec.load_daq_config(path) -> DaqConfig
//
// The bulk per-channel arrays (ChannelData.samples, ApvData.strips, etc.)
// are returned as numpy arrays — copies by default so the buffer stays
// valid after the next DecodeEvent call.

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "EvChannel.h"
#include "DaqConfig.h"
#include "load_daq_config.h"
#include "Fadc250Data.h"
#include "WaveAnalyzer.h"
#include "Fadc250FwAnalyzer.h"
#include "SspData.h"
#include "VtpData.h"
#include "TdcData.h"
#include "SyncData.h"

#include <cstdlib>
#include <memory>
#include <string>

#include "InstallPaths.h"

namespace py = pybind11;

#ifndef DATABASE_DIR
#define DATABASE_DIR "."
#endif

namespace {

std::string default_daq_config_path()
{
    // Same lookup policy as prad2py.cpp's copy (these two TUs each have
    // a `load_daq_config("")` overload that needs the default path).
    std::string dir = prad2::resolve_data_dir(
        "PRAD2_DATABASE_DIR",
        {"../../share/prad2evviewer/database",
         "../share/prad2evviewer/database"},
        DATABASE_DIR);
    return dir + "/daq_config.json";
}

// -------------------------------------------------------------------------
// FADC250 data types (fdec::)
// -------------------------------------------------------------------------
void bind_fadc(py::module_ &m)
{
    // ChannelData — bulk samples are exposed as a numpy view over the
    // first `nsamples` elements.  The buffer belongs to the owning event,
    // so callers must copy (default) to keep data across DecodeEvent calls.
    py::class_<fdec::ChannelData>(m, "ChannelData",
        "Per-channel FADC250 samples (nsamples uint16 values).")
        .def_readonly("nsamples", &fdec::ChannelData::nsamples)
        .def_property_readonly("samples",
            [](const fdec::ChannelData &c) {
                // Always copy — the owning EventData may be reused/overwritten.
                return py::array_t<uint16_t>(
                    {c.nsamples},
                    {static_cast<py::ssize_t>(sizeof(uint16_t))},
                    c.samples);
            },
            "16-bit ADC samples as a fresh numpy array (copy of nsamples values).");

    py::class_<fdec::SlotData>(m, "SlotData",
        "One FADC250 slot: trigger info plus per-channel samples.")
        .def_readonly("present",      &fdec::SlotData::present)
        .def_readonly("trigger",      &fdec::SlotData::trigger)
        .def_readonly("timestamp",    &fdec::SlotData::timestamp)
        .def_readonly("nchannels",    &fdec::SlotData::nchannels)
        .def_readonly("channel_mask", &fdec::SlotData::channel_mask)
        .def("channel",
            [](const fdec::SlotData &s, int ch) -> const fdec::ChannelData& {
                if (ch < 0 || ch >= fdec::MAX_CHANNELS)
                    throw py::index_error("channel out of range");
                return s.channels[ch];
            },
            py::arg("channel"),
            py::return_value_policy::reference_internal,
            "Access ChannelData by channel index (no presence check — use channel_mask).")
        .def("present_channels",
            [](const fdec::SlotData &s) {
                std::vector<int> out;
                out.reserve(s.nchannels);
                for (int c = 0; c < fdec::MAX_CHANNELS; ++c)
                    if (s.channel_mask & (1ULL << c)) out.push_back(c);
                return out;
            },
            "List of channel indices with nsamples > 0 this event.");

    py::class_<fdec::RocData>(m, "RocData",
        "One ROC crate worth of FADC data.")
        .def_readonly("present", &fdec::RocData::present)
        .def_readonly("tag",     &fdec::RocData::tag)
        .def_readonly("nslots",  &fdec::RocData::nslots)
        .def("slot",
            [](const fdec::RocData &r, int sl) -> const fdec::SlotData& {
                if (sl < 0 || sl >= fdec::MAX_SLOTS)
                    throw py::index_error("slot out of range");
                return r.slots[sl];
            },
            py::arg("slot"),
            py::return_value_policy::reference_internal)
        .def("present_slots",
            [](const fdec::RocData &r) {
                std::vector<int> out;
                for (int s = 0; s < fdec::MAX_SLOTS; ++s)
                    if (r.slots[s].present) out.push_back(s);
                return out;
            });

    py::class_<fdec::EventInfo>(m, "EventInfo",
        "Event-level metadata (TI + trigger bank).")
        .def_readwrite("type",            &fdec::EventInfo::type)
        .def_readwrite("trigger_type",    &fdec::EventInfo::trigger_type)
        .def_readwrite("trigger_bits",    &fdec::EventInfo::trigger_bits)
        .def_readwrite("event_tag",       &fdec::EventInfo::event_tag)
        .def_readwrite("event_number",    &fdec::EventInfo::event_number)
        .def_readwrite("trigger_number",  &fdec::EventInfo::trigger_number)
        .def_readwrite("timestamp",       &fdec::EventInfo::timestamp)
        .def_readwrite("run_number",      &fdec::EventInfo::run_number)
        .def_readwrite("unix_time",       &fdec::EventInfo::unix_time)
        .def("__repr__", [](const fdec::EventInfo &i) {
            char buf[256];
            std::snprintf(buf, sizeof(buf),
                "<EventInfo evt=%d trig=%d tag=0x%x bits=0x%x ts=%llu>",
                i.event_number, i.trigger_number, i.event_tag, i.trigger_bits,
                (unsigned long long)i.timestamp);
            return std::string(buf);
        });

    py::class_<psync::SyncInfo>(m, "SyncInfo",
        "Absolute-time / run-state snapshot populated from SYNC/EPICS 0xE112 "
        "HEAD banks and from PRESTART/GO/END control-event payloads.  Persists "
        "across events in EvChannel — see EvChannel.sync().")
        .def_readwrite("run_number",   &psync::SyncInfo::run_number)
        .def_readwrite("sync_counter", &psync::SyncInfo::sync_counter)
        .def_readwrite("unix_time",    &psync::SyncInfo::unix_time)
        .def_readwrite("event_tag",    &psync::SyncInfo::event_tag)
        .def_readwrite("run_type",     &psync::SyncInfo::run_type)
        .def("valid", &psync::SyncInfo::valid,
             "True if unix_time has been populated (any SYNC/control event seen).")
        .def("__repr__", [](const psync::SyncInfo &s) {
            char buf[160];
            std::snprintf(buf, sizeof(buf),
                "<SyncInfo run=%u counter=%u unix=%u tag=0x%x type=%u>",
                s.run_number, s.sync_counter, s.unix_time, s.event_tag,
                (unsigned)s.run_type);
            return std::string(buf);
        });

    py::class_<fdec::EventData, std::shared_ptr<fdec::EventData>>(m, "EventData",
        "Full decoded event: EventInfo + FADC ROCs.")
        .def(py::init<>())
        .def_readonly("info",   &fdec::EventData::info)
        .def_readonly("nrocs",  &fdec::EventData::nrocs)
        .def("roc",
            [](const fdec::EventData &e, int i) -> const fdec::RocData& {
                if (i < 0 || i >= e.nrocs)
                    throw py::index_error("ROC index out of range");
                return e.rocs[e.roc_index[i]];
            },
            py::arg("index"),
            py::return_value_policy::reference_internal,
            "Access the i-th active ROC by slot in roc_index[].")
        .def("find_roc",
            [](const fdec::EventData &e, uint32_t tag) -> const fdec::RocData* {
                return e.findRoc(tag);   // nullptr → None via pybind11
            },
            py::arg("tag"),
            py::return_value_policy::reference_internal,
            "Locate a ROC by its bank tag (e.g. 0x80). Returns None if absent.")
        .def("clear", &fdec::EventData::clear);

    // ----- WaveAnalyzer types --------------------------------------------
    py::class_<fdec::Peak>(m, "Peak",
        "One FADC pulse found by WaveAnalyzer.")
        .def(py::init<>())
        .def_readwrite("height",   &fdec::Peak::height)
        .def_readwrite("integral", &fdec::Peak::integral)
        .def_readwrite("time",     &fdec::Peak::time)
        .def_readwrite("pos",      &fdec::Peak::pos)
        .def_readwrite("left",     &fdec::Peak::left)
        .def_readwrite("right",    &fdec::Peak::right)
        .def_readwrite("overflow", &fdec::Peak::overflow);

    py::class_<fdec::Pedestal>(m, "Pedestal",
        "Per-channel pedestal (mean + RMS).")
        .def(py::init<>())
        .def_readwrite("mean", &fdec::Pedestal::mean)
        .def_readwrite("rms",  &fdec::Pedestal::rms);

    py::class_<fdec::WaveConfig>(m, "WaveConfig",
        "Knobs for WaveAnalyzer (smoothing, thresholds, pedestal window, ...).")
        .def(py::init<>())
        .def_readwrite("resolution",     &fdec::WaveConfig::resolution)
        .def_readwrite("threshold",      &fdec::WaveConfig::threshold)
        .def_readwrite("min_threshold",  &fdec::WaveConfig::min_threshold)
        .def_readwrite("min_peak_ratio", &fdec::WaveConfig::min_peak_ratio)
        .def_readwrite("int_tail_ratio", &fdec::WaveConfig::int_tail_ratio)
        .def_readwrite("ped_nsamples",   &fdec::WaveConfig::ped_nsamples)
        .def_readwrite("ped_flatness",   &fdec::WaveConfig::ped_flatness)
        .def_readwrite("ped_max_iter",   &fdec::WaveConfig::ped_max_iter)
        .def_readwrite("overflow",       &fdec::WaveConfig::overflow)
        .def_readwrite("clk_mhz",        &fdec::WaveConfig::clk_mhz);

    py::class_<fdec::WaveAnalyzer>(m, "WaveAnalyzer",
        "Fast C++ pedestal + peak finder used by the server.  Call "
        "analyze(samples_numpy) to get (ped_mean, ped_rms, [Peak, ...]).")
        .def(py::init<const fdec::WaveConfig &>(),
             py::arg("cfg") = fdec::WaveConfig{})
        .def_readwrite("cfg", &fdec::WaveAnalyzer::cfg)
        .def("analyze",
            [](const fdec::WaveAnalyzer &self, py::array_t<uint16_t> samples) {
                py::buffer_info buf = samples.request();
                if (buf.ndim != 1)
                    throw py::value_error("samples must be a 1-D uint16 array");
                fdec::WaveResult res;
                {
                    py::gil_scoped_release rel;
                    self.Analyze(static_cast<const uint16_t*>(buf.ptr),
                                 static_cast<int>(buf.shape[0]), res);
                }
                py::list peaks;
                for (int i = 0; i < res.npeaks; ++i)
                    peaks.append(res.peaks[i]);
                return py::make_tuple(res.ped.mean, res.ped.rms, peaks);
            },
            py::arg("samples"),
            "Analyze one channel's FADC waveform (uint16 numpy array). "
            "Returns (pedestal_mean, pedestal_rms, peaks_list).  Peaks "
            "above ``cfg.min_threshold`` / ``cfg.threshold × rms`` are "
            "kept; up to ``MAX_PEAKS`` per call.");

    // ----- Firmware-faithful (Mode 1/2/3) analyzer -----------------------
    // Quality bitmask constants — exposed at module scope so callers can
    // test pulse.quality & dec.Q_PEAK_AT_BOUNDARY etc.
    m.attr("Q_GOOD")             = py::int_(fdec::Q_GOOD);
    m.attr("Q_PEAK_AT_BOUNDARY") = py::int_(fdec::Q_PEAK_AT_BOUNDARY);
    m.attr("Q_NSB_TRUNCATED")    = py::int_(fdec::Q_NSB_TRUNCATED);
    m.attr("Q_NSA_TRUNCATED")    = py::int_(fdec::Q_NSA_TRUNCATED);
    m.attr("Q_VA_OUT_OF_RANGE")  = py::int_(fdec::Q_VA_OUT_OF_RANGE);

    py::class_<fdec::DaqPeak>(m, "DaqPeak",
        "One firmware-mode pulse (Mode 1 + Mode 2 + Mode 3 result).")
        .def(py::init<>())
        .def_readwrite("pulse_id",     &fdec::DaqPeak::pulse_id)
        .def_readwrite("vmin",         &fdec::DaqPeak::vmin)
        .def_readwrite("vpeak",        &fdec::DaqPeak::vpeak)
        .def_readwrite("va",           &fdec::DaqPeak::va)
        .def_readwrite("coarse",       &fdec::DaqPeak::coarse)
        .def_readwrite("fine",         &fdec::DaqPeak::fine)
        .def_readwrite("time_units",   &fdec::DaqPeak::time_units)
        .def_readwrite("time_ns",      &fdec::DaqPeak::time_ns)
        .def_readwrite("cross_sample", &fdec::DaqPeak::cross_sample)
        .def_readwrite("peak_sample",  &fdec::DaqPeak::peak_sample)
        .def_readwrite("integral",     &fdec::DaqPeak::integral)
        .def_readwrite("window_lo",    &fdec::DaqPeak::window_lo)
        .def_readwrite("window_hi",    &fdec::DaqPeak::window_hi)
        .def_readwrite("quality",      &fdec::DaqPeak::quality);

    py::class_<evc::DaqConfig::Fadc250FwConfig>(m, "Fadc250FwConfig",
        "Firmware Mode 1/2/3 emulation knobs.  Names follow the Hall-D V3 "
        "firmware API (faV3HallDSetProcMode).  NPEAK is an alias for "
        "MAX_PULSES.")
        .def(py::init<>())
        .def_readwrite("TET",        &evc::DaqConfig::Fadc250FwConfig::TET)
        .def_readwrite("NSB",        &evc::DaqConfig::Fadc250FwConfig::NSB)
        .def_readwrite("NSA",        &evc::DaqConfig::Fadc250FwConfig::NSA)
        .def_readwrite("MAX_PULSES", &evc::DaqConfig::Fadc250FwConfig::MAX_PULSES)
        .def_readwrite("NSAT",       &evc::DaqConfig::Fadc250FwConfig::NSAT)
        .def_readwrite("NPED",       &evc::DaqConfig::Fadc250FwConfig::NPED)
        .def_readwrite("MAXPED",     &evc::DaqConfig::Fadc250FwConfig::MAXPED)
        .def_readwrite("CLK_NS",     &evc::DaqConfig::Fadc250FwConfig::CLK_NS);

    py::class_<fdec::Fadc250FwAnalyzer>(m, "Fadc250FwAnalyzer",
        "Firmware-faithful FADC250 Mode 1/2/3 emulator.  Mirrors the on-board "
        "TDC + pulse-windowing exactly (per FADC250 User's Manual).  Use this "
        "to compare against firmware-reported pulse data, or for DAQ signal "
        "studies where the soft (local-maxima) analyzer would diverge.")
        .def(py::init<const evc::DaqConfig::Fadc250FwConfig &>(),
             py::arg("cfg") = evc::DaqConfig::Fadc250FwConfig{})
        .def_readwrite("cfg", &fdec::Fadc250FwAnalyzer::cfg)
        .def("analyze",
            [](const fdec::Fadc250FwAnalyzer &self,
               py::array_t<uint16_t> samples, float ped) {
                py::buffer_info buf = samples.request();
                if (buf.ndim != 1)
                    throw py::value_error("samples must be a 1-D uint16 array");
                fdec::DaqWaveResult res;
                {
                    py::gil_scoped_release rel;
                    self.Analyze(static_cast<const uint16_t*>(buf.ptr),
                                 static_cast<int>(buf.shape[0]), ped, res);
                }
                py::list peaks;
                for (int i = 0; i < res.npeaks; ++i)
                    peaks.append(res.peaks[i]);
                return py::make_tuple(res.vnoise, peaks);
            },
            py::arg("samples"), py::arg("ped"),
            "Run firmware Mode 3 (TDC) → Mode 1 + Mode 2 windowing on a "
            "uint16 sample array.  Returns (vnoise, [DaqPeak, ...]).  ``ped`` "
            "is the per-channel pedestal (firmware register or soft-analyzer "
            "estimate); samples have not been pedestal-subtracted.");
}

// -------------------------------------------------------------------------
// SSP / MPD / APV (ssp::)
// -------------------------------------------------------------------------
void bind_ssp(py::module_ &m)
{
    py::class_<ssp::ApvAddress>(m, "ApvAddress",
        "APV hardware identifier (crate, MPD/fiber, ADC channel).")
        .def(py::init<>())
        .def_readwrite("crate_id", &ssp::ApvAddress::crate_id)
        .def_readwrite("mpd_id",   &ssp::ApvAddress::mpd_id)
        .def_readwrite("adc_ch",   &ssp::ApvAddress::adc_ch)
        .def("__repr__", [](const ssp::ApvAddress &a) {
            char buf[96];
            std::snprintf(buf, sizeof(buf),
                "<ApvAddress crate=%d mpd=%d adc=%d>",
                a.crate_id, a.mpd_id, a.adc_ch);
            return std::string(buf);
        });

    py::class_<ssp::ApvData>(m, "ApvData",
        "One APV chip's strip/time-sample matrix.")
        .def_readonly("addr",          &ssp::ApvData::addr)
        .def_readonly("present",       &ssp::ApvData::present)
        .def_readonly("nstrips",       &ssp::ApvData::nstrips)
        .def_readonly("flags",         &ssp::ApvData::flags)
        .def_readonly("has_online_cm", &ssp::ApvData::has_online_cm)
        .def_property_readonly("strips",
            [](const ssp::ApvData &a) {
                // [APV_STRIP_SIZE][SSP_TIME_SAMPLES] int16 — fresh copy.
                return py::array_t<int16_t>(
                    {static_cast<py::ssize_t>(ssp::APV_STRIP_SIZE),
                     static_cast<py::ssize_t>(ssp::SSP_TIME_SAMPLES)},
                    {static_cast<py::ssize_t>(sizeof(int16_t) * ssp::SSP_TIME_SAMPLES),
                     static_cast<py::ssize_t>(sizeof(int16_t))},
                    &a.strips[0][0]);
            },
            "(128, 6) int16 numpy array of raw ADC samples (fresh copy).")
        .def_property_readonly("online_cm",
            [](const ssp::ApvData &a) {
                return py::array_t<int16_t>(
                    {static_cast<py::ssize_t>(ssp::SSP_TIME_SAMPLES)},
                    {static_cast<py::ssize_t>(sizeof(int16_t))},
                    a.online_cm);
            },
            "6-element online common-mode values (fresh copy).")
        .def("has_strip", &ssp::ApvData::hasStrip);

    py::class_<ssp::MpdData>(m, "MpdData",
        "One MPD card: a vector of APVs.")
        .def_readonly("crate_id", &ssp::MpdData::crate_id)
        .def_readonly("mpd_id",   &ssp::MpdData::mpd_id)
        .def_readonly("present",  &ssp::MpdData::present)
        .def_readonly("napvs",    &ssp::MpdData::napvs)
        .def("apv",
            [](const ssp::MpdData &m, int adc) -> const ssp::ApvData& {
                if (adc < 0 || adc >= ssp::MAX_APVS_PER_MPD)
                    throw py::index_error("APV index out of range");
                return m.apvs[adc];
            },
            py::return_value_policy::reference_internal);

    py::class_<ssp::SspEventData, std::shared_ptr<ssp::SspEventData>>(m, "SspEventData",
        "Per-event GEM strip data grouped by (crate, MPD).")
        .def(py::init<>())
        .def_readonly("nmpds", &ssp::SspEventData::nmpds)
        .def("mpd",
            [](const ssp::SspEventData &e, int i) -> const ssp::MpdData& {
                if (i < 0 || i >= e.nmpds)
                    throw py::index_error("MPD index out of range");
                return e.mpds[i];
            },
            py::return_value_policy::reference_internal)
        .def("find_apv",
            [](const ssp::SspEventData &e, int crate, int mpd, int adc)
                -> const ssp::ApvData* {
                // Return a raw pointer and let the .def handler do the
                // cast — manual py::cast(..., reference_internal) loses
                // the parent linkage and keep_alive<0,1> fails with
                // "Could not activate keep_alive".  nullptr → None is
                // handled automatically by pybind11.
                return e.findApv(crate, mpd, adc);
            },
            py::arg("crate"), py::arg("mpd"), py::arg("adc"),
            py::return_value_policy::reference_internal)
        .def("clear", &ssp::SspEventData::clear);
}

// -------------------------------------------------------------------------
// VTP (vtp::)
// -------------------------------------------------------------------------
void bind_vtp(py::module_ &m)
{
    py::class_<vtp::EcPeak>(m, "EcPeak")
        .def_readonly("roc_tag", &vtp::EcPeak::roc_tag)
        .def_readonly("inst",    &vtp::EcPeak::inst)
        .def_readonly("view",    &vtp::EcPeak::view)
        .def_readonly("time",    &vtp::EcPeak::time)
        .def_readonly("coord",   &vtp::EcPeak::coord)
        .def_readonly("energy",  &vtp::EcPeak::energy);

    py::class_<vtp::EcCluster>(m, "EcCluster")
        .def_readonly("roc_tag", &vtp::EcCluster::roc_tag)
        .def_readonly("inst",    &vtp::EcCluster::inst)
        .def_readonly("time",    &vtp::EcCluster::time)
        .def_readonly("energy",  &vtp::EcCluster::energy)
        .def_readonly("coordU",  &vtp::EcCluster::coordU)
        .def_readonly("coordV",  &vtp::EcCluster::coordV)
        .def_readonly("coordW",  &vtp::EcCluster::coordW);

    py::class_<vtp::VtpBlock>(m, "VtpBlock")
        .def_readonly("roc_tag",          &vtp::VtpBlock::roc_tag)
        .def_readonly("slot",             &vtp::VtpBlock::slot)
        .def_readonly("module_id",        &vtp::VtpBlock::module_id)
        .def_readonly("block_number",     &vtp::VtpBlock::block_number)
        .def_readonly("block_level",      &vtp::VtpBlock::block_level)
        .def_readonly("nwords",           &vtp::VtpBlock::nwords)
        .def_readonly("event_number",     &vtp::VtpBlock::event_number)
        .def_readonly("trigger_time",     &vtp::VtpBlock::trigger_time)
        .def_readonly("has_trailer",      &vtp::VtpBlock::has_trailer)
        .def_readonly("trailer_mismatch", &vtp::VtpBlock::trailer_mismatch);

    py::class_<vtp::VtpEventData, std::shared_ptr<vtp::VtpEventData>>(m, "VtpEventData")
        .def(py::init<>())
        .def_readonly("n_peaks",    &vtp::VtpEventData::n_peaks)
        .def_readonly("n_clusters", &vtp::VtpEventData::n_clusters)
        .def_readonly("n_blocks",   &vtp::VtpEventData::n_blocks)
        .def("peak",
            [](const vtp::VtpEventData &e, int i) -> const vtp::EcPeak& {
                if (i < 0 || i >= e.n_peaks)
                    throw py::index_error("peak index out of range");
                return e.peaks[i];
            },
            py::return_value_policy::reference_internal)
        .def("cluster",
            [](const vtp::VtpEventData &e, int i) -> const vtp::EcCluster& {
                if (i < 0 || i >= e.n_clusters)
                    throw py::index_error("cluster index out of range");
                return e.clusters[i];
            },
            py::return_value_policy::reference_internal)
        .def("block",
            [](const vtp::VtpEventData &e, int i) -> const vtp::VtpBlock& {
                if (i < 0 || i >= e.n_blocks)
                    throw py::index_error("block index out of range");
                return e.blocks[i];
            },
            py::return_value_policy::reference_internal)
        .def("clear", &vtp::VtpEventData::clear);
}

// -------------------------------------------------------------------------
// TDC (tdc::)
// -------------------------------------------------------------------------
void bind_tdc(py::module_ &m)
{
    py::class_<tdc::TdcHit>(m, "TdcHit")
        .def_readonly("roc_tag", &tdc::TdcHit::roc_tag)
        .def_readonly("slot",    &tdc::TdcHit::slot)
        .def_readonly("channel", &tdc::TdcHit::channel)
        .def_readonly("edge",    &tdc::TdcHit::edge)
        .def_readonly("value",   &tdc::TdcHit::value)
        .def("__repr__", [](const tdc::TdcHit &h) {
            char buf[128];
            std::snprintf(buf, sizeof(buf),
                "<TdcHit roc=0x%x slot=%u ch=%u edge=%u value=%u>",
                h.roc_tag, h.slot, h.channel, h.edge, h.value);
            return std::string(buf);
        });

    py::class_<tdc::TdcEventData, std::shared_ptr<tdc::TdcEventData>>(m, "TdcEventData")
        .def(py::init<>())
        .def_readonly("n_hits", &tdc::TdcEventData::n_hits)
        .def("hit",
            [](const tdc::TdcEventData &e, int i) -> const tdc::TdcHit& {
                if (i < 0 || i >= e.n_hits)
                    throw py::index_error("hit index out of range");
                return e.hits[i];
            },
            py::return_value_policy::reference_internal)
        .def_property_readonly("hits_numpy",
            [](const tdc::TdcEventData &e) {
                // Bulk accessor used by tight Python loops: returns the
                // first n_hits entries of the hits[] buffer as a numpy
                // structured array (copy).  Layout matches the in-memory
                // tdc::TdcHit struct exactly — 12 bytes per row, with one
                // byte of padding between ``edge`` and ``value`` so the
                // uint32 stays 4-byte aligned.
                py::list fields;
                fields.append(py::make_tuple("roc_tag", "<u4"));
                fields.append(py::make_tuple("slot",    "u1"));
                fields.append(py::make_tuple("channel", "u1"));
                fields.append(py::make_tuple("edge",    "u1"));
                fields.append(py::make_tuple("_pad",    "u1"));
                fields.append(py::make_tuple("value",   "<u4"));
                py::dtype dt = py::dtype::from_args(fields);
                // Pass the buffer address; py::array copies when no base
                // handle is given, so the returned array is independent of
                // the TdcEventData (safe across the next decode call).
                const void *src = (e.n_hits > 0) ? (const void*)e.hits
                                                  : nullptr;
                return py::array(dt, (py::ssize_t)e.n_hits, src);
            },
            "Hits as a numpy structured array (fresh copy, length n_hits).\n"
            "Fields: roc_tag <u4, slot u1, channel u1, edge u1, _pad u1, "
            "value <u4.  Use this for bulk per-event loops to avoid "
            "per-attribute Python-to-C++ call overhead.")
        .def("clear", &tdc::TdcEventData::clear);
}

// -------------------------------------------------------------------------
// DaqConfig + helpers
// -------------------------------------------------------------------------
void bind_config(py::module_ &m)
{
    py::class_<evc::DaqConfig::RocEntry>(m, "RocEntry")
        .def_readwrite("tag",   &evc::DaqConfig::RocEntry::tag)
        .def_readwrite("name",  &evc::DaqConfig::RocEntry::name)
        .def_readwrite("crate", &evc::DaqConfig::RocEntry::crate)
        .def_readwrite("type",  &evc::DaqConfig::RocEntry::type);

    py::class_<evc::DaqConfig>(m, "DaqConfig",
        "DAQ configuration (bank tags, ROC map, TI layout, ...)")
        .def(py::init<>())
        .def_readwrite("physics_tags",      &evc::DaqConfig::physics_tags)
        .def_readwrite("physics_base",      &evc::DaqConfig::physics_base)
        .def_readwrite("monitoring_tags",   &evc::DaqConfig::monitoring_tags)
        .def_readwrite("prestart_tag",      &evc::DaqConfig::prestart_tag)
        .def_readwrite("go_tag",            &evc::DaqConfig::go_tag)
        .def_readwrite("end_tag",           &evc::DaqConfig::end_tag)
        .def_readwrite("sync_tag",          &evc::DaqConfig::sync_tag)
        .def_readwrite("epics_tag",         &evc::DaqConfig::epics_tag)
        .def_readwrite("adc_format",        &evc::DaqConfig::adc_format)
        .def_readwrite("sparsify_sigma",    &evc::DaqConfig::sparsify_sigma)
        .def_readwrite("fadc_composite_tag",&evc::DaqConfig::fadc_composite_tag)
        .def_readwrite("adc1881m_bank_tag", &evc::DaqConfig::adc1881m_bank_tag)
        .def_readwrite("ti_bank_tag",       &evc::DaqConfig::ti_bank_tag)
        .def_readwrite("trigger_bank_tag",  &evc::DaqConfig::trigger_bank_tag)
        .def_readwrite("run_info_tag",      &evc::DaqConfig::run_info_tag)
        .def_readwrite("daq_config_tag",    &evc::DaqConfig::daq_config_tag)
        .def_readwrite("epics_bank_tag",    &evc::DaqConfig::epics_bank_tag)
        .def_readwrite("ssp_bank_tags",     &evc::DaqConfig::ssp_bank_tags)
        .def_readwrite("fadc_raw_tag",      &evc::DaqConfig::fadc_raw_tag)
        .def_readwrite("tdc_bank_tag",      &evc::DaqConfig::tdc_bank_tag)
        .def_readwrite("roc_tags",          &evc::DaqConfig::roc_tags)
        .def_readwrite("ti_master_tag",     &evc::DaqConfig::ti_master_tag)
        .def_readwrite("verbose_decode",    &evc::DaqConfig::verbose_decode)
        .def_readwrite("fadc250_fw",        &evc::DaqConfig::fadc250_fw)
        .def("is_physics",    &evc::DaqConfig::is_physics)
        .def("is_monitoring", &evc::DaqConfig::is_monitoring)
        .def("is_control",    &evc::DaqConfig::is_control)
        .def("is_sync",       &evc::DaqConfig::is_sync)
        .def("is_epics",      &evc::DaqConfig::is_epics)
        .def("is_ssp_bank",   &evc::DaqConfig::is_ssp_bank);

    m.def("load_daq_config",
        [](const std::string &path) {
            evc::DaqConfig cfg;
            std::string p = path.empty() ? default_daq_config_path() : path;
            if (!evc::load_daq_config(p, cfg))
                throw std::runtime_error("Failed to load DAQ config: " + p);
            return cfg;
        },
        py::arg("path") = std::string(""),
        "Load a DaqConfig from JSON. Empty path uses the installed default.");
}

// -------------------------------------------------------------------------
// EventType / Status enums
// -------------------------------------------------------------------------
void bind_enums(py::module_ &m)
{
    py::enum_<evc::EventType>(m, "EventType")
        .value("Unknown",  evc::EventType::Unknown)
        .value("Physics",  evc::EventType::Physics)
        .value("Sync",     evc::EventType::Sync)
        .value("Epics",    evc::EventType::Epics)
        .value("Prestart", evc::EventType::Prestart)
        .value("Go",       evc::EventType::Go)
        .value("End",      evc::EventType::End)
        .value("Control",  evc::EventType::Control);

    py::enum_<evc::status>(m, "Status")
        .value("failure",    evc::status::failure)
        .value("success",    evc::status::success)
        .value("incomplete", evc::status::incomplete)
        .value("empty",      evc::status::empty)
        .value("eof",        evc::status::eof);
}

// -------------------------------------------------------------------------
// EvChannel
// -------------------------------------------------------------------------
void bind_channel(py::module_ &m)
{
    py::class_<evc::EvChannel>(m, "EvChannel",
        "Evio event reader and scanner.")
        .def(py::init<size_t>(),
            py::arg("buflen") = 1024u * 2000u,
            "Construct with an internal buffer of `buflen` uint32 words.")
        .def("set_config", &evc::EvChannel::SetConfig, py::arg("cfg"))
        .def("get_config", &evc::EvChannel::GetConfig,
             py::return_value_policy::reference_internal)
        .def("open_sequential",
            [](evc::EvChannel &self, const std::string &path) {
                py::gil_scoped_release rel;
                return self.OpenSequential(path);
            },
            py::arg("path"),
            "Open an evio file in sequential mode.  Pairs with "
            "open_random_access() and open_auto() — most callers should "
            "prefer open_auto() which picks the best mode automatically.")
        .def("close", &evc::EvChannel::Close)
        .def("read",
            [](evc::EvChannel &self) {
                py::gil_scoped_release rel;
                return self.Read();
            },
            "Read the next record into the internal buffer. Returns Status.")
        // ---- Random-access mode (evio "ra") -------------------------------
        .def("open_random_access",
            [](evc::EvChannel &self, const std::string &path) {
                py::gil_scoped_release rel;
                return self.OpenRandomAccess(path);
            },
            py::arg("path"),
            "Open an evio file in random-access mode.  evio mmaps the file "
            "and builds an event pointer table during open — after this you "
            "can jump to any event via read_event_by_index().  Use "
            "get_random_access_event_count() to get the total.")
        .def("get_random_access_event_count",
            &evc::EvChannel::GetRandomAccessEventCount,
            "Total number of events in the random-access table.  Returns 0 "
            "if not opened in random-access mode.")
        .def("read_event_by_index",
            [](evc::EvChannel &self, int i) {
                py::gil_scoped_release rel;
                return self.ReadEventByIndex(i);
            },
            py::arg("index"),
            "Read the event at 0-based evio index `i` into the internal "
            "buffer.  scan() / select_event() / info() / fadc() / ... then "
            "work identically to the sequential path.  Returns Status.")
        .def("open_auto",
            [](evc::EvChannel &self, const std::string &path) {
                py::gil_scoped_release rel;
                return self.OpenAuto(path);
            },
            py::arg("path"),
            "Open ``path`` with random-access mode if the file supports it, "
            "otherwise fall back to the sequential mode.  After success, "
            "``is_random_access()`` reports which mode was selected; callers "
            "dispatch between ``read_event_by_index()`` (RA) and ``read()`` "
            "(sequential).  Recommended default when you don't care which "
            "mode you get.")
        .def("is_random_access", &evc::EvChannel::IsRandomAccess,
            "True iff the current handle was opened in random-access mode.")
        .def("scan", &evc::EvChannel::Scan,
            "Scan the currently-held record. Call after a successful Read().")
        .def("get_event_type", &evc::EvChannel::GetEventType)
        .def("get_n_events",   &evc::EvChannel::GetNEvents)
        // ---- Lazy per-product accessors (new API) -------------------------
        // select_event() picks the sub-event; info/fadc/gem/tdc/vtp each
        // decode on first call and return the cached result on repeat calls.
        .def("select_event", &evc::EvChannel::SelectEvent, py::arg("i") = 0,
            "Select the sub-event index for subsequent info/fadc/gem/tdc/vtp "
            "calls.  Clears the product cache if the index changed.  Must be "
            "called after scan(); for PRad-II single-event data use i=0.")
        .def("info",
            [](const evc::EvChannel &self) {
                fdec::EventInfo out;
                {
                    py::gil_scoped_release rel;
                    out = self.Info();
                }
                return out;
            },
            "Decode (or return cached) event info for the currently-selected "
            "sub-event.  Cheapest accessor — skips FADC/SSP/VTP/TDC work.")
        .def("fadc",
            [](const evc::EvChannel &self) {
                auto evt = std::make_shared<fdec::EventData>();
                {
                    py::gil_scoped_release rel;
                    *evt = self.Fadc();
                }
                return evt;
            },
            "Decode (or return cached) FADC250/ADC1881M waveforms as "
            "EventData (contains event info + per-ROC/slot/channel samples).")
        .def("gem",
            [](const evc::EvChannel &self) {
                auto evt = std::make_shared<ssp::SspEventData>();
                {
                    py::gil_scoped_release rel;
                    *evt = self.Gem();
                }
                return evt;
            },
            "Decode (or return cached) SSP/MPD GEM strip data as SspEventData.")
        .def("tdc",
            [](const evc::EvChannel &self) {
                auto evt = std::make_shared<tdc::TdcEventData>();
                {
                    py::gil_scoped_release rel;
                    *evt = self.Tdc();
                }
                return evt;
            },
            "Decode (or return cached) V1190 TDC tagger hits as TdcEventData.")
        .def("vtp",
            [](const evc::EvChannel &self) {
                auto evt = std::make_shared<vtp::VtpEventData>();
                {
                    py::gil_scoped_release rel;
                    *evt = self.Vtp();
                }
                return evt;
            },
            "Decode (or return cached) VTP ECAL peaks/clusters as VtpEventData.")
        .def("decode_event",
            [](const evc::EvChannel &self, int i,
               bool with_ssp, bool with_vtp, bool with_tdc) -> py::dict {
                auto evt = std::make_shared<fdec::EventData>();
                std::shared_ptr<ssp::SspEventData> ssp_ptr;
                std::shared_ptr<vtp::VtpEventData> vtp_ptr;
                std::shared_ptr<tdc::TdcEventData> tdc_ptr;
                if (with_ssp) ssp_ptr = std::make_shared<ssp::SspEventData>();
                if (with_vtp) vtp_ptr = std::make_shared<vtp::VtpEventData>();
                if (with_tdc) tdc_ptr = std::make_shared<tdc::TdcEventData>();
                bool ok;
                {
                    py::gil_scoped_release rel;
                    ok = self.DecodeEvent(i, *evt,
                                          ssp_ptr.get(), vtp_ptr.get(), tdc_ptr.get());
                }
                py::dict out;
                out["ok"]    = ok;
                out["event"] = py::cast(evt);
                out["ssp"]   = ssp_ptr ? py::cast(ssp_ptr) : py::none();
                out["vtp"]   = vtp_ptr ? py::cast(vtp_ptr) : py::none();
                out["tdc"]   = tdc_ptr ? py::cast(tdc_ptr) : py::none();
                return out;
            },
            py::arg("i") = 0,
            py::kw_only(),
            py::arg("with_ssp") = false,
            py::arg("with_vtp") = false,
            py::arg("with_tdc") = false,
            "Full decode. Returns {'ok': bool, 'event': EventData, "
            "'ssp': SspEventData|None, 'vtp': ..., 'tdc': ...}.")
        .def("sync",
            [](const evc::EvChannel &self) {
                psync::SyncInfo out;
                {
                    py::gil_scoped_release rel;
                    out = self.Sync();
                }
                return out;
            },
            "Absolute-time / run-state snapshot.  Persists across events — "
            "refreshed only when the channel lands on a SYNC/EPICS or "
            "control event (PRESTART/GO/END), otherwise returns the prior "
            "snapshot.  Diff `sync_counter` against your last-seen value to "
            "detect new SYNC ticks; for control events use `event_tag`.")
        .def("extract_epics_text", &evc::EvChannel::ExtractEpicsText,
            "Raw EPICS payload for the current event (empty if not EPICS).");
}

} // anonymous namespace

// -------------------------------------------------------------------------
// Entry point for the main module (prad2py.cpp calls this).
// -------------------------------------------------------------------------
void register_dec(py::module_ &m)
{
    auto dec = m.def_submodule("dec",
        "prad2dec bindings: evio reader, event data types, TDC/SSP/VTP accessors.");

    bind_enums(dec);
    bind_config(dec);
    bind_fadc(dec);
    bind_ssp(dec);
    bind_vtp(dec);
    bind_tdc(dec);
    bind_channel(dec);
}
