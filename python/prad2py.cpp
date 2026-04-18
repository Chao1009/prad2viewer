// prad2py.cpp — main module glue.
//
// The actual binding code lives in per-area translation units:
//   bind_dec.cpp   → prad2py.dec (evio reader, event data, TDC/SSP/VTP)
//   (future)       → prad2py.det (HyCal / GEM reconstruction)
//
// Each of those files defines a ``register_XXX(py::module_ &m)`` entry
// point that adds a submodule to the top-level module.
//
// No "do everything" helpers at module root — analyses should drive the
// per-event loop themselves via ``dec.EvChannel.select_event()`` plus
// ``info()``/``fadc()``/``gem()``/``tdc()``/``vtp()`` and accumulate into
// numpy / Python on their own terms.

#include <pybind11/pybind11.h>

#include <cstdlib>
#include <string>

namespace py = pybind11;

#ifndef DATABASE_DIR
#define DATABASE_DIR "."
#endif

namespace {

std::string default_daq_config_path()
{
    const char *env = std::getenv("PRAD2_DATABASE_DIR");
    std::string dir = env ? env : DATABASE_DIR;
    return dir + "/daq_config.json";
}

} // anonymous namespace

// Defined in bind_dec.cpp — registers the ``prad2py.dec`` submodule.
void register_dec(py::module_ &m);

PYBIND11_MODULE(prad2py, m)
{
    m.doc() = "PRad-II (prad2dec + prad2det) Python bindings.";

    m.attr("__version__")    = "0.3.0";
    m.attr("DATABASE_DIR")   = DATABASE_DIR;

    m.def("default_daq_config", &default_daq_config_path,
          "Return the default daq_config.json path used by analyses.");

    // Per-area submodules.  Phase 1: decoder.  Phase 2+: detector.
    register_dec(m);
}
