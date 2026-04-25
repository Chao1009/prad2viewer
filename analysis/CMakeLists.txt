# analysis/ — replay and physics analysis tools (requires ROOT)

find_package(ROOT 6.0 REQUIRED CONFIG)
include(${ROOT_USE_FILE})
message(STATUS "ROOT: ${ROOT_VERSION} at ${ROOT_DIR}")

include(GNUInstallDirs)

# ---------------------------------------------------------------------------
# Shared analysis library — Replay, PhysicsTools, MatchingTools.
#
# Was previously a list of .cpp files inlined into every executable's
# source list (each tool got its own copy of every shared symbol).  Now a
# real static library so:
#   - link is faster (one compile, many links)
#   - installed-mode ACLiC scripts can call analysis::* without
#     hitting unresolved-symbol errors at .x time.
# ---------------------------------------------------------------------------
add_library(prad2ana STATIC
    src/Replay.cpp
    src/PhysicsTools.cpp
    src/MatchingTools.cpp)

# PIC so this archive can be pulled into ACLiC-built shared objects.
set_target_properties(prad2ana PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_include_directories(prad2ana
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include/prad2ana>)

target_link_libraries(prad2ana PUBLIC
    prad2dec prad2det nlohmann_json::nlohmann_json
    ROOT::Core ROOT::Tree ROOT::RIO ROOT::Hist ROOT::Graf ROOT::Gpad)

set(ANA_DEPS prad2ana)   # tools just link the lib; transitives flow through

# ---------------------------------------------------------------------------
# helper: add an analysis tool executable.  The CMake target name stays
# clean for internal referencing; the built/installed binary gets a
# `prad2ana_` prefix so multiple projects can coexist in the same bindir.
# ---------------------------------------------------------------------------
function(add_analysis_tool name)
    add_executable(${name} ${ARGN})
    target_compile_definitions(${name} PRIVATE DATABASE_DIR="${DATABASE_DIR}")
    target_link_libraries(${name} PRIVATE ${ANA_DEPS})
    set_target_properties(${name} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
        OUTPUT_NAME "prad2ana_${name}")
endfunction()

add_analysis_tool(replay_rawdata tools/replay_rawdata.cpp)
add_analysis_tool(replay_recon tools/replay_recon.cpp)
add_analysis_tool(epCalib   tools/epCalib.cpp)
add_analysis_tool(cosmic_test tools/cosmic_test.cpp)
add_analysis_tool(quick_check tools/quick_check.cpp)
add_analysis_tool(matching tools/matching.cpp)
add_analysis_tool(det_calib tools/det_calib.cpp)
add_analysis_tool(gain_monitor tools/gain_monitor.cpp)
add_analysis_tool(gain_fitter tools/gain_fitter.cpp)
add_analysis_tool(sim2replay tools/sim2replay.cpp)
add_analysis_tool(cosmic_lg tools/cosmic_lg.cpp)
add_analysis_tool(gain_replay tools/gain_replay.cpp)

install(TARGETS replay_rawdata replay_recon epCalib
                cosmic_test quick_check matching det_calib cosmic_lg sim2replay gain_replay
                gain_monitor gain_fitter
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

# Install the static analysis library so installed-mode ACLiC scripts can
# link against analysis::* symbols (no shared lib needed since tools also
# link statically).
install(TARGETS prad2ana ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

# Install analysis/include headers so installed-tree ACLiC scripts can
# include "PhysicsTools.h" / "MatchingTools.h" / "Replay.h" / "ConfigSetup.h"
# the same way they do from the build tree.
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/prad2ana
        FILES_MATCHING PATTERN "*.h")

# Install the ACLiC macro scripts (rootlogon.C, gem_*, lms_*, tagger_*) so
# the install-tree workflow described in rootlogon.C actually has files to
# load: `root -l <prefix>/share/prad2evviewer/analysis/scripts/rootlogon.C`.
install(DIRECTORY scripts/ DESTINATION ${CMAKE_INSTALL_DATADIR}/prad2evviewer/analysis/scripts
        FILES_MATCHING PATTERN "*.C")
