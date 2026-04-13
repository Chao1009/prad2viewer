# analysis/ — replay and physics analysis tools (requires ROOT)

find_package(ROOT 6.0 REQUIRED CONFIG)
include(${ROOT_USE_FILE})
message(STATUS "ROOT: ${ROOT_VERSION} at ${ROOT_DIR}")

# shared sources + deps for all analysis tools
set(ANA_SOURCES src/Replay.cpp src/PhysicsTools.cpp src/MatchingTools.cpp)
set(ANA_DEPS prad2dec prad2det nlohmann_json::nlohmann_json
    ROOT::Core ROOT::Tree ROOT::RIO ROOT::Hist ROOT::Graf ROOT::Gpad)

# helper: add an analysis tool executable
function(add_analysis_tool name)
    add_executable(${name} ${ARGN} ${ANA_SOURCES})
    target_include_directories(${name} PRIVATE include)
    target_compile_definitions(${name} PRIVATE DATABASE_DIR="${DATABASE_DIR}")
    target_link_libraries(${name} PRIVATE ${ANA_DEPS})
    set_target_properties(${name} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
endfunction()

add_analysis_tool(replay_rawdata tools/replay_rawdata.cpp)
add_analysis_tool(replay_rawdata_m tools/replay_rawdata_m.cpp)
add_analysis_tool(replay_recon tools/replay_recon.cpp)
add_analysis_tool(replay_recon_m tools/replay_recon_m.cpp)
add_analysis_tool(epCalib tools/epCalib.cpp)
add_analysis_tool(cosmic_test tools/cosmic_test.cpp)
add_analysis_tool(quick_check tools/quick_check.cpp)
add_analysis_tool(matching tools/matching.cpp)
add_analysis_tool(det_calib tools/det_calib.cpp)
add_analysis_tool(cosmic_lg tools/cosmic_lg.cpp)
add_analysis_tool(sim2replay tools/sim2replay.cpp)

install(TARGETS replay_rawdata replay_rawdata_m replay_recon replay_recon_m epCalib
                cosmic_test quick_check matching det_calib cosmic_lg sim2replay
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
