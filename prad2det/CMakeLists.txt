# prad2det/CMakeLists.txt — static library: libprad2det.a
#
# Analysis library: HyCal clustering, event reconstruction, hit matching.
# Depends on nlohmann/json for configuration loading.

cmake_minimum_required(VERSION 3.14)

# If built standalone, set up project and standards
if(NOT CMAKE_PROJECT_NAME)
    project(prad2det VERSION 1.0.0 LANGUAGES CXX)
    set(CMAKE_CXX_STANDARD 17)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif()

# ---------------------------------------------------------------------------
# Dependencies
# ---------------------------------------------------------------------------
include(FetchContent)
set(FETCHCONTENT_QUIET OFF CACHE BOOL "" FORCE)

# nlohmann/json (header-only) — reuses top-level fetch if already available
if(NOT TARGET nlohmann_json::nlohmann_json)
    FetchContent_Declare(json
        GIT_REPOSITORY https://github.com/nlohmann/json.git
        GIT_TAG        v3.11.3
        GIT_SHALLOW    TRUE
    )
    FetchContent_MakeAvailable(json)
endif()

# ---------------------------------------------------------------------------
# Library target
# ---------------------------------------------------------------------------
add_library(prad2det STATIC
    src/HyCalSystem.cpp
    src/HyCalCluster.cpp
    src/EpicsStore.cpp
    src/GemSystem.cpp
    src/GemCluster.cpp
    src/GemPedestal.cpp
    src/PipelineBuilder.cpp
    src/TrackGeometry.cpp
    src/TrackMatcher.cpp
)

# PIC so the archive can be pulled into shared objects (pybind11 module,
# ACLiC-compiled ROOT macros).  See prad2dec/CMakeLists.txt for details.
set_target_properties(prad2det PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_include_directories(prad2det
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include/prad2det>
)

target_link_libraries(prad2det PUBLIC nlohmann_json::nlohmann_json)

# PipelineBuilder calls evc::load_daq_config (declared inline in
# prad2dec/include/load_daq_config.h), so prad2det now needs prad2dec at
# link time too — promoted from PRIVATE to PUBLIC so any caller of
# prad2det that uses the builder also gets prad2dec transitively.
# (GemSystem also needs SspData.h from prad2dec; this covers it.)
if(TARGET prad2dec)
    target_link_libraries(prad2det PUBLIC prad2dec)
endif()

# ---------------------------------------------------------------------------
# Install
# ---------------------------------------------------------------------------
include(GNUInstallDirs)

install(TARGETS prad2det ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/prad2det)

# ---------------------------------------------------------------------------
# Optional tests
# ---------------------------------------------------------------------------
option(PRAD2DET_BUILD_TESTS "Build prad2det unit tests" OFF)
if(PRAD2DET_BUILD_TESTS)
    enable_testing()
    add_executable(prad2det_test_track_matcher test/test_track_matcher.cpp)
    target_link_libraries(prad2det_test_track_matcher PRIVATE prad2det)
    add_test(NAME prad2det.track_matcher COMMAND prad2det_test_track_matcher)
endif()
