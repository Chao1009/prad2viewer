cmake_minimum_required(VERSION 3.14)
project(prad2evviewer VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# ---------------------------------------------------------------------------
# Compile-time paths (auto-discovered, overridable)
# ---------------------------------------------------------------------------
if(NOT DEFINED DATABASE_DIR)
    set(DATABASE_DIR "${CMAKE_SOURCE_DIR}/database")
endif()
get_filename_component(DATABASE_DIR "${DATABASE_DIR}" ABSOLUTE)

if(NOT DEFINED RESOURCE_DIR)
    set(RESOURCE_DIR "${CMAKE_SOURCE_DIR}/resources")
endif()
get_filename_component(RESOURCE_DIR "${RESOURCE_DIR}" ABSOLUTE)

# ---------------------------------------------------------------------------
# Decoder library (evio/et deps handled inside)
# ---------------------------------------------------------------------------
add_subdirectory(prad2dec)

# ---------------------------------------------------------------------------
# Analysis library (clustering, reconstruction)
# ---------------------------------------------------------------------------
add_subdirectory(prad2ana)

# ---------------------------------------------------------------------------
# Web server dependencies (json, websocketpp, asio)
# ---------------------------------------------------------------------------
include(FetchContent)
set(FETCHCONTENT_QUIET OFF CACHE BOOL "" FORCE)

FetchContent_Declare(json
    GIT_REPOSITORY https://github.com/nlohmann/json.git
    GIT_TAG        v3.11.3
    GIT_SHALLOW    TRUE
)
FetchContent_Declare(websocketpp
    GIT_REPOSITORY https://github.com/zaphoyd/websocketpp.git
    GIT_TAG        0.8.2
    GIT_SHALLOW    TRUE
)
FetchContent_Declare(asio
    GIT_REPOSITORY https://github.com/chriskohlhoff/asio.git
    GIT_TAG        asio-1-30-2
    GIT_SHALLOW    TRUE
)
FetchContent_MakeAvailable(json websocketpp asio)

if(NOT TARGET asio_standalone)
    add_library(asio_standalone INTERFACE)
    target_include_directories(asio_standalone INTERFACE ${asio_SOURCE_DIR}/asio/include)
    target_compile_definitions(asio_standalone INTERFACE ASIO_STANDALONE)
endif()

if(NOT TARGET websocketpp_lib)
    add_library(websocketpp_lib INTERFACE)
    target_include_directories(websocketpp_lib INTERFACE ${websocketpp_SOURCE_DIR})
    target_link_libraries(websocketpp_lib INTERFACE asio_standalone)
endif()

# ---------------------------------------------------------------------------
# Test tools
# ---------------------------------------------------------------------------
add_executable(evc_test test/test_main.cpp)
target_link_libraries(evc_test PRIVATE prad2dec)
set_target_properties(evc_test PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_executable(et_feeder test/et_feeder.cpp)
target_link_libraries(et_feeder PRIVATE prad2dec)
set_target_properties(et_feeder PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_executable(evchan_test test/evchan_test.cpp)
target_link_libraries(evchan_test PRIVATE prad2dec)
set_target_properties(evchan_test PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_executable(evio_dump test/evio_dump.cpp)
target_link_libraries(evio_dump PRIVATE prad2dec nlohmann_json::nlohmann_json)
set_target_properties(evio_dump PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_executable(ped_calc test/ped_calc.cpp)
target_link_libraries(ped_calc PRIVATE prad2dec nlohmann_json::nlohmann_json)
set_target_properties(ped_calc PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# ---------------------------------------------------------------------------
# Event viewer
# ---------------------------------------------------------------------------
add_executable(evc_viewer src/evc_viewer.cpp src/app_state.cpp)
target_compile_definitions(evc_viewer PRIVATE
    DATABASE_DIR="${DATABASE_DIR}"
    RESOURCE_DIR="${RESOURCE_DIR}"
)
target_link_libraries(evc_viewer PRIVATE
    prad2dec
    prad2ana
    nlohmann_json::nlohmann_json
    websocketpp_lib
    pthread
)
set_target_properties(evc_viewer PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# ---------------------------------------------------------------------------
# Online monitor
# ---------------------------------------------------------------------------
add_executable(evc_monitor src/evc_monitor.cpp src/app_state.cpp)
target_compile_definitions(evc_monitor PRIVATE
    DATABASE_DIR="${DATABASE_DIR}"
    RESOURCE_DIR="${RESOURCE_DIR}"
)
target_link_libraries(evc_monitor PRIVATE
    prad2dec
    prad2ana
    nlohmann_json::nlohmann_json
    websocketpp_lib
    pthread
)
set_target_properties(evc_monitor PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# ---------------------------------------------------------------------------
# Analysis tools (optional, requires ROOT)
# ---------------------------------------------------------------------------
option(BUILD_ANALYSIS "Build analysis tools (requires ROOT)" OFF)
if(BUILD_ANALYSIS)
    add_subdirectory(analysis)
endif()

# ---------------------------------------------------------------------------
# Qt WebEngine monitor client (optional)
# ---------------------------------------------------------------------------
option(BUILD_GUI "Build the Qt5 WebEngine monitor client" OFF)
if(BUILD_GUI)
    add_subdirectory(prad2qtmon)
endif()

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
message(STATUS "----------------------------------------------------")
message(STATUS "  ${PROJECT_NAME} ${PROJECT_VERSION}")
message(STATUS "  Build type   : ${CMAKE_BUILD_TYPE}")
message(STATUS "  Database     : ${DATABASE_DIR}")
message(STATUS "  Resources    : ${RESOURCE_DIR}")
message(STATUS "  Analysis     : ${BUILD_ANALYSIS}")
message(STATUS "  Qt client    : ${BUILD_GUI}")
message(STATUS "----------------------------------------------------")
