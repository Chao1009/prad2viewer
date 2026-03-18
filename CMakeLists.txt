cmake_minimum_required(VERSION 3.14)
project(prad2dec VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# ---------------------------------------------------------------------------
# Database and resource directories (auto-discovered, overridable)
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
# Dependencies: evio and et
# ---------------------------------------------------------------------------
set(EVIO_SOURCE "fetch" CACHE STRING "evio provider: prebuilt | fetch | system")
set(ET_SOURCE   "fetch" CACHE STRING "et provider:   prebuilt | fetch | system")

set(CODA_ROOT "/usr/clas12/release/2.0.0/coda/Linux_x86_64_RHEL9"
    CACHE PATH "Root of the prebuilt CODA installation")
set(CODA_LIB_DIR     "${CODA_ROOT}/lib"     CACHE PATH "Prebuilt lib dir")
set(CODA_INCLUDE_DIR "/usr/clas12/release/2.0.0/coda/common/include"
    CACHE PATH "CODA include dir")

set(FETCHCONTENT_QUIET OFF CACHE BOOL "" FORCE)
include(FetchContent)

# --- evio -----------------------------------------------------------------
if(EVIO_SOURCE STREQUAL "fetch")
    set(MAKE_EXAMPLES OFF CACHE BOOL "" FORCE)
    FetchContent_Declare(evio
        GIT_REPOSITORY https://github.com/JeffersonLab/evio.git
        GIT_TAG        evio-6.0
        GIT_PROGRESS   TRUE
    )
    FetchContent_MakeAvailable(evio)
elseif(EVIO_SOURCE STREQUAL "system")
    find_package(evio REQUIRED)
else()
    find_library(EVIO_LIB NAMES evio PATHS ${CODA_LIB_DIR} NO_DEFAULT_PATH REQUIRED)
    find_library(EXPAT_LIB NAMES expat REQUIRED)
    add_library(evio STATIC IMPORTED GLOBAL)
    set_target_properties(evio PROPERTIES
        IMPORTED_LOCATION             "${EVIO_LIB}"
        INTERFACE_INCLUDE_DIRECTORIES "${CODA_INCLUDE_DIR}"
        INTERFACE_LINK_LIBRARIES      "${EXPAT_LIB}"
    )
endif()

# --- et -------------------------------------------------------------------
if(ET_SOURCE STREQUAL "fetch")
    FetchContent_Declare(et
        GIT_REPOSITORY https://github.com/JeffersonLab/et.git
        GIT_TAG        main
        GIT_PROGRESS   TRUE
    )
    FetchContent_MakeAvailable(et)
    set(ET_INCLUDE_DIR "${et_SOURCE_DIR}/src/libsrc")
elseif(ET_SOURCE STREQUAL "system")
    find_package(coda_et REQUIRED)
else()
    find_library(ET_LIB NAMES et PATHS ${CODA_LIB_DIR} NO_DEFAULT_PATH REQUIRED)
    add_library(et STATIC IMPORTED GLOBAL)
    set_target_properties(et PROPERTIES
        IMPORTED_LOCATION             "${ET_LIB}"
        INTERFACE_INCLUDE_DIRECTORIES "${CODA_INCLUDE_DIR}"
    )
    set(ET_INCLUDE_DIR "${CODA_INCLUDE_DIR}")
endif()

# --- nlohmann/json --------------------------------------------------------
FetchContent_Declare(json
    GIT_REPOSITORY https://github.com/nlohmann/json.git
    GIT_TAG        v3.11.3
    GIT_SHALLOW    TRUE
)

# --- websocketpp ----------------------------------------------------------
FetchContent_Declare(websocketpp
    GIT_REPOSITORY https://github.com/zaphoyd/websocketpp.git
    GIT_TAG        0.8.2
    GIT_SHALLOW    TRUE
)

# --- standalone asio ------------------------------------------------------
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
# Static library: libprad2dec.a  (decoder + waveform analysis)
# ---------------------------------------------------------------------------
add_library(prad2dec STATIC
    prad2dec/src/EvChannel.cpp
    prad2dec/src/EtChannel.cpp
    prad2dec/src/Fadc250Decoder.cpp
    prad2dec/src/WaveAnalyzer.cpp
)

target_include_directories(prad2dec
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/prad2dec/include>
        $<INSTALL_INTERFACE:include/prad2dec>
        $<BUILD_INTERFACE:${ET_INCLUDE_DIR}>
)

target_link_libraries(prad2dec PUBLIC evio et)

# ---------------------------------------------------------------------------
# Test executable
# ---------------------------------------------------------------------------
add_executable(evc_test test/test_main.cpp)
target_link_libraries(evc_test PRIVATE prad2dec)
set_target_properties(evc_test PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# ---------------------------------------------------------------------------
# Event viewer (HTTP server)
# ---------------------------------------------------------------------------
add_executable(evc_viewer src/evc_viewer.cpp)

target_compile_definitions(evc_viewer PRIVATE
    DATABASE_DIR="${DATABASE_DIR}"
    RESOURCE_DIR="${RESOURCE_DIR}"
)

target_link_libraries(evc_viewer PRIVATE
    prad2dec
    nlohmann_json::nlohmann_json
    websocketpp_lib
    pthread
)
set_target_properties(evc_viewer PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# ---------------------------------------------------------------------------
# Online monitor (ET reader + WebSocket push)
# ---------------------------------------------------------------------------
add_executable(evc_monitor src/evc_monitor.cpp)

target_compile_definitions(evc_monitor PRIVATE
    DATABASE_DIR="${DATABASE_DIR}"
    RESOURCE_DIR="${RESOURCE_DIR}"
)

target_link_libraries(evc_monitor PRIVATE
    prad2dec
    nlohmann_json::nlohmann_json
    websocketpp_lib
    pthread
)
set_target_properties(evc_monitor PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# ---------------------------------------------------------------------------
# Install
# ---------------------------------------------------------------------------
include(GNUInstallDirs)

install(TARGETS prad2dec ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY prad2dec/include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/prad2dec)

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
message(STATUS "----------------------------------------------------")
message(STATUS "  ${PROJECT_NAME} ${PROJECT_VERSION}")
message(STATUS "  Build type   : ${CMAKE_BUILD_TYPE}")
message(STATUS "  evio source  : ${EVIO_SOURCE}")
message(STATUS "  et source    : ${ET_SOURCE}")
message(STATUS "  Database     : ${DATABASE_DIR}")
message(STATUS "  Resources    : ${RESOURCE_DIR}")
message(STATUS "----------------------------------------------------")
