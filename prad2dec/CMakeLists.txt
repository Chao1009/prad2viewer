# prad2dec/CMakeLists.txt — static library: libprad2dec.a
#
# Can be used standalone:  add_subdirectory(prad2dec)
# Or as part of the top-level project.

cmake_minimum_required(VERSION 3.14)

# If built standalone, set up project and standards
if(NOT CMAKE_PROJECT_NAME)
    project(prad2dec VERSION 1.0.0 LANGUAGES CXX)
    set(CMAKE_CXX_STANDARD 17)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif()

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

include(FetchContent)
set(FETCHCONTENT_QUIET OFF CACHE BOOL "" FORCE)

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
    set(ET_INCLUDE_DIR "${et_SOURCE_DIR}/src/libsrc" CACHE INTERNAL "")
elseif(ET_SOURCE STREQUAL "system")
    find_package(coda_et REQUIRED)
else()
    find_library(ET_LIB NAMES et PATHS ${CODA_LIB_DIR} NO_DEFAULT_PATH REQUIRED)
    add_library(et STATIC IMPORTED GLOBAL)
    set_target_properties(et PROPERTIES
        IMPORTED_LOCATION             "${ET_LIB}"
        INTERFACE_INCLUDE_DIRECTORIES "${CODA_INCLUDE_DIR}"
    )
    set(ET_INCLUDE_DIR "${CODA_INCLUDE_DIR}" CACHE INTERNAL "")
endif()

# ---------------------------------------------------------------------------
# Library target
# ---------------------------------------------------------------------------
add_library(prad2dec STATIC
    src/EvChannel.cpp
    src/EtChannel.cpp
    src/Fadc250Decoder.cpp
    src/Adc1881mDecoder.cpp
    src/WaveAnalyzer.cpp
    src/replay.cpp
)

target_include_directories(prad2dec
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include/prad2dec>
        $<BUILD_INTERFACE:${ET_INCLUDE_DIR}>
)

target_link_libraries(prad2dec PUBLIC evio et)

# ---------------------------------------------------------------------------
# Install
# ---------------------------------------------------------------------------
include(GNUInstallDirs)

install(TARGETS prad2dec ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/prad2dec)
