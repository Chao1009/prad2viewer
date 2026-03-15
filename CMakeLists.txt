cmake_minimum_required(VERSION 3.14)
project(evc VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# ---------------------------------------------------------------------------
# Paths to prebuilt CODA/ET libraries and headers
# Override with -DCODA_LIB_DIR=... and -DCODA_INCLUDE_DIR=... if needed
# ---------------------------------------------------------------------------
set(CODA_ROOT "/usr/clas12/release/1.4.0/coda/Linux_x86_64_RHEL9"
    CACHE PATH "Root of the CODA installation")

set(CODA_LIB_DIR     "${CODA_ROOT}/lib"
    CACHE PATH "Directory containing libevio.a and libet.a")

set(CODA_INCLUDE_DIR "/usr/clas12/release/1.4.0/coda/common/include"
    CACHE PATH "Directory containing evio.h and et.h")

# ---------------------------------------------------------------------------
# Locate libevio.a
# ---------------------------------------------------------------------------
find_library(EVIO_LIB
    NAMES evio
    PATHS ${CODA_LIB_DIR}
    NO_DEFAULT_PATH
    REQUIRED
)
message(STATUS "Found libevio : ${EVIO_LIB}")

add_library(evio STATIC IMPORTED GLOBAL)
set_target_properties(evio PROPERTIES
    IMPORTED_LOCATION             "${EVIO_LIB}"
    INTERFACE_INCLUDE_DIRECTORIES "${CODA_INCLUDE_DIR}"
)

# ---------------------------------------------------------------------------
# Locate libet.a
# ---------------------------------------------------------------------------
find_library(ET_LIB
    NAMES et
    PATHS ${CODA_LIB_DIR}
    NO_DEFAULT_PATH
    REQUIRED
)
message(STATUS "Found libet   : ${ET_LIB}")

add_library(et STATIC IMPORTED GLOBAL)
set_target_properties(et PROPERTIES
    IMPORTED_LOCATION             "${ET_LIB}"
    INTERFACE_INCLUDE_DIRECTORIES "${CODA_INCLUDE_DIR}"
)

# ---------------------------------------------------------------------------
# Static library: libevc.a
# ---------------------------------------------------------------------------
add_library(evc STATIC
    src/EvChannel.cpp
    src/EtChannel.cpp
)

target_include_directories(evc
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

target_link_libraries(evc
    PUBLIC
        evio
        et
)

# ---------------------------------------------------------------------------
# Test executable
# ---------------------------------------------------------------------------
add_executable(evc_test
    test/test_main.cpp
)

target_link_libraries(evc_test
    PRIVATE
        evc
)

# ---------------------------------------------------------------------------
# Install
# ---------------------------------------------------------------------------
include(GNUInstallDirs)

install(TARGETS evc
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(FILES
    include/EvStruct.h
    include/EvChannel.h
    include/EtChannel.h
    include/EtConfigWrapper.h
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/evc
)

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
message(STATUS "----------------------------------------------------")
message(STATUS "  evc ${PROJECT_VERSION}")
message(STATUS "  Build type   : ${CMAKE_BUILD_TYPE}")
message(STATUS "  CODA lib dir : ${CODA_LIB_DIR}")
message(STATUS "  CODA inc dir : ${CODA_INCLUDE_DIR}")
message(STATUS "----------------------------------------------------")
