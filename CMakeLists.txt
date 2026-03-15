cmake_minimum_required(VERSION 3.14)
project(evc VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# ---------------------------------------------------------------------------
# Dependencies: evio and et
# Three modes, selected by -DEVIO_SOURCE=... and -DET_SOURCE=...
#   prebuilt  (default) - use .a files from CODA_LIB_DIR
#   fetch               - download and build from GitHub with FetchContent
#   system              - let find_package locate an installed version
# ---------------------------------------------------------------------------
set(EVIO_SOURCE "fetch" CACHE STRING "evio provider: prebuilt | fetch | system")
set(ET_SOURCE   "fetch" CACHE STRING "et provider:   prebuilt | fetch | system")

set(CODA_ROOT "/usr/clas12/release/1.4.0/coda/Linux_x86_64_RHEL9"
    CACHE PATH "Root of the prebuilt CODA installation")
set(CODA_LIB_DIR     "${CODA_ROOT}/lib"
    CACHE PATH "Directory containing prebuilt libevio.a and libet.a")
set(CODA_INCLUDE_DIR "/usr/clas12/release/1.4.0/coda/common/include"
    CACHE PATH "Directory containing evio.h and et.h")

# --- evio -----------------------------------------------------------------
if(EVIO_SOURCE STREQUAL "fetch")
    include(FetchContent)
    FetchContent_Declare(evio
        GIT_REPOSITORY https://github.com/JeffersonLab/evio.git
        GIT_TAG        evio-6.0
    )
    FetchContent_MakeAvailable(evio)
    message(STATUS "evio: building from source (FetchContent)")

elseif(EVIO_SOURCE STREQUAL "system")
    find_package(evio REQUIRED)
    message(STATUS "evio: using system installation")

else()  # prebuilt
    find_library(EVIO_LIB NAMES evio PATHS ${CODA_LIB_DIR} NO_DEFAULT_PATH REQUIRED)
    message(STATUS "evio: prebuilt at ${EVIO_LIB}")

    find_library(EXPAT_LIB NAMES expat REQUIRED)
    message(STATUS "Found libexpat: ${EXPAT_LIB}")

    add_library(evio STATIC IMPORTED GLOBAL)
    set_target_properties(evio PROPERTIES
        IMPORTED_LOCATION             "${EVIO_LIB}"
        INTERFACE_INCLUDE_DIRECTORIES "${CODA_INCLUDE_DIR}"
        INTERFACE_LINK_LIBRARIES      "${EXPAT_LIB}"
    )
endif()

# --- et -------------------------------------------------------------------
if(ET_SOURCE STREQUAL "fetch")
    include(FetchContent)
    FetchContent_Declare(et
        GIT_REPOSITORY https://github.com/JeffersonLab/et.git
        GIT_TAG        main
    )
    FetchContent_MakeAvailable(et)
    message(STATUS "et: building from source (FetchContent)")

elseif(ET_SOURCE STREQUAL "system")
    find_package(coda_et REQUIRED)
    message(STATUS "et: using system installation")

else()  # prebuilt
    find_library(ET_LIB NAMES et PATHS ${CODA_LIB_DIR} NO_DEFAULT_PATH REQUIRED)
    message(STATUS "et: prebuilt at ${ET_LIB}")

    add_library(et STATIC IMPORTED GLOBAL)
    set_target_properties(et PROPERTIES
        IMPORTED_LOCATION             "${ET_LIB}"
        INTERFACE_INCLUDE_DIRECTORIES "${CODA_INCLUDE_DIR}"
    )
endif()

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
message(STATUS "  evio source  : ${EVIO_SOURCE}")
message(STATUS "  et source    : ${ET_SOURCE}")
message(STATUS "----------------------------------------------------")
