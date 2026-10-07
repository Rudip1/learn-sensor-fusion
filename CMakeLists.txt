cmake_minimum_required(VERSION 3.20)
project(sensor_fusion VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# A Python build (pip install .) only needs the library and the extension module.
if(DEFINED SKBUILD)
    set(_sf_default_dev OFF)
    set(_sf_default_python ON)
else()
    set(_sf_default_dev ON)
    set(_sf_default_python OFF)
endif()
option(SF_BUILD_TESTS "Build the Catch2 test suite" ${_sf_default_dev})
option(SF_BUILD_EXAMPLES "Build the command-line examples" ${_sf_default_dev})
option(SF_BUILD_PYTHON "Build the pybind11 module" ${_sf_default_python})

include(FetchContent)

# Eigen 3.4: use the system package when present, otherwise fetch the headers.
find_package(Eigen3 3.4 QUIET NO_MODULE)
if(NOT TARGET Eigen3::Eigen)
    message(STATUS "Eigen 3.4 not found on the system; fetching it")
    FetchContent_Declare(
        eigen
        URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
        URL_HASH SHA256=8586084f71f9bde545ee7fa6d00288b264a2b7ac3607b974e54d13e7162c1c72
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE
        SOURCE_SUBDIR headers-only)  # no CMakeLists there: download without configuring Eigen's build
    FetchContent_MakeAvailable(eigen)
    add_library(Eigen3::Eigen INTERFACE IMPORTED)
    set_target_properties(Eigen3::Eigen PROPERTIES INTERFACE_INCLUDE_DIRECTORIES ${eigen_SOURCE_DIR})
endif()

add_subdirectory(cpp)

if(SF_BUILD_TESTS)
    enable_testing()
    add_subdirectory(cpp/tests)
endif()

if(SF_BUILD_EXAMPLES)
    add_subdirectory(cpp/examples)
endif()

if(SF_BUILD_PYTHON)
    add_subdirectory(python)
endif()
