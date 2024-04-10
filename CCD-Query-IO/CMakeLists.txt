# Detects whether this is a top-level project
get_directory_property(HAS_PARENT PARENT_DIRECTORY)
if(HAS_PARENT)
  set(CCD_IO_TOPLEVEL_PROJECT OFF)
else()
  set(CCD_IO_TOPLEVEL_PROJECT ON)
endif()

# Check required CMake version
set(REQUIRED_CMAKE_VERSION "3.18.0")
if(CCD_IO_TOPLEVEL_PROJECT)
  cmake_minimum_required(VERSION ${REQUIRED_CMAKE_VERSION})
else()
  # Don't use cmake_minimum_required here to avoid implicitly overriding parent policies
  if(${CMAKE_VERSION} VERSION_LESS ${REQUIRED_CMAKE_VERSION})
    message(FATAL_ERROR "CMake required version to build IPC Toolkit is ${REQUIRED_CMAKE_VERSION}")
  endif()
endif()

# Include user-provided default options if available. We do that before the main
# `project()` so that we can define the C/C++ compilers from the option file.
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/CCDQueryIOOptions.cmake)
  message(STATUS "Using local options file: ${CMAKE_CURRENT_SOURCE_DIR}/CCDQueryIOOptions.cmake")
  include(${CMAKE_CURRENT_SOURCE_DIR}/CCDQueryIOOptions.cmake)
endif()

# Enable ccache if available
find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
  option(CCD_IO_WITH_CCACHE "Enable ccache when building CCDQueryIO" ${CCD_IO_TOPLEVEL_PROJECT})
else()
  option(CCD_IO_WITH_CCACHE "Enable ccache when building CCDQueryIO" OFF)
endif()
if(CCD_IO_WITH_CCACHE AND CCACHE_PROGRAM)
  message(STATUS "Enabling Ccache support (${CCACHE_PROGRAM})")
  set(ccacheEnv
    CCACHE_BASEDIR=${CMAKE_BINARY_DIR}
    CCACHE_SLOPPINESS=clang_index_store,include_file_ctime,include_file_mtime,locale,pch_defines,time_macros
  )
  foreach(lang IN ITEMS C CXX)
    set(CMAKE_${lang}_COMPILER_LAUNCHER
      ${CMAKE_COMMAND} -E env ${ccacheEnv} ${CCACHE_PROGRAM}
    )
  endforeach()
endif()

################################################################################
# CMake Policies
################################################################################

cmake_policy(SET CMP0054 NEW) # Only interpret if() arguments as variables or keywords when unquoted.
cmake_policy(SET CMP0076 NEW) # target_sources() command converts relative paths to absolute.
if(CMAKE_VERSION VERSION_GREATER_EQUAL "3.24")
  cmake_policy(SET CMP0135 NEW) # Set the timestamps of all extracted contents to the time of the extraction.
endif()

################################################################################

project(CCDQueryIO LANGUAGES CXX)

option(CCD_IO_BUILD_TESTS             "Build unit-tests"            ${CCD_IO_TOPLEVEL_PROJECT})
option(CCD_IO_DOWNLOAD_SAMPLE_QUERIES "Download sample CCD queries"                        OFF)

# Set default minimum C++ standard
if(CCD_IO_TOPLEVEL_PROJECT)
  set(CMAKE_CXX_STANDARD 17)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  set(CMAKE_CXX_EXTENSIONS OFF)
endif()

### Configuration
set(CCD_IO_SOURCE_DIR "${PROJECT_SOURCE_DIR}/src/ccd_io")
set(CCD_IO_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/src")

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/ccd_io/")
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/recipes/")
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/find/")

# General CMake utils
include(ccd_io_cpm_cache)
include(ccd_io_use_colors)

# Generate position-independent code by default
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

################################################################################
# IPC Toolkit Library
################################################################################

# Add an empty library and fill in the list of sources in `src/CMakeLists.txt`.
add_library(ccd_io_ccd_io)
add_library(ccd_io::ccd_io ALIAS ccd_io_ccd_io)

# Fill in configuration options
# configure_file(
#   "${CCD_IO_SOURCE_DIR}/config.hpp.in"
#   "${CCD_IO_SOURCE_DIR}/config.hpp")

# Add source and header files to ccd_io
add_subdirectory("${CCD_IO_SOURCE_DIR}")

# Public include directory for IPC Toolkit
target_include_directories(ccd_io_ccd_io PUBLIC "${CCD_IO_INCLUDE_DIR}")

################################################################################
# Optional Definitions
################################################################################

# For MSVC, do not use the min and max macros.
target_compile_definitions(ccd_io_ccd_io PUBLIC NOMINMAX)

################################################################################
# Dependencies
################################################################################

# rational-cpp (requires GMP)
include(rational_cpp)
target_link_libraries(ccd_io_ccd_io PUBLIC rational::rational)

# JSON
include(json)
target_link_libraries(ccd_io_ccd_io PUBLIC nlohmann_json::nlohmann_json)

# Logger
include(spdlog)
target_link_libraries(ccd_io_ccd_io PUBLIC spdlog::spdlog)

if(CCD_IO_DOWNLOAD_SAMPLE_QUERIES)
  include(sample_queries)
  target_link_libraries(ccd_io_ccd_io PUBLIC ccd_io::sample_queries)
endif()

# Warnings
include(ccd_io_warnings)
target_link_libraries(ccd_io_ccd_io PRIVATE ccd_io::warnings)

################################################################################
# Compiler options
################################################################################

# Use C++17
target_compile_features(ccd_io_ccd_io PUBLIC cxx_std_17)

################################################################################
# Tests
################################################################################

# Enable unit testing at the root level
if(CCD_IO_TOPLEVEL_PROJECT AND CCD_IO_BUILD_TESTS)
  include(CTest)
  enable_testing()
  add_subdirectory(tests)
endif()