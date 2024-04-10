# Detects whether this is a top-level project
get_directory_property(HAS_PARENT PARENT_DIRECTORY)
if(HAS_PARENT)
  set(RATIONAL_TOPLEVEL_PROJECT OFF)
else()
  set(RATIONAL_TOPLEVEL_PROJECT ON)
endif()

# Check required CMake version
set(REQUIRED_CMAKE_VERSION "3.18.0")
if(RATIONAL_TOPLEVEL_PROJECT)
  cmake_minimum_required(VERSION ${REQUIRED_CMAKE_VERSION})
else()
  # Don't use cmake_minimum_required here to avoid implicitly overriding parent policies
  if(${CMAKE_VERSION} VERSION_LESS ${REQUIRED_CMAKE_VERSION})
    message(FATAL_ERROR "CMake required version to build IPC Toolkit is ${REQUIRED_CMAKE_VERSION}")
  endif()
endif()

# Include user-provided default options if available. We do that before the main
# `project()` so that we can define the C/C++ compilers from the option file.
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/RationalOptions.cmake)
  message(STATUS "Using local options file: ${CMAKE_CURRENT_SOURCE_DIR}/RationalOptions.cmake")
  include(${CMAKE_CURRENT_SOURCE_DIR}/RationalOptions.cmake)
endif()

# Enable ccache if available
find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
  option(RATIONAL_WITH_CCACHE "Enable ccache when building Rational" ${RATIONAL_TOPLEVEL_PROJECT})
else()
  option(RATIONAL_WITH_CCACHE "Enable ccache when building Rational" OFF)
endif()
if(RATIONAL_WITH_CCACHE AND CCACHE_PROGRAM)
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

project(Rational LANGUAGES CXX)

option(RATIONAL_BUILD_TESTS "Build unit-tests" ${RATIONAL_TOPLEVEL_PROJECT})

# Set default minimum C++ standard
if(RATIONAL_TOPLEVEL_PROJECT)
  set(CMAKE_CXX_STANDARD 17)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  set(CMAKE_CXX_EXTENSIONS OFF)
endif()

### Configuration
set(RATIONAL_SOURCE_DIR "${PROJECT_SOURCE_DIR}/src/rational")
set(RATIONAL_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/src")

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/rational/")
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/recipes/")
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/find/")

# General CMake utils
include(rational_cpm_cache)
include(rational_use_colors)

# Generate position-independent code by default
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

################################################################################
# IPC Toolkit Library
################################################################################

# Add an empty library and fill in the list of sources in `src/CMakeLists.txt`.
add_library(rational_rational)
add_library(rational::rational ALIAS rational_rational)

# Add source and header files to rational
add_subdirectory("${RATIONAL_SOURCE_DIR}")

# Public include directory for IPC Toolkit
target_include_directories(rational_rational PUBLIC "${RATIONAL_INCLUDE_DIR}")

################################################################################
# Optional Definitions
################################################################################

# For MSVC, do not use the min and max macros.
target_compile_definitions(rational_rational PUBLIC NOMINMAX)

################################################################################
# Dependencies
################################################################################

# GMP
find_package(GMP REQUIRED)
target_link_libraries(rational_rational PUBLIC gmp::gmp)

# Warnings
include(rational_warnings)
target_link_libraries(rational_rational PRIVATE rational::warnings)

################################################################################
# Compiler options
################################################################################

# Use C++17
target_compile_features(rational_rational PUBLIC cxx_std_17)

################################################################################
# Tests
################################################################################

# Enable unit testing at the root level
if(RATIONAL_TOPLEVEL_PROJECT AND RATIONAL_BUILD_TESTS)
  include(CTest)
  enable_testing()
  add_subdirectory(tests)
endif()