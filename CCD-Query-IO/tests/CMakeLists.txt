################################################################################
# Tests
################################################################################

add_executable(ccd_io_tests
  test_read_ccd_queries.cpp
)

################################################################################
# Required Libraries
################################################################################

target_link_libraries(ccd_io_tests PUBLIC ccd_io::ccd_io)

include(catch2)
target_link_libraries(ccd_io_tests PUBLIC Catch2::Catch2WithMain)

include(sample_queries)
target_link_libraries(ccd_io_tests PUBLIC ccd_io::sample_queries)

# Extra warnings (link last for highest pccd_iority)
include(ccd_io_warnings)
target_link_libraries(ccd_io_tests PRIVATE ccd_io::warnings)

################################################################################
# Compiler options
################################################################################

# Test data directory definition
if(WIN32)
  target_compile_definitions(ccd_io_tests PUBLIC CCD_IO_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}\\\\data\\\\")
else()
  target_compile_definitions(ccd_io_tests PUBLIC CCD_IO_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/data/")
endif()

################################################################################
# Register tests
################################################################################

FetchContent_GetProperties(catch2)
list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
include(Catch)

# Register tests
set(PARSE_CATCH_TESTS_ADD_TO_CONFIGURE_DEPENDS ON)
catch_discover_tests(ccd_io_tests)
