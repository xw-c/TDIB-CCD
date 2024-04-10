################################################################################
# Tests
################################################################################

add_executable(rational_tests test_rational.cpp)

################################################################################
# Required Libraries
################################################################################

target_link_libraries(rational_tests PUBLIC rational::rational)

include(catch2)
target_link_libraries(rational_tests PUBLIC Catch2::Catch2WithMain)

# Extra warnings (link last for highest prationalrity)
include(rational_warnings)
target_link_libraries(rational_tests PRIVATE rational::warnings)

################################################################################
# Register tests
################################################################################

FetchContent_GetProperties(catch2)
list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
include(Catch)

# Register tests
set(PARSE_CATCH_TESTS_ADD_TO_CONFIGURE_DEPENDS ON)
catch_discover_tests(rational_tests)
