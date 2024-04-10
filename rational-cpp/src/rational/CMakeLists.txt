set(SOURCES
  rational.cpp
  rational.hpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(rational_rational PRIVATE ${SOURCES})