set(SOURCES
  ccd_query.hpp
  logger.cpp
  logger.hpp
  read_ccd_queries.cpp
  read_ccd_queries.hpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(ccd_io_ccd_io PRIVATE ${SOURCES})