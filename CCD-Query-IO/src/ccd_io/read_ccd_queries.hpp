#pragma once

#include <ccd_io/ccd_query.hpp>

#include <vector>
#include <string>

namespace ccd_io {

/// @brief Read CCD Queries from a file in rational CSV format.
/// @param filename The name of the file to read.
/// @return A vector of CCD queries.
std::vector<CCDQuery> read_ccd_queries(const std::string& filename);

/// @brief Read CCD Queries from a file in rational CSV format.
/// @param vertices_filename The name of the file containing the vertices.
/// @param ground_truth_filename The name of the file containing the ground
/// @return A vector of CCD queries.
std::vector<CCDQuery> read_ccd_queries(
    const std::string& vertices_filename,
    const std::string& ground_truth_filename);

} // namespace ccd_io
