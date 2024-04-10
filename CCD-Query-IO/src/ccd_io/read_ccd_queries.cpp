#include "read_ccd_queries.hpp"

#include <ccd_io/logger.hpp>

#include <rational/rational.hpp>
#include <nlohmann/json.hpp>

#include <array>
#include <fstream>
#include <sstream>

namespace ccd_io {

std::vector<CCDQuery> read_ccd_queries(const std::string& filename)
{
    // NOTE: If the file contains N lines, N % 8 == 0 because every 8 lines are
    // a single query.
    std::vector<CCDQuery> queries;

    // std::vector<std::array<double, 3>> vs;
    // vs.clear();
    std::ifstream file(filename);
    if (!file.is_open()) {
        log_and_throw_error("Unable to open file {}!", filename);
    }

    std::string line;
    bool gt_set = false;
    int i;
    for (i = 0; file && std::getline(file, line); i++) {
        if (line[0] == '#')
            continue;

        if (i % 8 == 0) {
            // New query
            queries.emplace_back();
            gt_set = false;
        }

        std::istringstream line_stream(line);
        // the first six are one vetex, the seventh is the result
        std::array<std::string, 7> line_items;
        for (int j = 0; j < 6; j++) {
            if (!std::getline(line_stream, line_items[j], ',')) {
                log_and_throw_error(
                    "Line {} in {} contains only {} numbers!", i, filename, j);
            }
        }
        std::getline(line_stream, line_items[6], ','); // optional ground truth

        for (int j = 0; j < 3; j++) {
            queries.back().vertices[i % 8][j] = rational::Rational(
                line_items[2 * j + 0], line_items[2 * j + 1]);
            if (!std::isfinite(queries.back().vertices[i % 8][j])) {
                log_and_throw_error(
                    "Line {} in {} contains a non-integer!", i, filename);
            }
        }

        if (!line_items[6].empty()) {
            if (!gt_set) {
                queries.back().ground_truth = std::stoi(line_items[6]);
                gt_set = true;
            } else if (
                bool(std::stoi(line_items[6])) != queries.back().ground_truth) {
                log_and_throw_error(
                    "Ground truth mismatch in file {} line {}!", filename, i);
            }
        }
    }

    if (i % 8 != 0) {
        log_and_throw_error(
            "File {} has {} lines, which is not a multiple of 8!", filename, i);
    }

    if (!file.eof()) {
        log_and_throw_error("Could not read file {}!", filename);
    }

    return queries;
}

std::vector<CCDQuery> read_ccd_queries(
    const std::string& vertices_filename,
    const std::string& ground_truth_filename)
{
    std::vector<CCDQuery> queries = read_ccd_queries(vertices_filename);

    nlohmann::json ground_truths;
    {
        std::ifstream file(ground_truth_filename);
        if (!file.is_open()) {
            log_and_throw_error(
                "Unable to open file {}!", ground_truth_filename);
        }
        file >> ground_truths;
    }

    if (!ground_truths.is_array()) {
        log_and_throw_error(
            "Ground truth file {} does not contain a JSON array!",
            ground_truth_filename);
    }

    if (ground_truths.size() != queries.size()) {
        log_and_throw_error(
            "Ground truth file {} has {} entries, but {} has {} entries!",
            ground_truth_filename, ground_truths.size(), vertices_filename,
            queries.size());
    }

    for (size_t i = 0; i < queries.size(); i++) {
        if (!ground_truths[i].is_boolean()) {
            log_and_throw_error(
                "Ground truth file {} entry {} is not a boolean!",
                ground_truth_filename, i);
        }

        queries[i].ground_truth = ground_truths[i];
    }

    return queries;
}

} // namespace ccd_io
