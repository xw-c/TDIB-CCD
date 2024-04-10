#pragma once

#include <array>

namespace ccd_io {
struct CCDQuery {
    /// @brief The vertices of the query.
    /// Order:
    ///   Point-Triangle:
    ///     p_t0, t0_t0, t1_t0, t2_t0, p_t1, t0_t1, t1_t1, t2_t1
    ///   Edge-Edge:
    ///     ea0_t0, ea1_t0, eb0_t0, eb1_t0, ea0_t1, ea1_t1, eb0_t1, eb1_t1
    std::array<std::array<double, 3>, 8> vertices;

    /// @brief The ground truth result of the query.
    bool ground_truth;
};
} // namespace ccd_io