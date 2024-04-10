# CCD Query IO

Read CCD queries in rational CSV data format using GMP.

### Existing Dataset

This reader is designed to read the CCD Dataset available [here](https://archive.nyu.edu/handle/2451/61518). You can also find a sample set of these queries on [GitHub](https://github.com/Continuous-Collision-Detection/Sample-Queries).

This dataset was generated as part of our paper [A Large Scale Benchmark and an Inclusion-Based Algorithm for Continuous Collision Detection [Wang et al. 2021]](https://continuous-collision-detection.github.io/tight_inclusion/).

## Build

The easiest way to add the reader to an existing CMake project is to download it through CMake.
CMake provides functionality for doing this called [FetchContent](https://cmake.org/cmake/help/latest/module/FetchContent.html) (requires CMake ≥ 3.14).
We use this same process to download all external dependencies.

For example,

```cmake
include(FetchContent)
FetchContent_Declare(
    ccd_io
    GIT_REPOSITORY https://github.com/Continuous-Collision-Detection/CCD-Query-IO.git
    GIT_TAG ${CCD_IO_GIT_TAG}
)
FetchContent_MakeAvailable(ccd_io)
```

where `CCD_IO_GIT_TAG` is set to the version of the toolkit you want to use. This will download and add the reader to CMake. The reader can then be linked against using

```cmake
target_link_libraries(${PROJECT_NAME} PUBLIC ccd_io::ccd_io)
```

where `PROJECT_NAME` is the name of your library/binary.

### Dependencies

Almost all required dependencies are downloaded through CMake depending on the build options.

The following libraries are used in this project:

* [GMP](https://gmplib.org/): reading rational numbers from strings
    * GMP must be installed at a system level
* [Nlohmann JSON](https://github.com/nlohmann/json): reading JSON files
* [spdlog](https://github.com/gabime/spdlog): logging information

## Usage

```c++
#include <ccd_io/read_ccd_queries.hpp>

// If the ground truth is in the same CSV file
std::vector<ccd_io::CCDQuery> queries = ccd_io::read_ccd_queries("queries.csv");
// or if the ground truth is in a separate JSON file
std::vector<ccd_io::CCDQuery> queries = ccd_io::read_ccd_queries("queries.csv", "ground_truth.json");
```

CCDQuery is defined as follows:

```c++
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
```

## File Format

The queries are stored in a CSV file with the following format:
* `8*n` rows where `n` is the number of queries
* every 8 rows cooresponds to a single query
* each row cooresponds to a vertex's position (and optionally the ground truth)
* the order of the rows is:
    * Edge-edge queries: `ea0_t0, ea1_t0, eb0_t0, eb1_t0, ea0_t1, ea1_t1, eb0_t1, eb1_t1` where `a` or `b` indicates which edge, `0` or `1` indicates the end point of the edge, and `t0` or `t1` indicates the starting or ending positions of the edge.
    * Vertex-face queries: `v_t0, f0_t0, f1_t0, f2_t0, v_t1, f0_t1, f1_t1, f2_t1` where `v` indicates the vertex, `f0`, `f1`, or `f2` indicates the face vertices, and `t0` or `t1` indicates the starting or ending positions of the vertex.
* `6` columns of integers (and optionally a seventh column for the boolean ground truth)
    * If the ground truth is included, the last column is `ground_truth` as `0` or `1`. Note, the ground truth for a single query should be identical for all 8 rows.
* each pair of columns (`1`/`2`, `3`/`4`, and `5`/`6`) comprise a single rational number for the `x`, `y`, and `z` coordinates of the vertex, respectively.

### Optional Separate Ground Truth File

Optionally, the ground truth can be stored in a separate JSON file. The JSON file should be an array of boolean values where each value cooresponds to the ground truth of a single query. The order of the values should match the order of the queries in the CSV file.

## Contributing

Contributions are welcome! Please submit a pull request with your changes.

## License

This project is licensed under the MIT License - see the [LICENSE](https://github.com/Continuous-Collision-Detection/CCD-Query-IO/blob/main/LICENSE) file for details.

