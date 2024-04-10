#include <catch2/catch_test_macros.hpp>

#include <ccd_io/read_ccd_queries.hpp>
#include <ccd_io/logger.hpp>

#include <fstream>
#include <filesystem>

using namespace ccd_io;

void check_fails_to_read(std::string filename)
{
    try {
        read_ccd_queries(filename);
        FAIL();
    } catch (...) {
        SUCCEED();
    }
}

TEST_CASE("Nonexistent input", "[read]")
{
    logger().set_level(spdlog::level::off);
    check_fails_to_read("nonexistent_file.csv");
}

TEST_CASE("Bad input", "[read]")
{
    logger().set_level(spdlog::level::off);

    std::function<std::string(int)> line;
    SECTION("Too few lines")
    {
        line = [](int i) { return "0,1,0,1,0,1,0"; };
    }
    SECTION("Too few items per line")
    {
        line = [](int i) { return "0,1,0,1,0"; };
    }
    SECTION("Ground truth mismatch")
    {
        line = [](int i) {
            return fmt::format("0,1,0,1,0,1,{}", i % 8 == 0 ? 1 : 0);
        };
    }
    SECTION("Not a number")
    {
        line = [](int i) { return "a,1,0,1,0,1,1"; };
    }

    {
        std::ofstream of("bad_input.csv");
        for (int i = 0; i < 7; i++) {
            of << line(i) << "\n";
        }
    }

    check_fails_to_read("bad_input.csv");
}

TEST_CASE("Read CCD Query", "[read]")
{
    constexpr int N = 2;
    constexpr double r = 0.25;
    {
        std::ofstream of("test_read_ccd_queries.csv");
        for (int i = 0; i < 8 * N; i++) {
            of << "1,4,1,4,1,4,1\n";
        }
    }

    std::vector<CCDQuery> queries =
        read_ccd_queries("test_read_ccd_queries.csv");

    REQUIRE(queries.size() == N);
    for (int q = 0; q < N; q++) {
        CHECK(queries[q].vertices.size() == 8);
        CHECK(queries[q].ground_truth == 1);
        for (int i = 0; i < 8; i++) {
            for (int j = 0; j < 3; j++) {
                CHECK(queries[q].vertices[i][j] == r);
            }
        }
    }
}

TEST_CASE("Read Sample Queries", "[read][sample]")
{
    namespace fs = std::filesystem;

    static const fs::path root_path(std::string(CCD_IO_SAMPLE_QUERIES_DIR));

    static const std::array<std::string, 15> folders { {
        "chain",
        "cow-heads",
        "erleben-sliding-spike",
        "erleben-spike-wedge",
        "erleben-sliding-wedge",
        "erleben-wedge-crack",
        "erleben-spike-crack",
        "erleben-wedges",
        "erleben-cube-cliff-edges",
        "erleben-spike-hole",
        "erleben-cube-internal-edges",
        "erleben-spikes",
        "golf-ball",
        "mat-twist",
        "unit-tests",
    } };

    static const std::array<std::string, 2> subfolders { {
        "edge-edge",
        "vertex-face",
    } };

    for (const std::string& folder : folders) {
        for (const std::string& subfolder : subfolders) {
            auto dir = root_path / folder / subfolder;
            for (const auto& csv : fs::directory_iterator(dir)) {
                try {
                    read_ccd_queries(csv.path().string());
                    SUCCEED();
                } catch (std::runtime_error& e) {
                    logger().error(
                        "Failed to read {}: {}", csv.path().string(), e.what());
                    FAIL();
                }
            }
        }
    }
}

TEST_CASE("Read Split Sample Queries", "[read][sample]")
{
    std::vector<CCDQuery> queries = read_ccd_queries(
        std::string(CCD_IO_TEST_DATA_DIR) + "/split-ground-truth/queries.csv",
        std::string(CCD_IO_TEST_DATA_DIR)
            + "/split-ground-truth/mma_bool.json");
    REQUIRE(queries.size() == 2);
    CHECK(queries[0].ground_truth == false);
    CHECK(queries[1].ground_truth == false);
}