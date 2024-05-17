#include <ccd_io/read_ccd_queries.hpp>
#include <fstream>
#include <cassert>
#include <filesystem>
#include <iostream>
#include "linearSolver.h"
int main(){
    static const std::filesystem::path root_path(std::string("D:\\CCD\\ccd-queries-handcrafted"));

    static const std::array<std::string, 11> folders { {
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
        "unit-tests",
    } };

    static const std::array<std::string, 2> subfolders { {
        "edge-edge",
        "vertex-face",
    } };

    for (const std::string& folder : folders) {
        for (const std::string& subfolder : subfolders) {
            auto dir = root_path / folder / subfolder;
            for (const auto& csv : std::filesystem::directory_iterator(dir)) {
                std::vector<ccd_io::CCDQuery> queries = ccd_io::read_ccd_queries(csv.path().string());
                std::cout<<queries[0].vertices[0][0]<<std::endl;
            }
        }
    }
}