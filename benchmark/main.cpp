#include <ccd_io/read_ccd_queries.hpp>
#include <chrono>
#include <fstream>
#include <cassert>
#include <filesystem>
#include <iostream>
#include "linearSolver.h"
#include "linearSolverTrad.h"

static std::ofstream outfile("log.txt");

void count_fp_fn(double t, bool gt, int& fn, int& fp, const std::filesystem::directory_entry& csv, int i){
    bool collide = t >= 0;
    if(gt && !collide) {
        fn++;
        // std::cout<<"FN occured in "<<csv.path().string()<<" case #"<<i<<std::endl;
        outfile<<"FN occured in "<<csv.path().string()<<" case #"<<i<<std::endl;
    }
    if(!gt && collide) {
        fp++;
        // std::cout<<"FP occured in "<<csv.path().string()<<" case #"<<i<<std::endl;
        outfile<<"FP occured in "<<csv.path().string()<<" case #"<<i<<std::endl;
    }
}

double do_EE_test(const std::array<std::array<double, 3>, 8>& vertices, double& time){
    Vector3d ea0_t0(vertices[0][0],vertices[0][1],vertices[0][2]);
    Vector3d ea1_t0(vertices[1][0],vertices[1][1],vertices[1][2]);
    Vector3d eb0_t0(vertices[2][0],vertices[2][1],vertices[2][2]);
    Vector3d eb1_t0(vertices[3][0],vertices[3][1],vertices[3][2]);
    Vector3d ea0_t1(vertices[4][0],vertices[4][1],vertices[4][2]);
    Vector3d ea1_t1(vertices[5][0],vertices[5][1],vertices[5][2]);
    Vector3d eb0_t1(vertices[6][0],vertices[6][1],vertices[6][2]);
    Vector3d eb1_t1(vertices[7][0],vertices[7][1],vertices[7][2]);
    Edge a(ea0_t0, ea1_t0), b(eb0_t0, eb1_t0), 
            va(ea0_t1-ea0_t0, ea1_t1-ea1_t0), vb(eb0_t1-eb0_t0, eb1_t1-eb1_t0);
    double u1 = 0, u2 = 0;
    const auto initialTime = std::chrono::steady_clock::now();
    double t = LinearSolverTD::solveEETest(a, va, b, vb, u1, u2, BoundingBoxType::OBB, 1, 1e-6);
    // double t = LinearSolverTrad::solveEETest(a, va, b, vb, u1, u2, BoundingBoxType::AABB, 1, 1e-6);
    const auto endTime = std::chrono::steady_clock::now();
    double used_time = std::chrono::duration<double>(endTime - initialTime).count();
    time += used_time;
    return t;
}

double do_VF_test(const std::array<std::array<double, 3>, 8>& vertices, double& time){
    Vector3d p_t0(vertices[0][0],vertices[0][1],vertices[0][2]);
    Vector3d t0_t0(vertices[1][0],vertices[1][1],vertices[1][2]);
    Vector3d t1_t0(vertices[2][0],vertices[2][1],vertices[2][2]);
    Vector3d t2_t0(vertices[3][0],vertices[3][1],vertices[3][2]);
    Vector3d p_t1(vertices[4][0],vertices[4][1],vertices[4][2]);
    Vector3d t0_t1(vertices[5][0],vertices[5][1],vertices[5][2]);
    Vector3d t1_t1(vertices[6][0],vertices[6][1],vertices[6][2]);
    Vector3d t2_t1(vertices[7][0],vertices[7][1],vertices[7][2]);
    std::array<Vector3d, 3> f{t0_t0, t1_t0, t2_t0}, vf{t0_t1 - t0_t0, t1_t1 - t1_t0, t2_t1 - t2_t0};
    Face face(f), vface(vf);
    Array2d uv(0, 0);
    const auto initialTime = std::chrono::steady_clock::now();
    double t = LinearSolverTD::solveVFTest(p_t0, p_t1-p_t0, face, vface, uv, BoundingBoxType::OBB, 1, 1e-6);
    // double t = LinearSolverTrad::solveVFTest(p_t0, p_t1-p_t0, face, vface, uv, BoundingBoxType::AABB, 1, 1e-6);
    const auto endTime = std::chrono::steady_clock::now();
    double used_time = std::chrono::duration<double>(endTime - initialTime).count();
    time += used_time;
    return t;
}

static const std::filesystem::path root_path(std::string("D:\\CCD\\wang21benchmark"));

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
    "unit-tests"
} };

static const std::array<std::string, 4> simulation_folders { {
    "chain",
    "cow-heads",
    "golf-ball",
    "mat-twist"
} };

static const std::array<std::string, 2> subfolders { {
    "edge-edge",
    "vertex-face"
} };

void handcrafted_tests(){
    double time_ee = 0, time_vf = 0;
    int fp_ee = 0, fn_ee = 0, fp_vf = 0, fn_vf = 0, cnt_ee = 0, cnt_vf = 0;
    for (const std::string& folder : folders) {
        for (const std::string& subfolder : subfolders) {
            auto dir = root_path / "ccd-queries-handcrafted" / folder / subfolder;
            for (const auto& csv : std::filesystem::directory_iterator(dir)) {
                std::vector<ccd_io::CCDQuery> queries = ccd_io::read_ccd_queries(csv.path().string());
                int len = queries.size();
                for(int i = 0; i < len; ++i) {
                    // std::cout<<"testing "<<csv.path().string()<<" case #"<< i <<std::endl;
                    const auto& query = queries[i];
                    const auto& vertices = query.vertices;
                    bool gt = query.ground_truth;
                    if(subfolder == "edge-edge") {
                        cnt_ee++;
                        double t = do_EE_test(vertices, time_ee);
                        count_fp_fn(t, gt, fn_ee, fp_ee, csv, i);
                    }
                    else if(subfolder == "vertex-face") {
                        cnt_vf++;
                        double t = do_VF_test(vertices, time_vf);
                        count_fp_fn(t, gt, fn_vf, fp_vf, csv, i);
                    }
                }
            }
        }
    }
    std::cout << std::fixed << std::setprecision(10);
    std::cout<<"FN EE:"<<fn_ee<<std::endl;
    std::cout<<"FP EE:"<<fp_ee<<std::endl;
    std::cout<<"FN VF:"<<fn_vf<<std::endl;
    std::cout<<"FP VF:"<<fp_vf<<std::endl;
    std::cout<<"total cases EE:"<<cnt_ee<<std::endl;
    std::cout<<"total cases VF:"<<cnt_vf<<std::endl;
    std::cout<<"total time EE:"<<time_ee<<std::endl;
    std::cout<<"total time VF:"<<time_vf<<std::endl;
}

void simulation_tests(){
    double time_ee = 0, time_vf = 0;
    int fp_ee = 0, fn_ee = 0, fp_vf = 0, fn_vf = 0, cnt_ee = 0, cnt_vf = 0;
    for (const std::string& simulation_folder : simulation_folders) {
        for (const std::string& subfolder : subfolders) {
            auto dir = root_path / simulation_folder / subfolder;
            for (const auto& csv : std::filesystem::directory_iterator(dir)) {
                std::vector<ccd_io::CCDQuery> queries = ccd_io::read_ccd_queries(csv.path().string());
                int len = queries.size();
                // std::cout<<"testing "<< csv.path().string() << std::endl;
                for(int i = 0; i < len; ++i) {
                    const auto& query = queries[i];
                    const auto& vertices = query.vertices;
                    bool gt = query.ground_truth;
                    if(subfolder == "edge-edge") {
                        cnt_ee++;
                        // outfile<<"testing "<< csv.path().string() <<"case #"<<i<<" cur time "<< time_ee << std::endl;
                        double t = do_EE_test(vertices, time_ee);
                        count_fp_fn(t, gt, fn_ee, fp_ee, csv, i);
                    }
                    else if(subfolder == "vertex-face") {
                        cnt_vf++;
                        double t = do_VF_test(vertices, time_vf);
                        count_fp_fn(t, gt, fn_vf, fp_vf, csv, i);
                    }
                }
            }
        }
    }
    std::cout << std::fixed << std::setprecision(10);
    std::cout<<"FN EE:"<<fn_ee<<std::endl;
    std::cout<<"FP EE:"<<fp_ee<<std::endl;
    std::cout<<"FN VF:"<<fn_vf<<std::endl;
    std::cout<<"FP VF:"<<fp_vf<<std::endl;
    std::cout<<"total cases EE:"<<cnt_ee<<std::endl;
    std::cout<<"total cases VF:"<<cnt_vf<<std::endl;
    std::cout<<"total time EE:"<<time_ee<<std::endl;
    std::cout<<"total time VF:"<<time_vf<<std::endl;
}

int main(int argc, char *argv[]){
    handcrafted_tests();
    simulation_tests();
}

// int main(int argc, char *argv[]){
//     double time = 0;
//     int fp = 0, fn = 0, cnt = 0;
//     std::string csv = "D:\\CCD\\ccd-queries-handcrafted\\fn_obb.csv";
//     std::vector<ccd_io::CCDQuery> queries = ccd_io::read_ccd_queries(csv);
//     int len = queries.size();
//     // std::cout<<len<<std::endl;
//     for(int i = 0; i < len; ++i) {
//         cnt++;
//         // std::cout<<"testing "<<csv.path().string()<<" case #"<< i <<std::endl;
//         const auto& query = queries[i];
//         const auto& vertices = query.vertices;
//         // std::cout << std::fixed << std::setprecision(18);
//         // for(int k = 0; k < 8; ++k){
//         //     for(int l = 0; l < 3; ++l){
//         //         std::cout<<vertices[k][l]<<" ";
//         //     }
//         //     std::cout<<std::endl;
//         // }
//         bool gt = query.ground_truth;
//         double t = do_EE_test(vertices, time);
//         bool collide = t >= 0;
//         if(gt && !collide) {
//             fn++;
//             std::cout<<"FN occured in "<<csv<<" case #"<<i<<std::endl;
//         }
//     }
//     std::cout << std::fixed << std::setprecision(10);
//     std::cout<<"FN :"<<fn<<std::endl;
//     std::cout<<"FP :"<<fp<<std::endl;
//     std::cout<<"total cases:"<<cnt<<std::endl;
//     std::cout<<"total time:"<<time<<std::endl;
// }