# pragma once
#include <algorithm>
#include <array>
#include <string>
#include <cstdint>
#include <cfloat>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <limits>
#include <queue>
#include <set>
#include <span>
#include <chrono>
#include <cmath>
#include <random>

#include <Eigen/Dense>

using Eigen::Array2d;
using Eigen::Vector3d;
using Eigen::Vector4d;
using Eigen::MatrixXd;
using Array2dError = std::pair<Array2d, Array2d>;
const double INFT = std::numeric_limits<double>::infinity();
const double PI = std::acos(-1);

static constexpr int KaseDefault = 100;
static constexpr double MinL1Dist = 1e-6; // 
static constexpr double DeltaT = 1;

// for multi-point detection
static constexpr double MeantimeEpsilon = 1e-3;
static constexpr double SeparationEucDist = 1e-0;
// static constexpr double SeparationUVDist = 1e-0;

std::normal_distribution<double> randNormal(0.0, 1.0);
std::default_random_engine randGenerator(0);

enum class BoundingBoxType { AABB, OBB };
enum class SolverType { TDIntv, TradIntv };
// BoundingBoxType BBDefault = BoundingBoxType::OBB;
// SolverType SolverDefault = SolverType::BaseIntv;
