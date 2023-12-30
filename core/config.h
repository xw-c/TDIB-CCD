# pragma once
#include <algorithm>
#include <array>
#include <string>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <limits>
#include <queue>
#include <span>
#include <chrono>
#include <cmath>
#include <random>

#include <Eigen/Dense>

using Eigen::Array2d;
using Eigen::Vector3d;
using Eigen::MatrixXd;

static constexpr double MinDeltaUV = 1e-6;
static constexpr double MinSquaredDist = 1e-12;
static constexpr double Epsilon = 1e-6;
static constexpr double DeltaT = 1;
bool DEBUG = 0;
static constexpr bool SHOWANS = 0;
enum class BoundingBoxType { AABB, OBB, DOP14 };
enum class PatchType { TriBezier = 10, RecBezier = 16};

std::normal_distribution<double> randNormal(0.0, 1.0); // 均值为0，标准差为1的正态分布

std::uint64_t cnt;
