#include <Eigen/Dense>

using Eigen::Array2d;
using Eigen::Vector3d;

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <queue>
#include <span>

static constexpr double MinDeltaUV = 1e-4;

template <typename T> static T lerp(double t, T const &t0, T const &t1) { return (1 - t) * t0 + t * t1; }

struct Bounds2d {
	Array2d pMin, pMax;

	Bounds2d() = default;
	Bounds2d(Array2d const &p1, Array2d const &p2) :
		pMin { std::min(p1[0], p2[0]), std::min(p1[1], p2[1]) },
		pMax { std::max(p1[0], p2[0]), std::max(p1[1], p2[1]) } {
	}

	Array2d operator[](int i) const { return i == 0 ? pMin : pMax; }

	Bounds2d operator&(Bounds2d const &o) const {
		Bounds2d ret;
		ret.pMin = pMin.max(o.pMin);
		ret.pMax = pMax.min(o.pMax);
		return ret;
	}

	bool isDegenerate() const { return pMin.x() > pMax.x() || pMin.y() > pMax.y(); }
	bool isInside(Array2d const &o) const { return (o.x() >= pMin.x() && o.x() <= pMax.x() && o.y() >= pMin.y() && o.y() <= pMax.y()); }

	Array2d diagonal() const { return pMax - pMin; }
	Array2d corner(int i) const { return Array2d((*this)[(i & 1)][0], (*this)[(i & 2) ? 1 : 0][1]); }
};

struct DividedPatch {
	Bounds2d uvB1;
	Bounds2d uvB2;
	double dLower;

	DividedPatch(Bounds2d const &uvB1, Bounds2d const &uvB2, double dLower = 0) : uvB1 { uvB1 }, uvB2 { uvB2 }, dLower { dLower } { }
	bool operator<(DividedPatch const &o) const { return dLower > o.dLower; }
};

// Bicubic Bezier Functions
template <typename P>
static P blossomCubicBezier(std::span<P const> p, double u0, double u1, double u2) {
	P a[3] = { lerp(u0, p[0], p[1]), lerp(u0, p[1], p[2]), lerp(u0, p[2], p[3]) };
	P b[2] = { lerp(u1, a[0], a[1]), lerp(u1, a[1], a[2]) };
	return lerp(u2, b[0], b[1]);
}

template <typename P>
static P blossomBicubicBezier(std::span<P const> cp, Array2d const &uv0, Array2d const &uv1, Array2d const &uv2) {
	std::array<P, 4> q;
	for (int i = 0; i < 4; i++) {
		q[i] = blossomCubicBezier(cp.subspan(i * 4, 4), uv0.y(), uv1.y(), uv2.y());
	}
	return blossomCubicBezier<P>(q, uv0.x(), uv1.x(), uv2.x());
}

static Vector3d evaluateBicubicBezier(std::span<Vector3d const> cp, Array2d const &uv, Vector3d *dpdu = nullptr, Vector3d *dpdv = nullptr) {
	if (dpdu && dpdv) {
		std::array<Vector3d, 4> cpU, cpV;
		for (int i = 0; i < 4; i++) {
			cpU[i] = blossomCubicBezier(cp.subspan(i * 4, 4), uv.y(), uv.y(), uv.y());
			std::array<Vector3d, 3> a = { lerp(uv.x(), cp[i], cp[i + 4]),
										 lerp(uv.x(), cp[i + 4], cp[i + 8]),
										 lerp(uv.x(), cp[i + 8], cp[i + 12]) };
			std::array<Vector3d, 2> b = { lerp(uv.x(), a[0], a[1]),
										 lerp(uv.x(), a[1], a[2]) };
			cpV[i] = lerp(uv.x(), b[0], b[1]);
		}
		std::array<Vector3d, 3> cpU1 = { lerp(uv.x(), cpU[0], cpU[1]),
										lerp(uv.x(), cpU[1], cpU[2]),
										lerp(uv.x(), cpU[2], cpU[3]) };
		std::array<Vector3d, 2> cpU2 = { lerp(uv.x(), cpU1[0], cpU1[1]),
										lerp(uv.x(), cpU1[1], cpU1[2]) };
		if ((cpU2[1] - cpU2[0]).squaredNorm() > 0)
			*dpdu = 3 * (cpU2[1] - cpU2[0]);
		else {
			*dpdu = cpU[3] - cpU[0];
		}
		std::array<Vector3d, 3> cpV1 = { lerp(uv.y(), cpV[0], cpV[1]),
										lerp(uv.y(), cpV[1], cpV[2]),
										lerp(uv.y(), cpV[2], cpV[3]) };
		std::array<Vector3d, 2> cpV2 = { lerp(uv.y(), cpV1[0], cpV1[1]),
										lerp(uv.y(), cpV1[1], cpV1[2]) };
		if ((cpV2[1] - cpV2[0]).squaredNorm() > 0)
			*dpdv = 3 * (cpV2[1] - cpV2[0]);
		else {
			*dpdv = cpV[3] - cpV[0];
		}
		return blossomCubicBezier<Vector3d>(cpU, uv.x(), uv.x(), uv.x());
	}
	else {
		return blossomBicubicBezier(cp, uv, uv, uv);
	}
}

// Patch Functions
template <typename P>
static std::array<P, 16> divideBezierPatch(std::span<P const> cp, Bounds2d const &uvB) {
	std::array<P, 16> divCp;
	divCp[0] = blossomBicubicBezier(cp, uvB.corner(0), uvB.corner(0), uvB.corner(0));
	divCp[1] = blossomBicubicBezier(cp, uvB.corner(0), uvB.corner(0), uvB.corner(2));
	divCp[2] = blossomBicubicBezier(cp, uvB.corner(0), uvB.corner(2), uvB.corner(2));
	divCp[3] = blossomBicubicBezier(cp, uvB.corner(2), uvB.corner(2), uvB.corner(2));
	divCp[4] = blossomBicubicBezier(cp, uvB.corner(0), uvB.corner(0), uvB.corner(1));
	divCp[5] = blossomBicubicBezier(cp, uvB.corner(0), uvB.corner(0), uvB.corner(3));
	divCp[6] = blossomBicubicBezier(cp, uvB.corner(0), uvB.corner(2), uvB.corner(3));
	divCp[7] = blossomBicubicBezier(cp, uvB.corner(2), uvB.corner(2), uvB.corner(3));
	divCp[8] = blossomBicubicBezier(cp, uvB.corner(0), uvB.corner(1), uvB.corner(1));
	divCp[9] = blossomBicubicBezier(cp, uvB.corner(0), uvB.corner(1), uvB.corner(3));
	divCp[10] = blossomBicubicBezier(cp, uvB.corner(0), uvB.corner(3), uvB.corner(3));
	divCp[11] = blossomBicubicBezier(cp, uvB.corner(2), uvB.corner(3), uvB.corner(3));
	divCp[12] = blossomBicubicBezier(cp, uvB.corner(1), uvB.corner(1), uvB.corner(1));
	divCp[13] = blossomBicubicBezier(cp, uvB.corner(1), uvB.corner(1), uvB.corner(3));
	divCp[14] = blossomBicubicBezier(cp, uvB.corner(1), uvB.corner(3), uvB.corner(3));
	divCp[15] = blossomBicubicBezier(cp, uvB.corner(3), uvB.corner(3), uvB.corner(3));
	return divCp;
}

std::array<Vector3d, 16> g_Cp1 = {
	Vector3d(0, 0, 0), Vector3d(1, 0, 0), Vector3d(2, 0, 0), Vector3d(3, 0, 0),
	Vector3d(0, 1, 0), Vector3d(1, 1, 0), Vector3d(2, 1, 0), Vector3d(3, 1, 0),
	Vector3d(0, 2, 0), Vector3d(1, 2, 0), Vector3d(2, 2, 0), Vector3d(3, 2, 0),
	Vector3d(0, 3, 0), Vector3d(1, 3, 0), Vector3d(2, 3, 0), Vector3d(3, 3, 0),
};
std::array<Vector3d, 16> g_Cp2 = {
	Vector3d(0, 0, 2), Vector3d(1, 0, 2), Vector3d(2, 0, 2), Vector3d(3, 0, 2),
	Vector3d(0, 1, 2), Vector3d(1, 1, 2), Vector3d(2, 1, 2), Vector3d(3, 1, 2),
	Vector3d(0, 2, 2), Vector3d(1, 2, 2), Vector3d(2, 2, 2), Vector3d(3, 2, 2),
	Vector3d(0, 3, 2), Vector3d(1, 3, 2), Vector3d(2, 3, 2), Vector3d(3, 3, 2),
};

std::uint64_t cnt;

static double getDist(double min1, double max1, double min2, double max2) {
	if (max1 <= min2) {
		return min2 - max1;
	} else if (max2 <= min1) {
		return min1 - max2;
	} else {
		return 0.;
	}
}

static double getLowerBound(Bounds2d const &divUvB1, Bounds2d const &divUvB2) {
	auto const pt1 = divideBezierPatch<Vector3d>(g_Cp1, divUvB1);
	auto const pt2 = divideBezierPatch<Vector3d>(g_Cp2, divUvB2);
	Vector3d min1 = pt1[0];
	Vector3d max1 = pt1[0];
	Vector3d min2 = pt2[0];
	Vector3d max2 = pt2[0];
	for (int i = 1; i < 16; i++) {
		min1 = min1.cwiseMin(pt1[i]);
		max1 = max1.cwiseMax(pt1[i]);
		min2 = min2.cwiseMin(pt2[i]);
		max2 = max2.cwiseMax(pt2[i]);
	}
	double const dx = getDist(min1.x(), max1.x(), min2.x(), max2.x()); 
	double const dy = getDist(min1.y(), max1.y(), min2.y(), max2.y()); 
	double const dz = getDist(min1.z(), max1.z(), min2.z(), max2.z()); 
	return dx * dx + dy * dy + dz * dz;
}

static double solve() {
	std::priority_queue<DividedPatch> heap;
	heap.emplace(Bounds2d(Array2d(0, 0), Array2d(1, 1)), Bounds2d(Array2d(0, 0), Array2d(1, 1)));

	while (!heap.empty()) {
		auto const cur = heap.top();
		heap.pop();
		cnt++;
		// Set uv of the middle point
		Array2d uvMid1 = (cur.uvB1.pMin + cur.uvB1.pMax) / 2;
		Array2d uvMid2 = (cur.uvB2.pMin + cur.uvB2.pMax) / 2;

		// Decide whether the algorithm converges
		if (std::max(cur.uvB1.diagonal().maxCoeff(), cur.uvB2.diagonal().maxCoeff()) < MinDeltaUV) {
			std::cout << ((cur.uvB1.pMin + cur.uvB1.pMax) / 2).transpose() << std::endl;
			std::cout << ((cur.uvB2.pMin + cur.uvB2.pMax) / 2).transpose() << std::endl;
			return cur.dLower;
		}

		// Divide the current patch into four-to-four pieces
		for (int i = 0; i < 4; i++) {
			Bounds2d divUvB1(cur.uvB1.corner(i), uvMid1);
			for (int j = 0; j < 4; j++) {
				Bounds2d divUvB2(cur.uvB2.corner(j), uvMid2);

				double const dLower = getLowerBound(divUvB1, divUvB2);
				if (dLower < std::numeric_limits<double>::infinity())
					heap.emplace(divUvB1, divUvB2, dLower);
			}
		}
	}

	return false;
}



static double solveNaive() {
	constexpr double eps = .01;
	double ans = std::numeric_limits<double>::infinity();
	Array2d ap1, ap2;
	for (double x1 = 0; x1 <= 1; x1 += eps) {
		for (double y1 = 0; y1 <= 1; y1 += eps) {
			Vector3d const p1 = evaluateBicubicBezier(g_Cp1, { x1, y1 });
			for (double x2 = 0; x2 <= 1; x2 += eps) {
				for (double y2 = 0; y2 <= 1; y2 += eps) {
					Vector3d const p2 = evaluateBicubicBezier(g_Cp2, { x2, y2 });
					double const dist2 = (p1 - p2).squaredNorm();
					if (dist2 < ans) {
						ans = dist2;
						ap1 = { x1, y1 };
						ap2 = { x2, y2 };
					}
				}
			}
		}
	}
	std::cout << ap1.transpose() << std::endl << ap2.transpose() << std::endl;
	return ans;
}

int main() {
	std::srand(std::time(nullptr));
	for (int i = 0; i < 16; i++) {
		g_Cp1[i] += Vector3d::Random() * .6 - Vector3d::Constant(.3);
		g_Cp2[i] += Vector3d::Random() * .6 - Vector3d::Constant(.3);
	}
	std::cout << solveNaive() << std::endl;
	std::cout << solve() << std::endl << cnt << std::endl;
	return 0;
}
