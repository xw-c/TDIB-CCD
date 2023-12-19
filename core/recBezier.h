#include <Eigen/Dense>

using Eigen::Array2d;
using Eigen::Vector3d;
using Eigen::MatrixXd;

#include<iostream>

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

struct DividedPatches {
	Bounds2d uvB1;
	Bounds2d uvB2;
	double tLower;

	DividedPatches(Bounds2d const &uvB1, Bounds2d const &uvB2, double tLower = std::numeric_limits<double>::infinity()) : uvB1 { uvB1 }, uvB2 { uvB2 }, tLower { tLower } { }
	bool operator<(DividedPatches const &o) const { return tLower > o.tLower; }
};

struct DividedPatch {
	Bounds2d uvB;
	double tLower;

	DividedPatch(Bounds2d const &uvB, double tLower = std::numeric_limits<double>::infinity()) : uvB { uvB }, tLower { tLower } { }
	bool operator<(DividedPatch const &o) const { return tLower > o.tLower; }
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
