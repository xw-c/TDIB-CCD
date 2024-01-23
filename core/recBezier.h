# pragma once
#include <Eigen/Dense>

using Eigen::Array2d;
using Eigen::Vector3d;
using Eigen::MatrixXd;

#include <iostream>
#include <span>

class RecParamBound {
public:
	Array2d pMin, pMax;

	// RecParamBound() = default;
	RecParamBound(Array2d const &p1 = Array2d(0,0), Array2d const &p2 = Array2d(1,1)) :
		pMin { std::min(p1[0], p2[0]), std::min(p1[1], p2[1]) },
		pMax { std::max(p1[0], p2[0]), std::max(p1[1], p2[1]) } {
	}

	Array2d operator[](int i) const { return i == 0 ? pMin : pMax; }

	RecParamBound operator&(RecParamBound const &o) const {
		RecParamBound ret;
		ret.pMin = pMin.max(o.pMin);
		ret.pMax = pMax.min(o.pMax);
		return ret;
	}

	bool isDegenerate() const { return pMin.x() > pMax.x() || pMin.y() > pMax.y(); }
	bool isInside(Array2d const &o) const { return (o.x() >= pMin.x() && o.x() <= pMax.x() && o.y() >= pMin.y() && o.y() <= pMax.y()); }

	Array2d diagonal() const { return pMax - pMin; }
	Array2d corner(int i) const { return Array2d((*this)[(i & 1)][0], (*this)[(i & 2) ? 1 : 0][1]); }

	RecParamBound interpSubpatchParam(const int id) const {
		Array2d pMid = 0.5 * (pMin + pMax);
		return RecParamBound(corner(id), pMid);
	}

	Array2d centerParam() const { return 0.5 * (pMin + pMax); }
	double diameter() const {
		return (pMax - pMin).maxCoeff();
	}
};

class RecLinearBezier{
public:
	static const int cntCp = 4;
	std::array<Vector3d, 4> ctrlp;
	Vector3d lerp(double t, Vector3d const &t0, Vector3d const &t1) const { return (1 - t) * t0 + t * t1; }
	Vector3d blossomLinearBezier(std::span<Vector3d const> p, double u0) const {
		return lerp(u0, p[0], p[1]);
	}

	Vector3d blossomBilinearBezier(std::span<Vector3d const> cp, Array2d const &uv0) const {
		std::array<Vector3d, 2> q;
		for (int i = 0; i < 2; i++) {
			q[i] = blossomLinearBezier(cp.subspan(i * 2, 2), uv0.y());
		}
		return blossomLinearBezier(q, uv0.x());
	}

	Vector3d evaluatePatchPoint(Array2d const &uv) const {
		return blossomBilinearBezier(ctrlp, uv);
	}

	static double feasibleUpperV(const double& u) { return 1; }
	static Vector3d axisU(const std::array<Vector3d, 4>& pt){
		return pt[2]-pt[0]+pt[3]-pt[1];
	}
	static Vector3d axisV(const std::array<Vector3d, 4>& pt){
		return pt[1]-pt[0]+pt[3]-pt[2];
	}
	
	// Patch Functions
	// 先v变再u变
	std::array<Vector3d, 4> divideBezierPatch(RecParamBound const &uvB) const {
		std::array<Vector3d, 4> divCp;
		divCp[0] = blossomBilinearBezier(ctrlp, uvB.corner(0));
		divCp[1] = blossomBilinearBezier(ctrlp, uvB.corner(2));
		divCp[2] = blossomBilinearBezier(ctrlp, uvB.corner(1));
		divCp[3] = blossomBilinearBezier(ctrlp, uvB.corner(3));
		return divCp;
	}
};

class RecQuadBezier{
public:
	static const int cntCp = 9;
	std::array<Vector3d, 9> ctrlp;
	Vector3d lerp(double t, Vector3d const &t0, Vector3d const &t1) const { return (1 - t) * t0 + t * t1; }
	Vector3d blossomQuadBezier(std::span<Vector3d const> p, double u0, double u1) const {
		Vector3d b[2] = { lerp(u0, p[0], p[1]), lerp(u0, p[1], p[2]) };
		return lerp(u1, b[0], b[1]);
	}

	Vector3d blossomBiquadBezier(std::span<Vector3d const> cp, Array2d const &uv0, Array2d const &uv1) const {
		std::array<Vector3d, 3> q;
		for (int i = 0; i < 3; i++) {
			q[i] = blossomQuadBezier(cp.subspan(i * 3, 3), uv0.y(), uv1.y());
		}
		return blossomQuadBezier(q, uv0.x(), uv1.x());
	}

	Vector3d evaluatePatchPoint(Array2d const &uv) const {
		return blossomBiquadBezier(ctrlp, uv, uv);
	}

	static double feasibleUpperV(const double& u) { return 1; }
	static Vector3d axisU(const std::array<Vector3d, 9>& pt){
		return pt[6]-pt[0]+pt[8]-pt[2];
	}
	static Vector3d axisV(const std::array<Vector3d, 9>& pt){
		return pt[2]-pt[0]+pt[8]-pt[6];
	}

	// Patch Functions
	// 先v变再u变
	std::array<Vector3d, 9> divideBezierPatch(RecParamBound const &uvB) const {
		std::array<Vector3d, 9> divCp;
		divCp[0] = blossomBiquadBezier(ctrlp, uvB.corner(0), uvB.corner(0));
		divCp[1] = blossomBiquadBezier(ctrlp, uvB.corner(0), uvB.corner(2));
		divCp[2] = blossomBiquadBezier(ctrlp, uvB.corner(2), uvB.corner(2));
		divCp[3] = blossomBiquadBezier(ctrlp, uvB.corner(0), uvB.corner(1));
		divCp[4] = blossomBiquadBezier(ctrlp, uvB.corner(0), uvB.corner(3));
		divCp[5] = blossomBiquadBezier(ctrlp, uvB.corner(2), uvB.corner(3));
		divCp[6] = blossomBiquadBezier(ctrlp, uvB.corner(1), uvB.corner(1));
		divCp[7] = blossomBiquadBezier(ctrlp, uvB.corner(1), uvB.corner(3));
		divCp[8] = blossomBiquadBezier(ctrlp, uvB.corner(3), uvB.corner(3));
		return divCp;
	}
};

// Bicubic Bezier Functions
class RecCubicBezier{
public:
	static const int cntCp = 16;
	std::array<Vector3d, 16> ctrlp;
	Vector3d lerp(double t, Vector3d const &t0, Vector3d const &t1) const { return (1 - t) * t0 + t * t1; }
	Vector3d blossomCubicBezier(std::span<Vector3d const> p, double u0, double u1, double u2) const {
		Vector3d a[3] = { lerp(u0, p[0], p[1]), lerp(u0, p[1], p[2]), lerp(u0, p[2], p[3]) };
		Vector3d b[2] = { lerp(u1, a[0], a[1]), lerp(u1, a[1], a[2]) };
		return lerp(u2, b[0], b[1]);
	}

	Vector3d blossomBicubicBezier(std::span<Vector3d const> cp, Array2d const &uv0, Array2d const &uv1, Array2d const &uv2) const {
		std::array<Vector3d, 4> q;
		for (int i = 0; i < 4; i++) {
			q[i] = blossomCubicBezier(cp.subspan(i * 4, 4), uv0.y(), uv1.y(), uv2.y());
		}
		return blossomCubicBezier(q, uv0.x(), uv1.x(), uv2.x());
	}

	Vector3d evaluatePatchPoint(Array2d const &uv) const {
		return blossomBicubicBezier(ctrlp, uv, uv, uv);
	}

	static double feasibleUpperV(const double& u) { return 1; }
	static Vector3d axisU(const std::array<Vector3d, 16>& pt){
		return pt[12]-pt[0]+pt[15]-pt[3];
	}
	static Vector3d axisV(const std::array<Vector3d, 16>& pt){
		return pt[3]-pt[0]+pt[15]-pt[12];
	}

	// Patch Functions
	// 先v变再u变
	std::array<Vector3d, 16> divideBezierPatch(RecParamBound const &uvB) const {
		std::array<Vector3d, 16> divCp;
		divCp[0] = blossomBicubicBezier(ctrlp, uvB.corner(0), uvB.corner(0), uvB.corner(0));
		divCp[1] = blossomBicubicBezier(ctrlp, uvB.corner(0), uvB.corner(0), uvB.corner(2));
		divCp[2] = blossomBicubicBezier(ctrlp, uvB.corner(0), uvB.corner(2), uvB.corner(2));
		divCp[3] = blossomBicubicBezier(ctrlp, uvB.corner(2), uvB.corner(2), uvB.corner(2));
		divCp[4] = blossomBicubicBezier(ctrlp, uvB.corner(0), uvB.corner(0), uvB.corner(1));
		divCp[5] = blossomBicubicBezier(ctrlp, uvB.corner(0), uvB.corner(0), uvB.corner(3));
		divCp[6] = blossomBicubicBezier(ctrlp, uvB.corner(0), uvB.corner(2), uvB.corner(3));
		divCp[7] = blossomBicubicBezier(ctrlp, uvB.corner(2), uvB.corner(2), uvB.corner(3));
		divCp[8] = blossomBicubicBezier(ctrlp, uvB.corner(0), uvB.corner(1), uvB.corner(1));
		divCp[9] = blossomBicubicBezier(ctrlp, uvB.corner(0), uvB.corner(1), uvB.corner(3));
		divCp[10] = blossomBicubicBezier(ctrlp, uvB.corner(0), uvB.corner(3), uvB.corner(3));
		divCp[11] = blossomBicubicBezier(ctrlp, uvB.corner(2), uvB.corner(3), uvB.corner(3));
		divCp[12] = blossomBicubicBezier(ctrlp, uvB.corner(1), uvB.corner(1), uvB.corner(1));
		divCp[13] = blossomBicubicBezier(ctrlp, uvB.corner(1), uvB.corner(1), uvB.corner(3));
		divCp[14] = blossomBicubicBezier(ctrlp, uvB.corner(1), uvB.corner(3), uvB.corner(3));
		divCp[15] = blossomBicubicBezier(ctrlp, uvB.corner(3), uvB.corner(3), uvB.corner(3));
		return divCp;
	}
};
