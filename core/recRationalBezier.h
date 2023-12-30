# pragma once
#include <Eigen/Dense>

using Eigen::Array2d;
using Eigen::Vector3d;
using Eigen::Vector4d;

#include <iostream>
#include <span>
#include "recBezier.h"

// Bicubic Bezier Functions
class RecRatBezierObj{
public:
	static const int cntCp = 16;
	// wx, wy, wz, w
	std::array<Vector4d, 16> ctrlp;

	Vector4d lerp(double t, Vector4d const &t0, Vector4d const &t1) const { return (1 - t) * t0 + t * t1; }
	
	Vector4d blossomCubicBezier(std::span<Vector4d const> p, double u0, double u1, double u2) const {
		Vector4d a[3] = { lerp(u0, p[0], p[1]), lerp(u0, p[1], p[2]), lerp(u0, p[2], p[3]) };
		Vector4d b[2] = { lerp(u1, a[0], a[1]), lerp(u1, a[1], a[2]) };
		return lerp(u2, b[0], b[1]);
	}

	Vector4d blossomBicubicBezier(std::span<Vector4d const> cp, Array2d const &uv0, Array2d const &uv1, Array2d const &uv2) const {
		std::array<Vector4d, 4> q;
		for (int i = 0; i < 4; i++) {
			q[i] = blossomCubicBezier(cp.subspan(i * 4, 4), uv0.y(), uv1.y(), uv2.y());
		}
		return blossomCubicBezier(q, uv0.x(), uv1.x(), uv2.x());
	}
	Vector3d get3dPos(const Vector4d& pt) const { return Vector3d(pt[0]/pt[3], pt[1]/pt[3], pt[2]/pt[3]); }
	Vector3d evaluatePatchPoint(Array2d const &uv) const {
		Vector4d pt = blossomBicubicBezier(ctrlp, uv, uv, uv);
		if(pt[3]==0){
			std::cout<<"why zero weight?\n";exit(-1);
		}
		return get3dPos(pt);
	}

	static double feasibleUpperV(const double& u) { return 1; }

	// Patch Functions
	// 先v变再u变
	std::array<Vector4d, 16> divideBezierPatch(RecParamBound const &uvB) const {
		std::array<Vector4d, 16> divCp;
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

void generateRationalPatches(RecRatBezierObj &CpPos1, RecRatBezierObj &CpPos2, RecRatBezierObj &CpVel1, RecRatBezierObj &CpVel2){
	std::srand(10);
	for (int i = 0; i < 16; i++) {
		CpPos1.ctrlp[i] = Vector4d::Random() - Vector4d::Constant(.6);
		CpVel1.ctrlp[i] = Vector4d::Random()*0.3 + Vector4d::Constant(.6);
		CpPos2.ctrlp[i] = Vector4d::Random() + Vector4d::Constant(.6);
		CpVel2.ctrlp[i] = Vector4d::Random()*0.3 - Vector4d::Constant(.6);
	}
}