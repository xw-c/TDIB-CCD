# pragma once
#include <Eigen/Dense>

using Eigen::Array2d;
using Eigen::Vector3d;
using Eigen::Vector4d;

#include <iostream>
#include <span>
#include "recBezier.h"

class RecQuadRatBezier{
public:
	static const int cntCp = 9;
	std::array<Vector4d, 9> ctrlp;
	RecQuadRatBezier(){}
	RecQuadRatBezier(const std::array<Vector3d, 9>& pos, const std::array<double, 9>& weight){
		for(int i = 0; i < 9; i++){
			ctrlp[i] = Vector4d(pos[i][0]*weight[i], pos[i][1]*weight[i], pos[i][2]*weight[i], weight[i]);
		}
	}
	RecQuadRatBezier(const std::array<double, 9>& weight){
		for(int i = 0; i < 9; i++){
			ctrlp[i] = Vector4d(0, 0, 0, weight[i]);
		}
	}
	Vector4d lerp(double t, Vector4d const &t0, Vector4d const &t1) const { return (1 - t) * t0 + t * t1; }
	Vector4d blossomQuadBezier(std::span<Vector4d const> p, double u0, double u1) const {
		Vector4d b[2] = { lerp(u0, p[0], p[1]), lerp(u0, p[1], p[2]) };
		return lerp(u1, b[0], b[1]);
	}

	Vector4d blossomBiquadBezier(std::span<Vector4d const> cp, Array2d const &uv0, Array2d const &uv1) const {
		std::array<Vector4d, 3> q;
		for (int i = 0; i < 3; i++) {
			q[i] = blossomQuadBezier(cp.subspan(i * 3, 3), uv0.y(), uv1.y());
		}
		return blossomQuadBezier(q, uv0.x(), uv1.x());
	}

	Vector3d get3dPos(const Vector4d& pt) const { return Vector3d(pt[0]/pt[3], pt[1]/pt[3], pt[2]/pt[3]); }
	Vector3d evaluatePatchPoint(Array2d const &uv) const {
		Vector4d pt = blossomBiquadBezier(ctrlp, uv, uv);
		if(pt[3]==0){
			std::cout<<"why zero weight?\n";exit(-1);
		}
		return get3dPos(pt);
	}

	static double feasibleUpperV(const double& u) { return 1; }
	static Vector3d axisU(const std::array<Vector3d, 9>& pt){
		return (pt[6]-pt[0]+pt[8]-pt[2]);
	}
	static Vector3d axisV(const std::array<Vector3d, 9>& pt){
		return (pt[2]-pt[0]+pt[8]-pt[6]);
	}

	// Patch Functions
	// 先v变再u变
	std::array<Vector3d, 9> divideBezierPatch(RecParamBound const &uvB) const {
		std::array<Vector3d, 9> divCp;
		divCp[0] = get3dPos(blossomBiquadBezier(ctrlp, uvB.corner(0), uvB.corner(0)));
		divCp[1] = get3dPos(blossomBiquadBezier(ctrlp, uvB.corner(0), uvB.corner(2)));
		divCp[2] = get3dPos(blossomBiquadBezier(ctrlp, uvB.corner(2), uvB.corner(2)));
		divCp[3] = get3dPos(blossomBiquadBezier(ctrlp, uvB.corner(0), uvB.corner(1)));
		divCp[4] = get3dPos(blossomBiquadBezier(ctrlp, uvB.corner(0), uvB.corner(3)));
		divCp[5] = get3dPos(blossomBiquadBezier(ctrlp, uvB.corner(2), uvB.corner(3)));
		divCp[6] = get3dPos(blossomBiquadBezier(ctrlp, uvB.corner(1), uvB.corner(1)));
		divCp[7] = get3dPos(blossomBiquadBezier(ctrlp, uvB.corner(1), uvB.corner(3)));
		divCp[8] = get3dPos(blossomBiquadBezier(ctrlp, uvB.corner(3), uvB.corner(3)));
		return divCp;
	}
};

// Bicubic Bezier Functions
class RecCubicRatBezier{
public:
	static const int cntCp = 16;
	// wx, wy, wz, w
	std::array<Vector4d, 16> ctrlp;

	// RecCubicRatBezier(const std::array<Vector3d, 16>& pos, const std::array<double, 16>& weight){
	// 	for(int i = 0; i < 16; i++){
	// 		ctrlp[i] = Vector4d(pos[i][0]*weight[i], pos[i][1]*weight[i], pos[i][2]*weight[i], weight[i]);
	// 	}
	// }
	// RecCubicRatBezier(const std::array<double, 16>& weight){
	// 	for(int i = 0; i < 16; i++){
	// 		ctrlp[i] = Vector4d(0, 0, 0, weight[i]);
	// 	}
	// }

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
	static Vector3d axisU(const std::array<Vector3d, 16>& pt){
		return (pt[12]-pt[0]+pt[15]-pt[3]);
	}
	static Vector3d axisV(const std::array<Vector3d, 16>& pt){
		return (pt[3]-pt[0]+pt[15]-pt[12]);
	}
	// Patch Functions
	// 先v变再u变
	std::array<Vector3d, 16> divideBezierPatch(RecParamBound const &uvB) const {
		std::array<Vector3d, 16> divCp;
		divCp[0] = get3dPos(blossomBicubicBezier(ctrlp, uvB.corner(0), uvB.corner(0), uvB.corner(0)));
		divCp[1] = get3dPos(blossomBicubicBezier(ctrlp, uvB.corner(0), uvB.corner(0), uvB.corner(2)));
		divCp[2] = get3dPos(blossomBicubicBezier(ctrlp, uvB.corner(0), uvB.corner(2), uvB.corner(2)));
		divCp[3] = get3dPos(blossomBicubicBezier(ctrlp, uvB.corner(2), uvB.corner(2), uvB.corner(2)));
		divCp[4] = get3dPos(blossomBicubicBezier(ctrlp, uvB.corner(0), uvB.corner(0), uvB.corner(1)));
		divCp[5] = get3dPos(blossomBicubicBezier(ctrlp, uvB.corner(0), uvB.corner(0), uvB.corner(3)));
		divCp[6] = get3dPos(blossomBicubicBezier(ctrlp, uvB.corner(0), uvB.corner(2), uvB.corner(3)));
		divCp[7] = get3dPos(blossomBicubicBezier(ctrlp, uvB.corner(2), uvB.corner(2), uvB.corner(3)));
		divCp[8] = get3dPos(blossomBicubicBezier(ctrlp, uvB.corner(0), uvB.corner(1), uvB.corner(1)));
		divCp[9] = get3dPos(blossomBicubicBezier(ctrlp, uvB.corner(0), uvB.corner(1), uvB.corner(3)));
		divCp[10] = get3dPos(blossomBicubicBezier(ctrlp, uvB.corner(0), uvB.corner(3), uvB.corner(3)));
		divCp[11] = get3dPos(blossomBicubicBezier(ctrlp, uvB.corner(2), uvB.corner(3), uvB.corner(3)));
		divCp[12] = get3dPos(blossomBicubicBezier(ctrlp, uvB.corner(1), uvB.corner(1), uvB.corner(1)));
		divCp[13] = get3dPos(blossomBicubicBezier(ctrlp, uvB.corner(1), uvB.corner(1), uvB.corner(3)));
		divCp[14] = get3dPos(blossomBicubicBezier(ctrlp, uvB.corner(1), uvB.corner(3), uvB.corner(3)));
		divCp[15] = get3dPos(blossomBicubicBezier(ctrlp, uvB.corner(3), uvB.corner(3), uvB.corner(3)));
		return divCp;
	}
};
