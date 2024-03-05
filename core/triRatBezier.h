# pragma once
#include"paramBound.h"

class TriQuadRatBezier{
public:
	static const int cntCp = 6;
	std::array<Vector4d, 6> ctrlp;
	TriQuadRatBezier() {}
	TriQuadRatBezier(const std::array<Vector3d, 6>& pos, const std::array<double, 6>& weight){
		for(int i = 0; i < 6; i++){
			if(weight[i]==0){
				std::cerr<<"zero weight!\n";
				exit(-1);
			}
			ctrlp[i] = Vector4d(pos[i][0]*weight[i], pos[i][1]*weight[i], pos[i][2]*weight[i], weight[i]);
		}
	}
	TriQuadRatBezier(const std::array<double, 6>& weight){
		for(int i = 0; i < 6; i++){
			if(weight[i]==0){
				std::cerr<<"zero weight!\n";
				exit(-1);
			}
			ctrlp[i] = Vector4d(0, 0, 0, weight[i]);
		}
	}

	Vector4d triLerp(const Vector4d& b0, const Vector4d& b1, const Vector4d& b2, const BaryCoord& coord) const {
		if(std::abs(coord.u + coord.v + coord.w - 1) > 1e-12){
			std::cerr<<coord.u + coord.v + coord.w<<"param coord wrong!\n";
			exit(-1);
		}
		return coord.w * b0 + coord.u * b1 + coord.v * b2;
	}

	Vector4d blossomBiquadBezier(const std::array<Vector4d, 6>& a, const BaryCoord& coord0, const BaryCoord& coord1) const {
		Vector4d b[3] = { triLerp(a[0], a[1], a[3], coord0), triLerp(a[1], a[2], a[4], coord0), triLerp(a[3], a[4], a[5], coord0) };
		return triLerp(b[0], b[1], b[2], coord1);
	}

	Vector3d get3dPos(const Vector4d& pt) const { return Vector3d(pt[0]/pt[3], pt[1]/pt[3], pt[2]/pt[3]); }
	Vector3d evaluatePatchPoint(const Array2d& uv) const {
		BaryCoord coord(uv);
		return get3dPos(blossomBiquadBezier(ctrlp, coord, coord));
	}

	double feasibleUpperV(const double& u) const { return 1 - u; }
	static Vector3d axisU(const std::array<Vector3d, 6>& pt) {
		return pt[2]-pt[0];
	}
	static Vector3d axisV(const std::array<Vector3d, 6>& pt) {
		return pt[5]-pt[0];
	}

	// 100,010,001
	std::array<Vector3d, 6> divideBezierPatch(const TriParamBound& coords) const {
		std::array<Vector3d, 6> divCp;
		divCp[0] = get3dPos(blossomBiquadBezier(ctrlp, coords.nodes[2], coords.nodes[2]));
		divCp[1] = get3dPos(blossomBiquadBezier(ctrlp, coords.nodes[0], coords.nodes[2]));
		divCp[2] = get3dPos(blossomBiquadBezier(ctrlp, coords.nodes[0], coords.nodes[0]));
		divCp[3] = get3dPos(blossomBiquadBezier(ctrlp, coords.nodes[1], coords.nodes[2]));
		divCp[4] = get3dPos(blossomBiquadBezier(ctrlp, coords.nodes[0], coords.nodes[1]));
		divCp[5] = get3dPos(blossomBiquadBezier(ctrlp, coords.nodes[1], coords.nodes[1]));
		return divCp;
	}
};

class TriCubicRatBezier {
public:
	static const int cntCp = 10;
	// 003,102,201,300,012,111,210,021,120,030
	std::array<Vector4d, 10> ctrlp;

	TriCubicRatBezier() {}
	TriCubicRatBezier(const std::array<Vector3d, 10>& pos, const std::array<double, 10>& weight){
		for(int i = 0; i < 10; i++){
			if(weight[i]==0){
				std::cerr<<"zero weight!\n";
				exit(-1);
			}
			ctrlp[i] = Vector4d(pos[i][0]*weight[i], pos[i][1]*weight[i], pos[i][2]*weight[i], weight[i]);
		}
	}
	TriCubicRatBezier(const std::array<double, 10>& weight){
		for(int i = 0; i < 10; i++){
			if(weight[i]==0){
				std::cerr<<"zero weight!\n";
				exit(-1);
			}
			ctrlp[i] = Vector4d(0, 0, 0, weight[i]);
		}
	}

	TriCubicRatBezier(const std::array<Vector4d, 10>& p): ctrlp(p) {}

	Vector4d triLerp(const Vector4d& b0, const Vector4d& b1, const Vector4d& b2, const BaryCoord& coord) const {
		if(std::abs(coord.u + coord.v + coord.w - 1) > 1e-12){
			std::cerr<<coord.u + coord.v + coord.w<<"param coord wrong!\n";
			exit(-1);
		}
		return coord.w * b0 + coord.u * b1 + coord.v * b2;
	}

	Vector4d blossomBiquadBezier(const std::array<Vector4d, 6>& a, const BaryCoord& coord0, const BaryCoord& coord1) const {
		Vector4d b[3] = { triLerp(a[0], a[1], a[3], coord0), triLerp(a[1], a[2], a[4], coord0), triLerp(a[3], a[4], a[5], coord0) };
		return triLerp(b[0], b[1], b[2], coord1);
	}

	Vector4d blossomBicubicBezier(const std::array<Vector4d, 10>& p, const BaryCoord& coord0, const BaryCoord& coord1, const BaryCoord& coord2) const {
		std::array<Vector4d,6> a = { triLerp(p[0], p[1], p[4], coord0), triLerp(p[1], p[2], p[5], coord0), triLerp(p[2], p[3], p[6], coord0),
						 triLerp(p[4], p[5], p[7], coord0), triLerp(p[5], p[6], p[8], coord0), triLerp(p[7], p[8], p[9], coord0) };
		return blossomBiquadBezier(a, coord1, coord2);
	}

	Vector3d get3dPos(const Vector4d& pt) const { return Vector3d(pt[0]/pt[3], pt[1]/pt[3], pt[2]/pt[3]); }
	Vector3d evaluatePatchPoint(const Array2d& uv) const {
		BaryCoord coord(uv);
		return get3dPos(blossomBicubicBezier(ctrlp, coord, coord, coord));
	}

	double feasibleUpperV(const double& u) const { return 1 - u; }
	static Vector3d axisU(const std::array<Vector3d, 10>& pt) {
		return pt[3]-pt[0];
	}
	static Vector3d axisV(const std::array<Vector3d, 10>& pt) {
		return pt[9]-pt[0];
	}

	// 100,010,001
	std::array<Vector3d, 10> divideBezierPatch(const TriParamBound& coords) const {
		std::array<Vector3d, 10> divCp;
		divCp[0] = get3dPos(blossomBicubicBezier(ctrlp, coords.nodes[2], coords.nodes[2], coords.nodes[2]));
		divCp[1] = get3dPos(blossomBicubicBezier(ctrlp, coords.nodes[0], coords.nodes[2], coords.nodes[2]));
		divCp[2] = get3dPos(blossomBicubicBezier(ctrlp, coords.nodes[0], coords.nodes[0], coords.nodes[2]));
		divCp[3] = get3dPos(blossomBicubicBezier(ctrlp, coords.nodes[0], coords.nodes[0], coords.nodes[0]));
		divCp[4] = get3dPos(blossomBicubicBezier(ctrlp, coords.nodes[1], coords.nodes[2], coords.nodes[2]));
		divCp[5] = get3dPos(blossomBicubicBezier(ctrlp, coords.nodes[0], coords.nodes[1], coords.nodes[2]));
		divCp[6] = get3dPos(blossomBicubicBezier(ctrlp, coords.nodes[0], coords.nodes[0], coords.nodes[1]));
		divCp[7] = get3dPos(blossomBicubicBezier(ctrlp, coords.nodes[1], coords.nodes[1], coords.nodes[2]));
		divCp[8] = get3dPos(blossomBicubicBezier(ctrlp, coords.nodes[0], coords.nodes[1], coords.nodes[1]));
		divCp[9] = get3dPos(blossomBicubicBezier(ctrlp, coords.nodes[1], coords.nodes[1], coords.nodes[1]));
		return divCp;
	}
};

