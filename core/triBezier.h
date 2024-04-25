# pragma once
#include"paramBound.h"

class TriLinearBezier{
public:
	static const int order = 1;
	static const int cntCp = 3;
	std::array<Vector3d, 3> ctrlp;
	TriLinearBezier() {}
	TriLinearBezier(const std::array<Vector3d, 3>& p): ctrlp(p) {}

	Vector3d blossomBilinearBezier(const std::array<Vector3d, 3> &b, const BaryCoord& coord) const {
		if(std::abs(coord.u + coord.v + coord.w - 1) > 1e-12){
			std::cerr<<coord.u + coord.v + coord.w<<"param coord wrong!\n";
			exit(-1);
		}
		return coord.w * b[0] + coord.u * b[1] + coord.v * b[2];
	}

	Vector3d evaluatePatchPoint(const Array2d& uv) const {
		BaryCoord coord(uv);
		return blossomBilinearBezier(ctrlp, coord);
	}

	double feasibleUpperV(const double& u) const { return 1 - u; }
	static Vector3d axisU(const std::array<Vector3d, 3>& pt) {
		return pt[1]-pt[0];
	}
	static Vector3d axisV(const std::array<Vector3d, 3>& pt) {
		return pt[2]-pt[0];
	}

	// 100,010,001
	std::array<Vector3d, 3> divideBezierPatch(const TriParamBound& coords) const {
		std::array<Vector3d, 3> divCp;
		divCp[0] = blossomBilinearBezier(ctrlp, coords.nodes[2]);
		divCp[1] = blossomBilinearBezier(ctrlp, coords.nodes[0]);
		divCp[2] = blossomBilinearBezier(ctrlp, coords.nodes[1]);
		return divCp;
	}
};

class TriQuadBezier{
public:
	static const int order = 2;
	static const int cntCp = 6;
	std::array<Vector3d, 6> ctrlp;
	TriQuadBezier() {}
	TriQuadBezier(const std::array<Vector3d, 6>& p): ctrlp(p) {}

	Vector3d triLerp(const Vector3d& b0, const Vector3d& b1, const Vector3d& b2, const BaryCoord& coord) const {
		if(std::abs(coord.u + coord.v + coord.w - 1) > 1e-12){
			std::cerr<<coord.u + coord.v + coord.w<<"param coord wrong!\n";
			exit(-1);
		}
		return coord.w * b0 + coord.u * b1 + coord.v * b2;
	}

	Vector3d blossomBiquadBezier(const std::array<Vector3d, 6>& a, const BaryCoord& coord0, const BaryCoord& coord1) const {
		Vector3d b[3] = { triLerp(a[0], a[1], a[3], coord0), triLerp(a[1], a[2], a[4], coord0), triLerp(a[3], a[4], a[5], coord0) };
		return triLerp(b[0], b[1], b[2], coord1);
	}

	Vector3d evaluatePatchPoint(const Array2d& uv) const {
		BaryCoord coord(uv);
		return blossomBiquadBezier(ctrlp, coord, coord);
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
		divCp[0] = blossomBiquadBezier(ctrlp, coords.nodes[2], coords.nodes[2]);
		divCp[1] = blossomBiquadBezier(ctrlp, coords.nodes[0], coords.nodes[2]);
		divCp[2] = blossomBiquadBezier(ctrlp, coords.nodes[0], coords.nodes[0]);
		divCp[3] = blossomBiquadBezier(ctrlp, coords.nodes[1], coords.nodes[2]);
		divCp[4] = blossomBiquadBezier(ctrlp, coords.nodes[0], coords.nodes[1]);
		divCp[5] = blossomBiquadBezier(ctrlp, coords.nodes[1], coords.nodes[1]);
		return divCp;
	}
};

class TriCubicBezier {
public:
	static const int order = 3;
	static const int cntCp = 10;
	// 003,102,201,300,012,111,210,021,120,030
	std::array<Vector3d, 10> ctrlp;

	TriCubicBezier() {}
	TriCubicBezier(int randSeed){
		if(randSeed<0) std::srand(std::time(nullptr));
		else std::srand(randSeed);
		for (int i = 0; i < 10; i++)
			ctrlp[i] = Vector3d::Random();
	}

	TriCubicBezier(const std::array<Vector3d, 10>& p): ctrlp(p) {}

	Vector3d triLerp(const Vector3d& b0, const Vector3d& b1, const Vector3d& b2, const BaryCoord& coord) const {
		if(std::abs(coord.u + coord.v + coord.w - 1) > 1e-12){
			std::cerr<<coord.u + coord.v + coord.w<<"param coord wrong!\n";
			exit(-1);
		}
		return coord.w * b0 + coord.u * b1 + coord.v * b2;
	}

	Vector3d blossomBiquadBezier(const std::array<Vector3d, 6>& a, const BaryCoord& coord0, const BaryCoord& coord1) const {
		Vector3d b[3] = { triLerp(a[0], a[1], a[3], coord0), triLerp(a[1], a[2], a[4], coord0), triLerp(a[3], a[4], a[5], coord0) };
		return triLerp(b[0], b[1], b[2], coord1);
	}

	Vector3d blossomBicubicBezier(const std::array<Vector3d, 10>& p, const BaryCoord& coord0, const BaryCoord& coord1, const BaryCoord& coord2) const {
		std::array<Vector3d,6> a = { triLerp(p[0], p[1], p[4], coord0), triLerp(p[1], p[2], p[5], coord0), triLerp(p[2], p[3], p[6], coord0),
						 triLerp(p[4], p[5], p[7], coord0), triLerp(p[5], p[6], p[8], coord0), triLerp(p[7], p[8], p[9], coord0) };
		return blossomBiquadBezier(a, coord1, coord2);
	}

	Vector3d evaluatePatchPoint(const Array2d& uv) const {
		BaryCoord coord(uv);
		return blossomBicubicBezier(ctrlp, coord, coord, coord);
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
		divCp[0] = blossomBicubicBezier(ctrlp, coords.nodes[2], coords.nodes[2], coords.nodes[2]);
		divCp[1] = blossomBicubicBezier(ctrlp, coords.nodes[0], coords.nodes[2], coords.nodes[2]);
		divCp[2] = blossomBicubicBezier(ctrlp, coords.nodes[0], coords.nodes[0], coords.nodes[2]);
		divCp[3] = blossomBicubicBezier(ctrlp, coords.nodes[0], coords.nodes[0], coords.nodes[0]);
		divCp[4] = blossomBicubicBezier(ctrlp, coords.nodes[1], coords.nodes[2], coords.nodes[2]);
		divCp[5] = blossomBicubicBezier(ctrlp, coords.nodes[0], coords.nodes[1], coords.nodes[2]);
		divCp[6] = blossomBicubicBezier(ctrlp, coords.nodes[0], coords.nodes[0], coords.nodes[1]);
		divCp[7] = blossomBicubicBezier(ctrlp, coords.nodes[1], coords.nodes[1], coords.nodes[2]);
		divCp[8] = blossomBicubicBezier(ctrlp, coords.nodes[0], coords.nodes[1], coords.nodes[1]);
		divCp[9] = blossomBicubicBezier(ctrlp, coords.nodes[1], coords.nodes[1], coords.nodes[1]);
		return divCp;
	}

	void checkTri() const {
		for(int i=0;i<10;i++){
			std::cout<<"cp "<<i<<": "<<ctrlp[i].transpose()<<"\n";
		}
	}
};

