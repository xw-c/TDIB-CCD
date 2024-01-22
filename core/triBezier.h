# pragma once
#include <Eigen/Dense>

using Eigen::Array2d;
using Eigen::Vector3d;
using Eigen::MatrixXd;

#include<iostream>
#include<array>

class BaryCoord {
public:
	double u,v,w;
	BaryCoord() {}

	BaryCoord(const double _u, const double _v): u(_u), v(_v), w(1-_u-_v) {}
	BaryCoord(const double _u, const double _v, const double _w): u(_u), v(_v), w(1-_u-_v) {}
	BaryCoord(const Array2d& bc): u(bc[0]), v(bc[1]), w(1-bc[0]-bc[1]) {}
	BaryCoord(const BaryCoord& bc): u(bc.u), v(bc.v), w(1-bc.u-bc.v) {}
	// BaryCoord(const double u, const double v, const double w): u(u), v(v), w(w) {}
	// BaryCoord(const BaryCoord& bc): u(bc.u), v(bc.v), w(bc.w) {}
	BaryCoord operator*(const double scalar) const {
        return BaryCoord(u * scalar, v * scalar, w * scalar);
    }
	BaryCoord operator+(const BaryCoord& bc) const {
        return BaryCoord(u + bc.u, v + bc.v, w + bc.w);
    }
	friend std::ostream& operator<<(std::ostream& os, const BaryCoord& coord) {
        os << "(" << coord.u << ", " << coord.v << ", " << coord.w << ")";
        return os;
    }
	double computeSquaredDist(const BaryCoord& bc) const {
		return (u - bc.u)*(u - bc.u)+ (v - bc.v)*(v - bc.v) + (w - bc.w)*(w - bc.w);
	}
};

class TriParamBound{
public:

    std::array<BaryCoord,3> nodes;

	inline static const std::array<BaryCoord,3> relativeSubpatch[4] = {
		{BaryCoord(0.5,0,0.5), BaryCoord(0,0.5,0.5), BaryCoord(0,0,1)}, 
		{BaryCoord(1,0,0), BaryCoord(0.5,0.5,0), BaryCoord(0.5,0,0.5)}, 
		{BaryCoord(0.5,0.5,0), BaryCoord(0,1,0), BaryCoord(0,0.5,0.5)}, 
		{BaryCoord(0,0.5,0.5), BaryCoord(0.5,0,0.5), BaryCoord(0.5,0.5,0)}
	};

    TriParamBound(const BaryCoord& n0 = BaryCoord(1,0,0), const BaryCoord& n1 = BaryCoord(0,1,0), const BaryCoord& n2 = BaryCoord(0,0,1)) 
				{ nodes[0] = n0, nodes[1] = n1, nodes[2] = n2; }
    TriParamBound(const std::array<BaryCoord,3>& n): nodes(n) {}

	BaryCoord interpPointParam(const BaryCoord& interp) const { return nodes[0]*interp.u + nodes[1]*interp.v + nodes[2]*interp.w; }
	TriParamBound interpSubpatchParam(const int id) const {
		std::array<BaryCoord,3> subParam;
		for(int i=0; i<3; i++){
			subParam[i]=interpPointParam(relativeSubpatch[id][i]);
		}
		return TriParamBound(subParam);
	}

	Array2d centerParam() const { return Array2d((nodes[0].u+nodes[1].u+nodes[2].u)/3., (nodes[0].v+nodes[1].v+nodes[2].v)/3.); }
	double diameter() const {
		const double e[3]={nodes[0].computeSquaredDist(nodes[1]), nodes[2].computeSquaredDist(nodes[1]), nodes[0].computeSquaredDist(nodes[2])};
		return std::max(std::max(e[0], e[1]), e[2]);
	}
};
class TriLinearBezier{
public:
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

