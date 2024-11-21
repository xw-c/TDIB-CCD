# pragma once
#include <Eigen/Dense>

using Eigen::Array2d;
using Eigen::Vector3d;
using Eigen::Vector4d;
using Eigen::MatrixXd;

#include <iostream>
#include <span>
#include <array>

// For rec patches
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
	double width() const {
		return (pMax - pMin).maxCoeff();
	}
};

// For tri patches
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
};

// For tri patches
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
	double width() const {
		double uMax=std::max(std::max(nodes[0].u, nodes[1].u), nodes[2].u), uMin=std::min(std::min(nodes[0].u, nodes[1].u), nodes[2].u); 
		double vMax=std::max(std::max(nodes[0].v, nodes[1].v), nodes[2].v), vMin=std::min(std::min(nodes[0].v, nodes[1].v), nodes[2].v); 
		return std::max(uMax-uMin, vMax-vMin);
	}
};