# pragma once
#include <algorithm>
#include <array>
#include <string>
#include <cstdint>
#include <cfloat>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <limits>
#include <queue>
#include <set>
#include <span>
#include <chrono>
#include <cmath>
#include <random>
#include "../core/triBezier.h"

#include <Eigen/Dense>

using Eigen::Array2d;
using Eigen::Vector3d;
const double INFT = std::numeric_limits<double>::infinity();
const double PI = std::acos(-1);
enum class BoundingBoxType { AABB, OBB };

struct Line
{
	double k, b;
	Line():k(0),b(0){}
	void set(const double& _k, const double &_b){k=_k,b=_b;}
	Line(const double& k,const double& b): k(k), b(b) {}
	bool operator<(const Line &l) const { 
		return k < l.k || (k == l.k && b > l.b); // 相同斜率的直线中只有截距最大的被留下来
	}
	bool operator==(const Line &l) const {return k == l.k;}
};

// class Vert{
// public:
// 	static const int cntCp = 1;
// 	std::array<Vector3d,1> ctrlp;
// 	Vert() {}
// 	Vert(const Vector3d&v) {ctrlp[0]=v;}
// 	Vert(const std::array<Vector3d,1>&v): ctrlp(v) {}
// };

class Edge{
public:
	static const int cntCp = 2;
	std::array<Vector3d,2> ctrlp;	
	Edge() {}
	Edge(const Vector3d&e0, const Vector3d&e1) {ctrlp[0]=e0, ctrlp[1]=e1;}
	Edge(const std::array<Vector3d,2>&e): ctrlp(e) {}

	static Vector3d direction(const std::array<Vector3d,2>& p) {return p[1]-p[0];}
	Vector3d lerp(const double& u)const {return u*ctrlp[0]+(1-u)*ctrlp[1];}

	std::array<Vector3d, 2> divideBezierPatch(const Array2d& coords) const {
		std::array<Vector3d, 2> pts={lerp(coords[0]), lerp(coords[1])};
		return pts;
	}
};

class Face: public TriLinearBezier{
public:
	Face():TriLinearBezier(){}
	Face(const std::array<Vector3d, 3>& p): TriLinearBezier(p) {}
};

class EEPair{
public:
	Array2d pb1;
	Array2d pb2;
	Array2d tIntv;
	EEPair(const Array2d& c1, const Array2d& c2, const Array2d& t): pb1(c1), pb2(c2), tIntv(t) {}
	bool operator<(EEPair const &o) const { return tIntv[0] > o.tIntv[0]; }
	double calcWidth() const{
		const double w1 = pb1[1]-pb1[0], w2 = pb2[1]-pb2[0], wt = tIntv[1]-tIntv[0];
		return std::max(std::max(w1,w2), wt);
	}
	// double calcL1Dist(const Edge &CpPos1, const Edge &CpVel1, 
	// 				const Edge &CpPos2, const Edge &CpVel2) const{
	// 	auto ptPos1 = CpPos1.divideBezierPatch(pb1);
	// 	auto ptVel1 = CpVel1.divideBezierPatch(pb1);
	// 	auto ptPos2 = CpPos2.divideBezierPatch(pb2);
	// 	auto ptVel2 = CpVel2.divideBezierPatch(pb2);
	// 	for(int i=0;i<Edge::cntCp;i++)
	// 		ptPos1[i]+=ptVel1[i]*tIntv[0];
	// 	for(int i=0;i<Edge::cntCp;i++)
	// 		ptPos2[i]+=ptVel2[i]*tIntv[0];
	// 	double d1=(ptPos1[0].cwiseMax(ptPos1[1])-ptPos1[0].cwiseMin(ptPos1[1])).maxCoeff();
	// 	double d2=(ptPos2[0].cwiseMax(ptPos2[1])-ptPos2[0].cwiseMin(ptPos2[1])).maxCoeff();
	// 	return std::max(d1, d2);
	// }
};

struct VFPair{
	TriParamBound pb;
	Array2d tIntv;
	VFPair(const TriParamBound& c, const Array2d& t): pb(c), tIntv(t) {}
	bool operator<(VFPair const &o) const { return tIntv[0] > o.tIntv[0]; }
	double calcWidth() const{
		const double wpb = pb.width(), wt = tIntv[1]-tIntv[0];
		return std::max(wpb, wt);
	}
	// double calcL1Dist(const Vector3d &CpPos1, const Vector3d &CpVel1, 
	// 				const Face &CpPos2, const Face &CpVel2) const{
	// 	auto ptPos2 = CpPos2.divideBezierPatch(pb);
	// 	auto ptVel2 = CpVel2.divideBezierPatch(pb);
	// 	for(int i=0;i<Face::cntCp;i++)
	// 		ptPos2[i]+=ptVel2[i]*tIntv[0];
	// 	Vector3d projMax = (ptPos2[0].cwiseMax(ptPos2[1])).cwiseMax(ptPos2[2]);
	// 	Vector3d projMin = (ptPos2[0].cwiseMin(ptPos2[1])).cwiseMin(ptPos2[2]);
	// 	return (projMax-projMin).maxCoeff();
	// }
};
