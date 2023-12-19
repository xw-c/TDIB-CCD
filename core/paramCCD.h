# pragma once
#include "triBezier.h"
#include "config.h"

class ParamCoord{
	double u, v;
	ParamCoord(const double _u, const double _v): u(_u), v(_v) {}
	ParamCoord(const ParamCoord& bc): u(bc.u), v(bc.v) {}
	ParamCoord operator*(const double scalar) const {
        return {u * scalar, v * scalar};
    }
	ParamCoord operator+(const ParamCoord& bc) const {
        return {u + bc.u, v + bc.v};
    }
	friend std::ostream& operator<<(std::ostream& os, const ParamCoord& coord) {
        os << "(" << coord.u << ", " << coord.v << ")";
        return os;
    }
};

class ParamBound{
public:
	ParamBound() {}
	virtual void initialize() = 0;
	virtual double diameter() const = 0;
	virtual ParamCoord centerParam() const = 0;
	virtual ParamBound interpSubpatchParam(const int id) const = 0;
	virtual ~ParamBound () = default;
};

//并不知道这个类的实现靠不靠谱
class ParamPatchPair{
	ParamBound *pb1, *pb2;
    double tLower;
    ParamPatchPair(const ParamBound* c1, const ParamBound* c2, double t = std::numeric_limits<double>::infinity()): 
		tLower(t) { *pb1=*c1, *pb2=*c2; }
	bool operator<(ParamPatchPair const &o) const { return tLower > o.tLower; }
	virtual ~ParamPatchPair() { delete pb1, pb2; }
};

template<const PatchType cntCp>
class ParamObj{
public:
	std::array<Vector3d, cntCp> ctrlp {};

	ParamObj() {}
	// ParamObj(const std::array<Vector3d, cntCp>& p): ctrlp(p) {}

	virtual Vector3d blossomBicubicBezier(const ParamCoord& coord) = 0;
	virtual std::array<Vector3d, cntCp> divideBezierPatch(const ParamBound& coords) = 0;
	virtual ~ParamBound () = default;
};

template<const PatchType cntCp>
class ParamCCD{
	ParamCoord uv1,uv2;
	const BoundingBoxType bbtype;
	ParamObj<cntCp> *CpPos1, *CpPos2, *CpVel1, *CpVel2;
	
	ParamCCD(const BoundingBoxType _bbtype);
	static double PrimitiveCheck(ParamBound const &divUvB1, ParamBound const &divUvB2);
	static double CCD();
	~ParamCCD() {delete CpPos1, CpPos2, CpVel1, CpVel2;}
};