# pragma once
#include " config.h"

template<typename ParamCoord>
class ParamBound{
public:
	ParamBound() {}
	virtual void initialize() = 0;
	virtual double diameter() const = 0;
	virtual ParamCoord centerParam() const = 0;
	virtual ParamBound interpSubpatchParam(const int id) const = 0;
	virtual ~ParamBound () = default;
};

template<const PatchType cntCp, typename ParamCoord>
class ParamObj{
public:
	std::array<Vector3d, cntCp> ctrlp {};

	ParamObj() {}
	// ParamObj(const std::array<Vector3d, cntCp>& p): ctrlp(p) {}

	virtual Vector3d evaluatePatchPoint(const ParamCoord& coord) = 0;
	virtual std::array<Vector3d, cntCp> divideBezierPatch(const ParamBound& coords) = 0;
	virtual ~ParamObj () = default;
};