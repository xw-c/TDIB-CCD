# pragma once
#include "mathOps.h"
#include "triBezier.h"
#include "recBezier.h"
static double point2BoxIntersect(const MatrixXd &x, const MatrixXd &u, const double tMax)
{
	int groups = int(x.rows()), items = int(x.cols());
	std::vector<Array2d> invlLists[2*x.rows()], finalList;
	for (int num = 0; num < groups; num++) {
		double left = 0. - Epsilon, right = tMax + Epsilon;
		for (int i = 0; i < items; i++) {
			if (u(num, i) > Epsilon) left = std::max(left, -x(num, i) / u(num, i));
			else if (-u(num, i) > Epsilon) right = std::min(right, -x(num, i) / u(num, i));
			else if (x(num, i) < Epsilon) { left = tMax; break; }
		}
		invlLists[num].clear();
		if (left >= right || left >= tMax || right <= 0.) invlLists[num].push_back(Array2d(0., tMax));
		else {
			if (left >= 0.) invlLists[num].push_back(Array2d(0., left));
			if (right <= tMax) invlLists[num].push_back(Array2d(right, tMax));
		}
	}

	auto intervalIntersection = [&] (const std::vector<Array2d> &list1, const std::vector<Array2d> &list2, std::vector<Array2d> &res)
	{
		res.clear();
		for (int idx1 = 0, idx2 = 0; idx1 < list1.size() && idx2 < list2.size();) {
			double start1 = list1[idx1][0], end1 = list1[idx1][1];
			double start2 = list2[idx2][0], end2 = list2[idx2][1];
			if (end1 >= start2 && end2 >= start1)
				res.push_back(Array2d(std::max(start1, start2), std::min(end1, end2)));
			if (end1 < end2) idx1++;
			else idx2++;
		}
	};

	finalList.clear();
	for (int num = 1; num < groups; num++) {
		intervalIntersection(invlLists[num - 1], invlLists[num], finalList);
		invlLists[num].clear();
		for (auto &item : finalList) invlLists[num].push_back(item);
	}
	if (finalList.empty())
		return std::numeric_limits<double>::infinity();
	return finalList[0][0];
}

template <typename ParamBound, typename Func>
static double greedyIntersect(Func &&checkIntersect, Array2d &intsctUV)
{
	struct DividedPatch{
		ParamBound pb;
		double tLower;
		DividedPatch(const ParamBound& c, double t = std::numeric_limits<double>::infinity()): pb(c), tLower(t) {}
		bool operator<(DividedPatch const &o) const { return tLower > o.tLower; }
	};


	// Initialize the heap
	std::priority_queue<DividedPatch> heap;
	{
		ParamBound uvB;
		const double tTent = checkIntersect(uvB);
		if (tTent < DeltaT)
			heap.emplace(uvB, tTent);
	}

	// Greedily search intersection points
	while (!heap.empty()) {
		const auto cur = heap.top();
		heap.pop();

		// Set uv of the middle point
		// Array2d uvMid = (cur.uvB.pMin + cur.uvB.pMax) / 2;

		// Decide whether the algorithm converges
		if (cur.pb.diameter() < MinDeltaUV) {
			intsctUV = cur.pb.centerParam();
			return cur.tLower;
		}

		// Divide the current patch into four pieces
		for (int i = 0; i < 4; ++i) {
			ParamBound divUvB(cur.pb.interpSubpatchParam(i));

			const double tTent = checkIntersect(divUvB);
			if (tTent < std::numeric_limits<double>::infinity())
				heap.emplace(divUvB, tTent);
		}
	}

	return -1;
}

template <typename ParamObj2, typename ParamBound2>
static double point2patch(const Vector3d& p1, const Vector3d& v1, 
						const ParamObj2 &CpPos2, const ParamObj2 &CpVel2,
						Array2d &intsctUV, const BoundingBoxType bbtype) {

	auto checkIntersection = [&](const ParamBound2 &b) {
		const auto cpPosB = CpPos2.divideBezierPatch(b);
		const auto cpVelB = CpVel2.divideBezierPatch(b);
		// Get the coefficients of inequalities.

		auto setAxes = [&] (std::vector<Vector3d>& axes) {
			if(bbtype==BoundingBoxType::AABB){
				axes = {Vector3d::Unit(0), Vector3d::Unit(1), Vector3d::Unit(2)};
			}
			else if(bbtype==BoundingBoxType::DOP14){
				axes = {Vector3d::Unit(0), Vector3d::Unit(1), Vector3d::Unit(2), 
						Vector3d(1,1,1).normalized(), Vector3d(-1,1,1).normalized(), Vector3d(1,-1,1).normalized(), Vector3d(1,1,-1).normalized()};
			}
			else if(bbtype==BoundingBoxType::OBB){
				Vector3d lu = (cpPosB[9]-cpPosB[0]).normalized();
				Vector3d lv = (cpPosB[9]-cpPosB[3]+cpPosB[0]-cpPosB[3]);
				lv = (lv-lv.dot(lu)*lu).normalized().eval();
				Vector3d ln = lu.cross(lv).normalized();
				axes = {lu, lv, ln};
			}
		};

		std::vector<Vector3d> axes;
		axes.clear();
		setAxes(axes);
		const int axesSize = axes.size();

		MatrixXd x(2 * axesSize, 16), u(2 * axesSize, 16);
		double x0[2 * axesSize], u0[2 * axesSize];
		//改这里对应BB的dim
		for (int dim = 0; dim < axesSize; dim++)
			x0[dim] = p1.dot(axes[dim]), x0[axesSize + dim] = -p1.dot(axes[dim]),
			u0[dim] = v1.dot(axes[dim]), u0[axesSize + dim] = -v1.dot(axes[dim]);
		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 4; j++)
				for (int dim = 0; dim < axesSize; dim++)
					x(dim, i << 2 | j) = cpPosB[i << 2 | j].dot(axes[dim])- x0[dim],
					x(axesSize + dim, i << 2 | j) = -cpPosB[i << 2 | j].dot(axes[dim]) - x0[axesSize + dim],
					u(dim, i << 2 | j) = cpVelB[i << 2 | j].dot(axes[dim]) - u0[dim],
					u(axesSize + dim, i << 2 | j) = -cpVelB[i << 2 | j].dot(axes[dim]) - u0[axesSize + dim];

		// Perform ray bounding box intersection.
		return point2BoxIntersect(x, u, DeltaT);
	};

	return greedyIntersect<ParamBound2>(std::move(checkIntersection), intsctUV);
}

template<typename ParamObj1, typename ParamObj2, typename ParamBound2>
static double sampleCCD(const ParamObj1 &CpPos1, const ParamObj1 &CpVel1, 
						const ParamObj2 &CpPos2, const ParamObj2 &CpVel2,
						Array2d& uv1, Array2d& uv2, 
						const BoundingBoxType bbtype) {
	using steady_clock = std::chrono::steady_clock;
	using duration = std::chrono::duration<double>;
	const auto initialTime = steady_clock::now();

	constexpr double eps = .01;
	double minT = std::numeric_limits<double>::infinity();
	Array2d ap1, ap2, temptUV;

	for (double x1 = 0; x1 < 1; x1 += eps) {
		for (double y1 = 0; y1 < ParamObj1::feasibleUpperV(x1); y1 += eps) {
			Vector3d const p1 = CpPos1.evaluatePatchPoint({ x1, y1 });
			Vector3d const v1 = CpVel1.evaluatePatchPoint({ x1, y1 });
			
			double temptT = point2patch<ParamObj2, ParamBound2>(p1, v1, CpPos2, CpVel2, temptUV, bbtype);
			if(temptT >= 0 && temptT <= minT) minT=temptT, ap1={ x1, y1 }, ap2=temptUV;
		}
	}
	if(minT<=DeltaT&&minT>=0){
		std::cout << ap1.transpose() << std::endl << ap2.transpose() << std::endl;
		const auto endTime = steady_clock::now();
		uv1 = ap1,
		uv2 = ap2;
		std::cout << "min time: "<<  minT << "used seconds: " <<
			duration(endTime - initialTime).count()
			<< std::endl;
		return minT;
	}
	else{
		const auto endTime = steady_clock::now();
		std::cout << "used seconds: " <<
			duration(endTime - initialTime).count()
			<< std::endl;
		return -1;
	}
}

//得改一下对三角形的采样
auto recSampleCCD = sampleCCD<RecBezierObj, RecBezierObj, RecParamBound>;