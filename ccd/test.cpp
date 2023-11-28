#include <Eigen/Dense>

using Eigen::Array2d;
using Eigen::Vector3d;
using Eigen::MatrixXd;

#include <algorithm>
#include <array>
#include <string>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <limits>
#include <queue>
#include <span>
#include <chrono>
static constexpr double MinDeltaUV = 1e-6;
static constexpr double Epsilon = 1e-6;
static constexpr double DeltaT = 1;
bool DEBUG = 0;
static constexpr bool SHOWANS = 0;
enum class BB { AABB, OBB };

std::array<Vector3d, 16> CpPos1 {};
std::array<Vector3d, 16> CpPos2 {};
std::array<Vector3d, 16> CpVel1 {};
std::array<Vector3d, 16> CpVel2 {};

std::uint64_t cnt;

template <typename T> static T lerp(double t, T const &t0, T const &t1) { return (1 - t) * t0 + t * t1; }

struct Bounds2d {
	Array2d pMin, pMax;

	Bounds2d() = default;
	Bounds2d(Array2d const &p1, Array2d const &p2) :
		pMin { std::min(p1[0], p2[0]), std::min(p1[1], p2[1]) },
		pMax { std::max(p1[0], p2[0]), std::max(p1[1], p2[1]) } {
	}

	Array2d operator[](int i) const { return i == 0 ? pMin : pMax; }

	Bounds2d operator&(Bounds2d const &o) const {
		Bounds2d ret;
		ret.pMin = pMin.max(o.pMin);
		ret.pMax = pMax.min(o.pMax);
		return ret;
	}

	bool isDegenerate() const { return pMin.x() > pMax.x() || pMin.y() > pMax.y(); }
	bool isInside(Array2d const &o) const { return (o.x() >= pMin.x() && o.x() <= pMax.x() && o.y() >= pMin.y() && o.y() <= pMax.y()); }

	Array2d diagonal() const { return pMax - pMin; }
	Array2d corner(int i) const { return Array2d((*this)[(i & 1)][0], (*this)[(i & 2) ? 1 : 0][1]); }
};

struct DividedPatches {
	Bounds2d uvB1;
	Bounds2d uvB2;
	double tLower;

	DividedPatches(Bounds2d const &uvB1, Bounds2d const &uvB2, double tLower = std::numeric_limits<double>::infinity()) : uvB1 { uvB1 }, uvB2 { uvB2 }, tLower { tLower } { }
	bool operator<(DividedPatches const &o) const { return tLower > o.tLower; }
};

struct DividedPatch {
	Bounds2d uvB;
	double tLower;

	DividedPatch(Bounds2d const &uvB, double tLower = std::numeric_limits<double>::infinity()) : uvB { uvB }, tLower { tLower } { }
	bool operator<(DividedPatch const &o) const { return tLower > o.tLower; }
};

// Bicubic Bezier Functions
template <typename P>
static P blossomCubicBezier(std::span<P const> p, double u0, double u1, double u2) {
	P a[3] = { lerp(u0, p[0], p[1]), lerp(u0, p[1], p[2]), lerp(u0, p[2], p[3]) };
	P b[2] = { lerp(u1, a[0], a[1]), lerp(u1, a[1], a[2]) };
	return lerp(u2, b[0], b[1]);
}

template <typename P>
static P blossomBicubicBezier(std::span<P const> cp, Array2d const &uv0, Array2d const &uv1, Array2d const &uv2) {
	std::array<P, 4> q;
	for (int i = 0; i < 4; i++) {
		q[i] = blossomCubicBezier(cp.subspan(i * 4, 4), uv0.y(), uv1.y(), uv2.y());
	}
	return blossomCubicBezier<P>(q, uv0.x(), uv1.x(), uv2.x());
}

static Vector3d evaluateBicubicBezier(std::span<Vector3d const> cp, Array2d const &uv, Vector3d *dpdu = nullptr, Vector3d *dpdv = nullptr) {
	if (dpdu && dpdv) {
		std::array<Vector3d, 4> cpU, cpV;
		for (int i = 0; i < 4; i++) {
			cpU[i] = blossomCubicBezier(cp.subspan(i * 4, 4), uv.y(), uv.y(), uv.y());
			std::array<Vector3d, 3> a = { lerp(uv.x(), cp[i], cp[i + 4]),
										 lerp(uv.x(), cp[i + 4], cp[i + 8]),
										 lerp(uv.x(), cp[i + 8], cp[i + 12]) };
			std::array<Vector3d, 2> b = { lerp(uv.x(), a[0], a[1]),
										 lerp(uv.x(), a[1], a[2]) };
			cpV[i] = lerp(uv.x(), b[0], b[1]);
		}
		std::array<Vector3d, 3> cpU1 = { lerp(uv.x(), cpU[0], cpU[1]),
										lerp(uv.x(), cpU[1], cpU[2]),
										lerp(uv.x(), cpU[2], cpU[3]) };
		std::array<Vector3d, 2> cpU2 = { lerp(uv.x(), cpU1[0], cpU1[1]),
										lerp(uv.x(), cpU1[1], cpU1[2]) };
		if ((cpU2[1] - cpU2[0]).squaredNorm() > 0)
			*dpdu = 3 * (cpU2[1] - cpU2[0]);
		else {
			*dpdu = cpU[3] - cpU[0];
		}
		std::array<Vector3d, 3> cpV1 = { lerp(uv.y(), cpV[0], cpV[1]),
										lerp(uv.y(), cpV[1], cpV[2]),
										lerp(uv.y(), cpV[2], cpV[3]) };
		std::array<Vector3d, 2> cpV2 = { lerp(uv.y(), cpV1[0], cpV1[1]),
										lerp(uv.y(), cpV1[1], cpV1[2]) };
		if ((cpV2[1] - cpV2[0]).squaredNorm() > 0)
			*dpdv = 3 * (cpV2[1] - cpV2[0]);
		else {
			*dpdv = cpV[3] - cpV[0];
		}
		return blossomCubicBezier<Vector3d>(cpU, uv.x(), uv.x(), uv.x());
	}
	else {
		return blossomBicubicBezier(cp, uv, uv, uv);
	}
}

// Patch Functions
template <typename P>
static std::array<P, 16> divideBezierPatch(std::span<P const> cp, Bounds2d const &uvB) {
	std::array<P, 16> divCp;
	divCp[0] = blossomBicubicBezier(cp, uvB.corner(0), uvB.corner(0), uvB.corner(0));
	divCp[1] = blossomBicubicBezier(cp, uvB.corner(0), uvB.corner(0), uvB.corner(2));
	divCp[2] = blossomBicubicBezier(cp, uvB.corner(0), uvB.corner(2), uvB.corner(2));
	divCp[3] = blossomBicubicBezier(cp, uvB.corner(2), uvB.corner(2), uvB.corner(2));
	divCp[4] = blossomBicubicBezier(cp, uvB.corner(0), uvB.corner(0), uvB.corner(1));
	divCp[5] = blossomBicubicBezier(cp, uvB.corner(0), uvB.corner(0), uvB.corner(3));
	divCp[6] = blossomBicubicBezier(cp, uvB.corner(0), uvB.corner(2), uvB.corner(3));
	divCp[7] = blossomBicubicBezier(cp, uvB.corner(2), uvB.corner(2), uvB.corner(3));
	divCp[8] = blossomBicubicBezier(cp, uvB.corner(0), uvB.corner(1), uvB.corner(1));
	divCp[9] = blossomBicubicBezier(cp, uvB.corner(0), uvB.corner(1), uvB.corner(3));
	divCp[10] = blossomBicubicBezier(cp, uvB.corner(0), uvB.corner(3), uvB.corner(3));
	divCp[11] = blossomBicubicBezier(cp, uvB.corner(2), uvB.corner(3), uvB.corner(3));
	divCp[12] = blossomBicubicBezier(cp, uvB.corner(1), uvB.corner(1), uvB.corner(1));
	divCp[13] = blossomBicubicBezier(cp, uvB.corner(1), uvB.corner(1), uvB.corner(3));
	divCp[14] = blossomBicubicBezier(cp, uvB.corner(1), uvB.corner(3), uvB.corner(3));
	divCp[15] = blossomBicubicBezier(cp, uvB.corner(3), uvB.corner(3), uvB.corner(3));
	return divCp;
}



static double point2Box(const MatrixXd &x, const MatrixXd &u, const double tMax)
{
	int groups = int(x.rows()), items = int(x.cols());
	std::vector<Array2d> invlLists[6], finalList;
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

template <typename Func>
static double greedyIntersect(Func &&checkIntersect, Array2d &intsctUV)
{
	// Initialize the heap
	std::priority_queue<DividedPatch> heap;
	{
		Bounds2d uvB(Array2d(0, 0), Array2d(1, 1));
		const double tTent = checkIntersect(uvB);
		if (tTent < DeltaT)
			heap.emplace(uvB, tTent);
	}

	// Greedily search intersection points
	while (!heap.empty()) {
		const auto cur = heap.top();
		heap.pop();

		// Set uv of the middle point
		Array2d uvMid = (cur.uvB.pMin + cur.uvB.pMax) / 2;

		// Decide whether the algorithm converges
		if (cur.uvB.diagonal().maxCoeff() < MinDeltaUV) {
			intsctUV = uvMid;
			return cur.tLower;
		}

		// Divide the current patch into four pieces
		for (int i = 0; i < 4; ++i) {
			Bounds2d divUvB(cur.uvB.corner(i), uvMid);

			const double tTent = checkIntersect(divUvB);
			if (tTent < std::numeric_limits<double>::infinity())
				heap.emplace(divUvB, tTent);
		}
	}

	return -1;
}

static double point2patch(const Vector3d& p1, const Vector3d& v1, Array2d &intsctUV) {
	MatrixXd x(6, 16), u(6, 16);
	double x0[6], u0[6];
	for (int dim = 0; dim < 3; dim++)
		x0[dim] = p1[dim], x0[3 + dim] = -p1[dim],
		u0[dim] = v1[dim], u0[3 + dim] = -v1[dim];

	auto checkIntersection = [&](const Bounds2d &b) {
		const auto cpPosB = divideBezierPatch<Vector3d>(CpPos2, b);
		const auto cpVelB = divideBezierPatch<Vector3d>(CpVel2, b);
		// Get the coefficients of inequalities.
		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 4; j++)
				for (int dim = 0; dim < 3; dim++)
					x(dim, i << 2 | j) = cpPosB[i << 2 | j][dim] - x0[dim],
					x(3 + dim, i << 2 | j) = -cpPosB[i << 2 | j][dim] - x0[3 + dim],
					u(dim, i << 2 | j) = cpVelB[i << 2 | j][dim] - u0[dim],
					u(3 + dim, i << 2 | j) = -cpVelB[i << 2 | j][dim] - u0[3 + dim];

		// Perform ray bounding box intersection.
		return point2Box(x, u, DeltaT);
	};

	return greedyIntersect(std::move(checkIntersection), intsctUV);
}

static double ccdSample() {
	using steady_clock = std::chrono::steady_clock;
	using duration = std::chrono::duration<double>;
	const auto initialTime = steady_clock::now();

	constexpr double eps = .01;
	double minT = std::numeric_limits<double>::infinity();
	Array2d ap1, ap2, temptUV;

	for (double x1 = 0; x1 <= 1; x1 += eps) {
		for (double y1 = 0; y1 <= 1; y1 += eps) {
			Vector3d const p1 = evaluateBicubicBezier(CpPos1, { x1, y1 });
			Vector3d const v1 = evaluateBicubicBezier(CpVel1, { x1, y1 });
			
			double temptT = point2patch(p1, v1, temptUV);
			if(temptT >= 0 && temptT <= minT) minT=temptT, ap1={ x1, y1 }, ap2=temptUV;
		}
	}
	if(minT<=DeltaT&&minT>=0){
		std::cout << ap1.transpose() << std::endl << ap2.transpose() << std::endl;
		const auto endTime = steady_clock::now();
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





struct Line
{
	double k, b;
	Line(const double& k,const double& b): k(k), b(b) {}
	bool operator<(const Line &l) const { 
		return k < l.k || (k == l.k && b > l.b); // 相同斜率的直线中只有截距最大的被留下来
	}
	bool operator==(const Line &l) const {return k == l.k;}

	double lineIntersect_x(const Line &l) const {
		if(k==l.k){
			std::cout<<"parallel lines do not intersect at a single point!\n";
			exit(-1);
		}
		return -(b-l.b)/(k-l.k);
	}

};

static void getCH(std::vector<Line>& lines, std::vector<Line>& ch, std::vector<double>& pts) {
	lines.erase(std::unique(lines.begin(), lines.end()), lines.end()); // 去重
	ch.clear();
	pts.clear();
	pts.push_back(0);
	ch.push_back(lines[0]);
	if(DEBUG) for(const auto&l:lines)std::cout<<"lines:\t"<<l.k<<" "<<l.b<<"\n";
	double id = 1, intsctX = 0;
	while(id < lines.size()){
		// std::cout<<id<<"  "<<pts.size()<<"\n";
		while(!ch.empty()){
			intsctX = lines[id].lineIntersect_x(ch.back());
			if(intsctX<=pts.back()){
			// if(ch.back().k*pts.back()+ch.back().b<=lines[id].k*pts.back()+lines[id].b){
				pts.pop_back();
				ch.pop_back();
			}
			else break;
		}
		ch.push_back(lines[id]);
		pts.push_back(std::max(0.,intsctX));
		id++;
	}
	while(pts.back()>=DeltaT){
		pts.pop_back();
		ch.pop_back();
	}
	pts.push_back(DeltaT);
	if(DEBUG) for(const auto&l:ch)std::cout<<"ch:\t"<<l.k<<" "<<l.b<<"\n";
	if(DEBUG) for(const auto&pt:pts)std::cout<<"pt:\t"<<pt<<"\n";
	if(ch.empty()){
		std::cout<<"empty CH!\n";
		exit(-1);
	}
	if(ch.size()+1!=pts.size()){
		std::cout<<"segments and inflections are not compatible!\n";
		exit(-1);
	}
}

static Array2d linearCHIntersect(const std::vector<Line>& ch1, const std::vector<Line>& ch2, 
								const std::vector<double>& pts1, const std::vector<double>& pts2) {
	int id1=1, id2=1;
	double intvL=-1, intvR=-1;
	double sweep=0, lastsweep=0;
	bool stopInAdv = false; // 还没写上，但是比如：intvL==-1&&ch1[id1-1].k>ch1[id1-1].k

	auto checkSweepLine = [&] (const int id1, const int id2) {
		double y1=ch1[id1-1].k*sweep+ch1[id1-1].b;
		double y2=ch2[id2-1].k*sweep+ch2[id2-1].b;
		if (y1>y2){
			if(intvL!=-1)
				intvR = ch1[id1-1].lineIntersect_x(ch2[id2-1]);
		} // 如果k1==k2，那么必然前一个节点必然已经满足y1<y2了
		else if(y1<y2){
			if(intvL==-1)
				intvL = ch1[id1-1].lineIntersect_x(ch2[id2-1]); 
		}// 如果k1==k2，那么必然前一个节点必然已经满足y1<y2了
		// else{
		// 	//y1==y2
		// 	// if(intvL==-1) intvR = intvL = sweep;//这个不对，如果intv一直持续到deltaT就会变成只有一个点
		// 	if(intvL!=-1) 
		// 		intvR = sweep;
		// }
		if (DEBUG) std::cout<<id1<<"  "<<id2<<" /  "<<y1<<"  "<<y2<<" /  "<<intvL<<" "<<intvR<<"\n";
		lastsweep = sweep;
	};

	if(ch1[0].b<ch2[0].b)
		intvL = 0;
	// else if(ch1[0].b==ch2[0].b)
	// 	intvL = intvR = 0;
	while(id1<pts1.size()&&id2<pts1.size()){
		sweep = std::min(pts1[id1], pts2[id2]);
		if(sweep!=lastsweep)//并不知道这样跳过能不能更快
			checkSweepLine(id1, id2);
		if(pts1[id1] < pts2[id2]) id1++;
		else id2++;
	}
	while(id1<pts1.size()){
		sweep = pts1[id1];
		if(sweep!=lastsweep)//并不知道这样跳过能不能更快
			checkSweepLine(id1, id2);
		id1++;
	}
	while(id2<pts2.size()){
		sweep = pts2[id2];
		if(sweep!=lastsweep)//并不知道这样跳过能不能更快axesOBB
			checkSweepLine(id1, id2);
		id2++;
	}
	if(intvL!=-1 && intvR==-1)intvR=DeltaT;
	return Array2d(intvL,intvR);
}

// maybe no need to contruct from scratch?
static double PrimitiveCheck(Bounds2d const &divUvB1, Bounds2d const &divUvB2, const BB bbtype){
	auto const ptPos1 = divideBezierPatch<Vector3d>(CpPos1, divUvB1);
	auto const ptVel1 = divideBezierPatch<Vector3d>(CpVel1, divUvB1);
	auto const ptPos2 = divideBezierPatch<Vector3d>(CpPos2, divUvB2);
	auto const ptVel2 = divideBezierPatch<Vector3d>(CpVel2, divUvB2);

	auto setAxes = [&] (std::vector<Vector3d>& axes) {
		if(bbtype==BB::AABB){
			axes = {Vector3d::Unit(0), Vector3d::Unit(1), Vector3d::Unit(2)};
		}
		else if(bbtype==BB::OBB) {
			Vector3d 
			lu1 = (ptPos1[3]-ptPos1[0]+ptPos1[15]-ptPos1[12]),
			lv1 = (ptPos1[12]-ptPos1[0]+ptPos1[15]-ptPos1[3]),
			ln1 = lu1.cross(lv1);
			Vector3d 
			lu2 = (ptPos2[3]-ptPos2[0]+ptPos2[15]-ptPos2[12]),
			lv2 = (ptPos2[12]-ptPos2[0]+ptPos2[15]-ptPos2[3]),
			ln2 = lu2.cross(lv2);
			axes = {lu1,lv1,ln1,lu2,lv2,ln2, 
				lu1.cross(lu2), lu1.cross(lv2), lu1.cross(ln2), 
				lv1.cross(lu2), lv1.cross(lv2), lv1.cross(ln2), 
				ln1.cross(lu2), ln1.cross(lv2), ln1.cross(ln2)};
		}
	};

	// std::cout<<"done!\n";
	std::vector<Vector3d> axes;
	axes.clear();
	setAxes(axes);
    std::vector<Array2d> feasibleIntvs;
	feasibleIntvs.clear();

	auto AxisCheck=[&](const std::array<Vector3d, 16>& p1, const std::array<Vector3d, 16>& v1, 
						const std::array<Vector3d, 16>& p2, const std::array<Vector3d, 16>& v2, 
						const Vector3d& axis){
		std::vector<Line> lines1, lines2;
        std::vector<Line> ch1, ch2;
        std::vector<double> pts1, pts2;
		lines1.clear(); lines2.clear();
		ch1.clear(); ch2.clear();
		pts1.clear(); pts2.clear();
		for(int i = 0; i < 16; i++) lines1.emplace_back(v1[i].dot(axis), p1[i].dot(axis));
		for(int i = 0; i < 16; i++) lines2.emplace_back(-v2[i].dot(axis), -p2[i].dot(axis));
        std::sort(lines1.begin(), lines1.end());
        std::sort(lines2.begin(), lines2.end());
        getCH(lines1, ch1, pts2);
        getCH(lines2, ch2, pts2);
		// std::cout<<"getCHOK!\n";
        for(auto & l:ch2)
            l.k = -l.k, l.b = -l.b;
        const auto intvT = linearCHIntersect(ch1, ch2, pts1, pts2);
		if(SHOWANS) std::cout<<intvT.transpose()<<"\n";
		if(DEBUG) std::cin.get();
        if(intvT[0]!=-1)feasibleIntvs.push_back(intvT);
	};

    for(const auto& axis:axes)
        AxisCheck(ptPos1, ptVel1, ptPos2, ptVel2, axis);
    for(const auto& axis:axes)
        AxisCheck(ptPos2, ptVel2, ptPos1, ptVel1, axis);

	if(SHOWANS) std::cout<<"done!\n";

	if (feasibleIntvs.size()==0) return 0; //这意味着整段时间都有碰撞

	//无碰撞发生的并，剩下的就是有碰撞发生的
	std::sort(feasibleIntvs.begin(), feasibleIntvs.end(), 
		[](const Array2d& intv1, const Array2d& intv2){
			return (intv1(0)<intv2(0));
		});
	if(feasibleIntvs.size()==0 || (feasibleIntvs[0](0)>0)) return 0;
	double minT = feasibleIntvs[0](1);
	for(int i=1;i<feasibleIntvs.size();i++)
		if(feasibleIntvs[i](0)<=minT)
			minT=std::max(minT, feasibleIntvs[i](1));
		else break;
	if(minT<DeltaT)return minT;
	else return -1;

	// double minT = std::numeric_limits<double>::infinity(), maxT = -std::numeric_limits<double>::infinity();
	// for(const auto intv: feasibleIntvs){
	// 	minT = std::max(minT, intv(0));
	// 	maxT = std::min(maxT, intv(0));
	// }
	// // std::cout<<"done!\n";
	// if(minT<maxT)return minT;
	// else return -1;

	// Vector3d min1 = ptPos1[0];
	// Vector3d max1 = ptPos1[0];
	// Vector3d min2 = ptPos2[0];
	// Vector3d max2 = ptPos2[0];
	// for (int i = 1; i < 16; i++) {
	// 	min1 = min1.cwiseMin(ptPos1[i]);
	// 	max1 = max1.cwiseMax(ptPos1[i]);
	// 	min2 = min2.cwiseMin(ptPos2[i]);
	// 	max2 = max2.cwiseMax(ptPos2[i]);
	// }
	// for(int i = 0; i < 3; i++)
	// 	if (min1[i]>max2[i]||min2[i]>max1[i]) return false;
	// return true;
}
Array2d uv1, uv2;
static double ccd(const BB bbtype) {
	using steady_clock = std::chrono::steady_clock;
	using duration = std::chrono::duration<double>;
	const auto initialTime = steady_clock::now();

	std::priority_queue<DividedPatches> heap;
	Bounds2d initUv(Array2d(0, 0), Array2d(1, 1));
    double colTime = PrimitiveCheck(initUv, initUv, bbtype);
	if (colTime>=0 && colTime<=DeltaT)
		heap.emplace(initUv, initUv, colTime);

	while (!heap.empty()) {
		auto const cur = heap.top();
		heap.pop();
		cnt++;
		if(DEBUG) std::cout<<cnt<<"\n";
		// Set uv of the middle point
		Array2d uvMid1 = (cur.uvB1.pMin + cur.uvB1.pMax) / 2;
		Array2d uvMid2 = (cur.uvB2.pMin + cur.uvB2.pMax) / 2;
		// Decide whether the algorithm converges
		if (std::max(cur.uvB1.diagonal().maxCoeff(), cur.uvB2.diagonal().maxCoeff()) < MinDeltaUV) {
			// DEBUG=1;
			// BBcheck(cur.uvB1, cur.uvB2);
			std::cout << ((cur.uvB1.pMin + cur.uvB1.pMax) / 2).transpose() << std::endl;
			std::cout << ((cur.uvB2.pMin + cur.uvB2.pMax) / 2).transpose() << std::endl;
			uv1 = (cur.uvB1.pMin + cur.uvB1.pMax) / 2;
			uv2 = (cur.uvB2.pMin + cur.uvB2.pMax) / 2;
			const auto endTime = steady_clock::now();
			std::cout << "min time: "<<  cur.tLower << "used seconds: " <<
				duration(endTime - initialTime).count()
				<< std::endl;
			return cur.tLower;
		}

		// Divide the current patch into four-to-four pieces
		for (int i = 0; i < 4; i++) {
			Bounds2d divUvB1(cur.uvB1.corner(i), uvMid1);
			for (int j = 0; j < 4; j++) {
				Bounds2d divUvB2(cur.uvB2.corner(j), uvMid2);
				colTime = PrimitiveCheck(divUvB1, divUvB2, bbtype);//maybe also need timeLB?
				if (colTime>=0 && colTime<=DeltaT){
					heap.emplace(divUvB1, divUvB2, colTime);
                }
			}
		}
	}

	const auto endTime = steady_clock::now();
	std::cout << "used seconds: " <<
		duration(endTime - initialTime).count()
		<< std::endl;
	return -1;
}

void readinDoFs(){
	std::ifstream readin("DoFs.txt");
	for (int i = 0; i < 16; i++)
		for(int k=0;k<3;k++)
			readin>>CpPos1[i](k);
	for (int i = 0; i < 16; i++)
		for(int k=0;k<3;k++)
			readin>>CpPos2[i](k);
	for (int i = 0; i < 16; i++)
		for(int k=0;k<3;k++)
			readin>>CpVel1[i](k);
	for (int i = 0; i < 16; i++)
		for(int k=0;k<3;k++)
			readin>>CpVel2[i](k);
	readin.close();
}
void saveDoFs(){
	std::ofstream f("DoFs.txt");
	for(auto item : CpPos1)
		f<<item.transpose()<<"\n";
	for(auto item : CpPos2)
		f<<item.transpose()<<"\n";
	for(auto item : CpVel1)
		f<<item.transpose()<<"\n";
	for(auto item : CpVel2)
		f<<item.transpose()<<"\n";
	f.close();
}
void generatePatches(){
	// std::srand(0);
	for (int i = 0; i < 16; i++) {
		CpPos1[i] = Vector3d::Random() - Vector3d::Constant(.6);
		CpVel1[i] = Vector3d::Random()*0.3 + Vector3d::Constant(.6);
		CpPos2[i] = Vector3d::Random() + Vector3d::Constant(.6);
		CpVel2[i] = Vector3d::Random()*0.3 - Vector3d::Constant(.6);
	}
}
void randomTest(){
	using steady_clock = std::chrono::steady_clock;
	using duration = std::chrono::duration<double>;
	int cntAABB=0, cntSample=0;
	const int Kase = 50;
	double ans[2][Kase];

	const auto initOBB = steady_clock::now();
	std::srand(0);
	for(int kase=0;kase<Kase;kase++){
		generatePatches();
		ans[0][kase]=ccd(BB::OBB);
	}
	const auto endOBB = steady_clock::now();
	std::cout<<"OBB used seconds: "<<duration(endOBB - initOBB).count()*0.01<<"\n";

	const auto initAABB = steady_clock::now();
	std::srand(0);
	for(int kase=0;kase<Kase;kase++){
		generatePatches();
		ans[0][kase]=ccd(BB::AABB);
		// if(ccd(AABBcheck)!=-1)cntAABB++;
		// if(dcd(OBBcheck))cntOBB++;
		// if(cntAABB!=cntOBB){saveDoFs();exit(-1);}
	}
	const auto endAABB = steady_clock::now();
	std::cout<<"AABB used seconds: "<<duration(endAABB - initAABB).count()*0.01<<"\n";

	std::cout<<"OBB used seconds: "<<duration(endOBB - initOBB).count()*0.01<<"\n";
	std::cout<<"OBB used seconds: "<<duration(endOBB - initOBB).count()*0.01<<"\n";


	// const auto initSample = steady_clock::now();
	// std::srand(0);
	// for(int kase=0;kase<Kase;kase++){
	// 	generatePatches();
	// 	saveDoFs();
	// 	ans[1][kase]=ccdSample();
	// 	// if(ccdSample()!=-1)cntSample++;
	// }
	// const auto endSample = steady_clock::now();
	// std::cout<<"Sample used seconds: "<<duration(endSample - initSample).count()*0.01<<"\n";

	// Array2d checkUV;
	// Vector3d const p1 = evaluateBicubicBezier(CpPos1, uv1);
	// Vector3d const v1 = evaluateBicubicBezier(CpVel1, uv1);
	// Vector3d const p2 = evaluateBicubicBezier(CpPos2, uv2);
	// Vector3d const v2 = evaluateBicubicBezier(CpVel2, uv2);
	// std::cout<<"trajectories: "<<v1.transpose()<<"   "<<p1.transpose()<<"\n";
	// std::cout<<"trajectories: "<<v2.transpose()<<"   "<<p2.transpose()<<"\n";
	// std::cout<<"\ncheck AABB ans: min time "<<point2patch(p1, v1, checkUV);
	// std::cout<<"at\n"<<uv1.transpose()<<"\n"<<checkUV.transpose();

	// std::cout<<"\n\n";
	// int acc=0;
	// for(int k=0;k<Kase;k++)
	// 	if(std::abs(ans[0][k]-ans[1][k])<1e-3)acc++;
	// 	else {std::cout<<ans[0][k]<<"  "<<ans[1][k]<<"\n";}

	// std::cout<<"AABB used seconds: "<<duration(endAABB - initAABB).count()*0.01<<"\n";
	// std::cout<<"Sample used seconds: "<<duration(endSample - initSample).count()*0.01<<"\n";
	// std::cout<<"acc: "<< acc/Kase;
	// std::cout<<"LPOBB: "<<cntLPOBB<<"used seconds: "<<duration(endLPOBB - initLPOBB).count()*0.01<<"\n";
}

void singleTest(){

{	std::srand(0);
	generatePatches();
	double t = ccd(BB::AABB);
	Array2d checkUV;
	Vector3d const p1 = evaluateBicubicBezier(CpPos1, uv1);
	Vector3d const v1 = evaluateBicubicBezier(CpVel1, uv1);
	Vector3d const p2 = evaluateBicubicBezier(CpPos2, uv2);
	Vector3d const v2 = evaluateBicubicBezier(CpVel2, uv2);
	Vector3d const pt1=(v1*t+p1), pt2=(v2*t+p2);
	std::cout<<"trajectories: "<<v1.transpose()<<"   "<<p1.transpose()<<"\n";
	std::cout<<"trajectories: "<<v2.transpose()<<"   "<<p2.transpose()<<"\n";
	std::cout<<"pos at: "<<pt1.transpose()<<"     "<<pt2.transpose()<<"\n";
	std::cout<<"delta: "<<(pt2-pt1).norm()<<"\n";
	std::cout<<"check AABB ans: min time "<<point2patch(p1, v1, checkUV)<<"\n\n";}

	{std::srand(0);
	generatePatches();
	double t = ccd(BB::OBB);

	Array2d checkUV;
	Vector3d const p1 = evaluateBicubicBezier(CpPos1, uv1);
	Vector3d const v1 = evaluateBicubicBezier(CpVel1, uv1);
	Vector3d const p2 = evaluateBicubicBezier(CpPos2, uv2);
	Vector3d const v2 = evaluateBicubicBezier(CpVel2, uv2);
	Vector3d const pt1=(v1*t+p1), pt2=(v2*t+p2);
	std::cout<<"trajectories: "<<v1.transpose()<<"   "<<p1.transpose()<<"\n";
	std::cout<<"trajectories: "<<v2.transpose()<<"   "<<p2.transpose()<<"\n";
	std::cout<<"pos at: "<<pt1.transpose()<<"     "<<pt2.transpose()<<"\n";
	std::cout<<"delta: "<<(pt2-pt1).norm()<<"\n";
	std::cout<<"check OBB ans: min time "<<point2patch(p1, v1, checkUV);}
}

int main(){
	// readinDoFs();
	// detectIntersect();
    // generatePatches();
    // std::cout<<"AABB: "<<ccd(AABBcheck)<<"\n";
    // std::cout<<"sample: "<<ccdSample()<<"\n";
	// const std::vector<Line> ch1 = {Line(0.602529, -0.824276)};
	// const std::vector<Line> ch2 = {Line(-0.382798, -0.143568)}; 
	// const std::vector<double> pts1 = {0, 1}, pts2=pts1; 
	// std::cout<<"final res: "<<(linearCHIntersect(ch1,ch2,pts1,pts2)).transpose();

	// Vector3d v1(0.587538, 0.647203, 0.591826),   p1(-0.431234, -0.555342, -0.578118),
	// v2(-0.80066, -0.359232, -0.719931),   p2( 0.226999, -0.0781256,  0.0438646);
	// double t=0.474162;
	// std::cout<<v1*t+p1<<"\n\n\n"<<v2*t+p2;

	// randomTest();
	singleTest();

}