#include <Eigen/Dense>

using Eigen::Array2d;
using Eigen::Vector3d;

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
static constexpr double MinDeltaUV = 1e-4;
static constexpr double DeltaT = 1e-4;

std::array<Vector3d, 16> CpPos1 = {
	Vector3d(0, 0, 0), Vector3d(1, 0, 0), Vector3d(2, 0, 0), Vector3d(3, 0, 0.1),
	Vector3d(0, 1, 0), Vector3d(1, 1, 0), Vector3d(2, 1, 0), Vector3d(3, 1, 0.1),
	Vector3d(0, 2, 0), Vector3d(1, 2, 0), Vector3d(2, 2, 0), Vector3d(3, 2, 0.1),
	Vector3d(0, 3, 0), Vector3d(1, 3, 0), Vector3d(2, 3, 0), Vector3d(3, 3, 0.1),
};
std::array<Vector3d, 16> CpPos2 = {
	Vector3d(0, 0, 0.1), Vector3d(1, 0, 0.1), Vector3d(2, 0, 0.1), Vector3d(3, 0, 0.2),
	Vector3d(0, 1, 0.1), Vector3d(1, 1, 0.1), Vector3d(2, 1, 0.1), Vector3d(3, 1, 0.2),
	Vector3d(0, 2, 0.1), Vector3d(1, 2, 0.1), Vector3d(2, 2, 0.1), Vector3d(3, 2, 0.2),
	Vector3d(0, 3, 0.1), Vector3d(1, 3, 0.1), Vector3d(2, 3, 0.1), Vector3d(3, 3, 0.2),
};
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

struct DividedPatch {
	Bounds2d uvB1;
	Bounds2d uvB2;
	double dLower;

	DividedPatch(Bounds2d const &uvB1, Bounds2d const &uvB2, double dLower = std::numeric_limits<double>::infinity()) : uvB1 { uvB1 }, uvB2 { uvB2 }, dLower { dLower } { }
	bool operator<(DividedPatch const &o) const { return dLower > o.dLower; }
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

// static double solveNaive() {
// 	using steady_clock = std::chrono::steady_clock;
// 	using duration = std::chrono::duration<double>;
// 	const auto initialTime = steady_clock::now();
// 	constexpr double eps = .01;
// 	double ans = std::numeric_limits<double>::infinity();
// 	Array2d ap1, ap2;
// 	for (double x1 = 0; x1 <= 1; x1 += eps) {
// 		for (double y1 = 0; y1 <= 1; y1 += eps) {
// 			Vector3d const p1 = evaluateBicubicBezier(CpPos1, { x1, y1 });
// 			for (double x2 = 0; x2 <= 1; x2 += eps) {
// 				for (double y2 = 0; y2 <= 1; y2 += eps) {
// 					Vector3d const p2 = evaluateBicubicBezier(CpPos2, { x2, y2 });
// 					double const dist2 = (p1 - p2).squaredNorm();
// 					if (dist2 < ans) {
// 						ans = dist2;
// 						ap1 = { x1, y1 };
// 						ap2 = { x2, y2 };
// 					}
// 				}
// 			}
// 		}
// 	}
// 	std::cout << ap1.transpose() << std::endl << ap2.transpose() << std::endl;
// 				const auto endTime = steady_clock::now();
// 			std::cout << "used seconds: " <<
// 				duration(endTime - initialTime).count()
// 				<< std::endl;
// 	return ans;
// }

// static bool LPOBBcheck(Bounds2d const &divUvB1, Bounds2d const &divUvB2){
// 	auto const pt1 = divideBezierPatch<Vector3d>(CpPos1, divUvB1);
// 	auto const pt2 = divideBezierPatch<Vector3d>(CpPos2, divUvB2);
// 	Vector3d 
// 	lu1 = (pt1[3]-pt1[0]+pt1[15]-pt1[12]),
// 	lv1 = (pt1[12]-pt1[0]+pt1[15]-pt1[3]),
// 	ln1 = lu1.cross(lv1);
// 	Vector3d 
// 	lu2 = (pt2[3]-pt2[0]+pt2[15]-pt2[12]),
// 	lv2 = (pt2[12]-pt2[0]+pt2[15]-pt2[3]),
// 	ln2 = lu2.cross(lv2);
// 	std::array<Vector3d,15> axes = {lu1,lv1,ln1,lu2,lv2,ln2, 
// 		lu1.cross(lu2), lu1.cross(lv2), lu1.cross(ln2), 
// 		lv1.cross(lu2), lv1.cross(lv2), lv1.cross(ln2), 
// 		ln1.cross(lu2), ln1.cross(lv2), ln1.cross(ln2)};
// 	for (auto& axis:axes){
// 		if(axis.norm()<1e-6)continue;
// 		else axis.normalize();
// 		double projMax1=-std::numeric_limits<double>::infinity(),projMin1=std::numeric_limits<double>::infinity();
// 		for(const auto &p:pt1){
// 			double proj=p.dot(axis);
// 			projMax1=std::max(proj, projMax1); projMin1=std::min(proj, projMin1);
// 		}
// 		double projMax2=-std::numeric_limits<double>::infinity(),projMin2=std::numeric_limits<double>::infinity();
// 		for(const auto &p:pt2){
// 			double proj=p.dot(axis);
// 			projMax2=std::max(proj, projMax2); projMin2=std::min(proj, projMin2);
// 		}
// 		if (projMin2>projMax1||projMin1>projMax2){
// 			// std::cout<<"no intersect: project "<<projMax1-projMin1<<"  "<<projMax2-projMin2<<" on axis "<< axis.transpose()<<"\n";
// 			return false;
// 		}
// 	}
// 	// std::cout<<"does intersect!\n";
// 	return true;
// }

// static void intervalIntersection(const std::vector<Vector2d> &list1, const std::vector<Vector2d> &list2, std::vector<Vector2d> &res)
// {
// 	res.clear();
// 	for (int idx1 = 0, idx2 = 0; idx1 < list1.size() && idx2 < list2.size();) {
// 		double start1 = list1[idx1][0], end1 = list1[idx1][1];
// 		double start2 = list2[idx2][0], end2 = list2[idx2][1];
// 		if (end1 >= start2 && end2 >= start1)
// 			res.push_back(Vector2d(std::max(start1, start2), std::min(end1, end2)));
// 		if (end1 < end2) idx1++;
// 		else idx2++;
// 	}
// }

// static double point2Box(const MatrixXd &x, const MatrixXd &u, const double tMax)
// {
// 	int groups = int(x.rows()), items = int(x.cols());
// 	std::vector<Vector2d> invlLists[6], finalList;
// 	for (int num = 0; num < groups; num++) {
// 		double left = 0. - Epsilon, right = tMax + Epsilon;
// 		for (int i = 0; i < items; i++) {
// 			if (u(num, i) > Epsilon) left = std::max(left, -x(num, i) / u(num, i));
// 			else if (-u(num, i) > Epsilon) right = std::min(right, -x(num, i) / u(num, i));
// 			else if (x(num, i) < Epsilon) { left = tMax; break; }
// 		}
// 		invlLists[num].clear();
// 		if (left >= right || left >= tMax || right <= 0.) invlLists[num].push_back(Vector2d(0., tMax));
// 		else {
// 			if (left >= 0.) invlLists[num].push_back(Vector2d(0., left));
// 			if (right <= tMax) invlLists[num].push_back(Vector2d(right, tMax));
// 		}
// 	}
// 	finalList.clear();
// 	for (int num = 1; num < groups; num++) {
// 		intervalIntersection(invlLists[num - 1], invlLists[num], finalList);
// 		invlLists[num].clear();
// 		for (auto &item : finalList) invlLists[num].push_back(item);
// 	}
// 	if (finalList.empty())
// 		return std::numeric_limits<double>::infinity();
// 	return finalList[0][0];
// }

// maybe no need to contruct from scratch?
static double AABBcheck(Bounds2d const &divUvB1, Bounds2d const &divUvB2){
	auto const ptPos1 = divideBezierPatch<Vector3d>(CpPos1, divUvB1);
	auto const ptVel1 = divideBezierPatch<Vector3d>(CpVel1, divUvB1);
	auto const ptPos2 = divideBezierPatch<Vector3d>(CpPos2, divUvB2);
	auto const ptVel2 = divideBezierPatch<Vector3d>(CpVel2, divUvB2);

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
            // double x=-(b-l.b)/(k-l.k), y=k*x+b;
            // return Array2s(x,y);
        }
    };

    auto getCH = [&] (std::vector<Line>& lines) {
        lines.erase(std::unique(lines.begin(), lines.end()), lines.end()); // 去重
        if(lines.size() < 2) return lines;

        std::vector<std::pair<Line, double>> res; // (直线，开始点)
        res.push_back(std::make_pair(lines[0], 0));
        double id = 1, intsctX = 0;
        while(id < lines.size()){
            while(!res.empty()){
                intsctX = -(lines[id].b-res.back().first.b)/(lines[id].k-res.back().first.k);
                if(intsctX<=res.back().second)res.pop_back();
            }
            res.push_back(std::make_pair(lines[id],intsctX));
            id++;
        }

        return res;
    };

    auto linearCHIntersect = [&] (const std::vector<Line>& ch1, const std::vector<Line>& ch2, 
                                    const std::vector<double>& pts1, const std::vector<double>& pts2) {
        if(ch1.empty()||ch1.empty()){
            std::cout<<"empty CH!\n";
            exit(-1);
        }
        if(ch1.size()!=pts1.size()||ch2.size()!=pts2.size()){
            std::cout<<"segments and inflections are not compatible!\n";
            exit(-1);
        }

        int id1=0, id2=0;
        double intvL=-1, intvR=-1;
        double sweep=-1, lastsweep=-1;

        auto checkSweepLine = [&] (const int id1, const int id2) {
            if(sweep==lastsweep)continue;//并不知道这样跳过能不能更快
            double y1=ch1[id1].k*sweep+ch1[id1].b;
            double y2=l2.k*sweep+l2.b;
            if(y1==y2){
                if(intvL==-1)intvL=sweep;
                else if()
            }
            if(y1<=y2&&intvL==-1){
                if(l1.k!=l2.k)
                    intvL = l1.lineIntsct_x(l2.first);
                else intvL = sweep;
            }
            else if (y1>y2&&intvL!=-1)
                intvR = l1.lineIntsct_x(l2.first);
                if(l1.k==l2.k){std::cout<<"how come!\n";exit(-1);}
                intvR = l1.lineIntsct_x(l2.first);
                else intvR = lastsweep;
            lastsweep = sweep;
        };

        while(id1<lines1.size()&&id2<lines2.size()){
            sweep = std::min(lines1[id1].second, lines2[id2].second);
            checkSweepLine();
            if(lines1[id1].second < lines2[id2].second) id1++;
        }
        while(id1<lines1.size()){
            sweep = lines1[id1].second;
            checkSweepLine(id1, id2-1);
            id1++;
        }
        while(id2<lines2.size()){
            sweep = lines2[id2].second;
            checkSweepLine(id1-1, id2);
            id2++;
        }
        if(intvL!=-1 && intvR==-1)intvR=DeltaT;
        return Array(intvL,intvR);


        // std::vector<double> sweepline;
        // for(const auto& l: lines1){
        //     if(l.second>=DeltaT) break;
        //     sweepline.push_back(l.second);
        // }
        // for(const auto& l: lines2){
        //     if(l.second>=DeltaT) break;
        //     sweepline.push_back(l.second);
        // }
        // std::sort(sweepline.begin(), sweepline.end());
        // sweepline.push_back(DeltaT);
        // sweepline.erase(std::unique(sweepline.begin(), sweepline.end()), sweepline.end()); // 去重
        // int id1 = 0, id2 = 0;
        // double intvL = -1, intvR = -1;
        // bool flag = false;
        // for(const double& x: sweepline){
        //     double nextX1 = id1<lines1.size()-1 ? lines1[id1+1].second : DeltaT;
        //     double nextX1 = id2<lines2.size()-1 ? lines2[id2+1].second : DeltaT;
        //     if(x > nextX1) id1++;
        //     else if(x > nextX2) id2++;
        //     double y1=lines1[id1].first.k*x+lines1[id1].first.b;
        //     double y2=lines2[id2].first.k*x+lines2[id2].first.b;
        //     if(y1<=y2&&!flag){
        //         // start intsct
        //         intvL=std::max(0, lineIntsct(lines1[id1].first,lines2[id2].first));
        //         flag=true;
        //     }
        //     else if (y1>=y2&&flag){
        //         // end intsct
        //         intvR=std::min(DeltaT, lineIntsct(lines1[id1].first,lines2[id2].first));
        //         flag=false;
        //         break;
        //     }
        // }//平行，只相交一个点????
        // if(flag)intvR=DeltaT;
        // return Array(intvL,intvR);
    }

    // std::array<Array2d, 6> feasibleInts = ;
    for(int dim = 0; dim < 3; dim++){
        std::vector<Line> lines1, lines2;
        std::vector<Line> ch1, ch2;
        std::vector<double> pts1, pts2;
        for(int i = 0; i < 16; i++) lines1.emplace_back(ptVel1(i), ptPos1(i));
        for(int i = 0; i < 16; i++) lines2.emplace_back(-ptVel2(i), -ptPos2(i));
        std::sort(lines1.begin(), lines1.end());
        std::sort(lines2.begin(), lines2.end());
        getCH(lines1, ch1, pts2);
        getCH(lines2, ch2, pts2);
        for(auto & l:ch2)
            l.first.k = -l.first.k, l.first.b = -l.first.b;
        const auto intvT = linearCHIntersect(ch1, ch2, pts1, pts2);
        feasibleIntvs.push_?
    }


	// Vector3d min1 = pt1[0];
	// Vector3d max1 = pt1[0];
	// Vector3d min2 = pt2[0];
	// Vector3d max2 = pt2[0];
	// for (int i = 1; i < 16; i++) {
	// 	min1 = min1.cwiseMin(pt1[i]);
	// 	max1 = max1.cwiseMax(pt1[i]);
	// 	min2 = min2.cwiseMin(pt2[i]);
	// 	max2 = max2.cwiseMax(pt2[i]);
	// }
	// for(int i = 0; i < 3; i++)
	// 	if (min1[i]>max2[i]||min2[i]>max1[i]) return false;
	// return true;
    return 0;
}

static double ccd(const std::function<double(const Bounds2d &, const Bounds2d &)> &BBcheck) {
	using steady_clock = std::chrono::steady_clock;
	using duration = std::chrono::duration<double>;
	const auto initialTime = steady_clock::now();

	std::priority_queue<DividedPatch> heap;
	Bounds2d initUv(Array2d(0, 0), Array2d(1, 1));
    double colTime = BBcheck(initUv, initUv);
	if (colTime>=0 && colTime<=DeltaT)
		heap.emplace(initUv, initUv, colTime);

	while (!heap.empty()) {
		auto const cur = heap.top();
		heap.pop();
		cnt++;
		// Set uv of the middle point
		Array2d uvMid1 = (cur.uvB1.pMin + cur.uvB1.pMax) / 2;
		Array2d uvMid2 = (cur.uvB2.pMin + cur.uvB2.pMax) / 2;
		// Decide whether the algorithm converges
		if (std::max(cur.uvB1.diagonal().maxCoeff(), cur.uvB2.diagonal().maxCoeff()) < MinDeltaUV) {
			std::cout << ((cur.uvB1.pMin + cur.uvB1.pMax) / 2).transpose() << std::endl;
			std::cout << ((cur.uvB2.pMin + cur.uvB2.pMax) / 2).transpose() << std::endl;
			const auto endTime = steady_clock::now();
			std::cout << "used seconds: " <<
				duration(endTime - initialTime).count()
				<< std::endl;
			return cur.dLower;
		}

		// Divide the current patch into four-to-four pieces
		for (int i = 0; i < 4; i++) {
			Bounds2d divUvB1(cur.uvB1.corner(i), uvMid1);
			for (int j = 0; j < 4; j++) {
				Bounds2d divUvB2(cur.uvB2.corner(j), uvMid2);
				colTime = BBcheck(divUvB1, divUvB2);//maybe also need timeLB?
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
	return false;
}

int computeMinDist() {
	std::srand(std::time(nullptr));
	for (int i = 0; i < 16; i++) {
		CpPos1[i] += Vector3d::Random() * .6 - Vector3d::Constant(.3);
		CpPos2[i] += Vector3d::Random() * .6 - Vector3d::Constant(.3);
	}
	// std::cout << solveNaive() << std::endl;
	std::cout << solve() << std::endl << cnt << std::endl;
	return 0;
}

int detectIntersect(){
	cnt=0;
	std::cout << dcd(AABBcheck) << std::endl << cnt << std::endl;
	cnt=0;
	std::cout << dcd(OBBcheck) << std::endl << cnt << std::endl;
	return 0;
}
void readinDoFs(){
	std::ifstream readin("DoFs.txt");
	for (int i = 0; i < 16; i++)
		for(int k=0;k<3;k++)
			readin>>CpPos1[i](k);
	for (int i = 0; i < 16; i++)
		for(int k=0;k<3;k++)
			readin>>CpPos2[i](k);
	readin.close();
}
void saveDoFs(){
	std::ofstream f("DoFs.txt");
	for(auto item : CpPos1)
		f<<item.transpose()<<"\n";
	for(auto item : CpPos2)
		f<<item.transpose()<<"\n";
	f.close();
}
void generatePatches(){
	std::srand(std::time(nullptr));
	for (int i = 0; i < 16; i++) {
		CpPos1[i] = Vector3d::Random() - Vector3d::Constant(.5);
		CpVel1[i] = Vector3d::Random() - Vector3d::Constant(.3);
		CpPos2[i] = Vector3d::Random() + Vector3d::Constant(.5);
		CpVel2[i] = Vector3d::Random() + Vector3d::Constant(.3);
	}
}
void randomTest(){
	using steady_clock = std::chrono::steady_clock;
	using duration = std::chrono::duration<double>;
	int cntAABB=0, cntOBB=0, cntLPOBB=0;
	const auto initAABB = steady_clock::now();
	std::srand(0);
	for(int kase=0;kase<100;kase++){
		generatePatches();
		if(dcd(AABBcheck))cntAABB++;
		// if(dcd(OBBcheck))cntOBB++;
		// if(cntAABB!=cntOBB){saveDoFs();exit(-1);}
	}
	const auto endAABB = steady_clock::now();
	const auto initOBB = steady_clock::now();
	std::srand(0);
	for(int kase=0;kase<100;kase++){
		generatePatches();
		if(dcd(OBBcheck))cntOBB++;
	}
	const auto endOBB = steady_clock::now();
	const auto initLPOBB = steady_clock::now();
	std::srand(0);
	for(int kase=0;kase<100;kase++){
		generatePatches();
		if(dcd(LPOBBcheck))cntLPOBB++;
	}
	const auto endLPOBB = steady_clock::now();
	std::cout<<"AABB: "<<cntAABB<<"used seconds: "<<duration(endAABB - initAABB).count()*0.01<<"\n";
	std::cout<<"OBB: "<<cntOBB<<"used seconds: "<<duration(endOBB - initOBB).count()*0.01<<"\n";
	std::cout<<"LPOBB: "<<cntLPOBB<<"used seconds: "<<duration(endLPOBB - initLPOBB).count()*0.01<<"\n";
}
int main(){
	// readinDoFs();
	// detectIntersect();
    generatePatches();
    std::cout<<
}