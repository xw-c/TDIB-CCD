# pragma once
#include "mathOps.h"
#include "triBezier.h"

class Edge{
public:
	static const int cntCp = 2;
	Vector3d e1,e2;	
	Edge() {}
	Edge(const Vector3d&a, const Vector3d&b):e1(a), e2(b){}
	Vector3d lerp(const double& u)const {return u*e2+(1-u)*e1;}
	std::array<Vector3d, 2> dividePatch(const Array2d& coords) const {
		std::array<Vector3d, 2> pts={lerp(coords[0]), lerp(coords[1])};
		return pts;
	}
};
class Point{
public:
	static const int cntCp = 1;
	Vector3d p;	
	Point() {}
	Point(const Vector3d&a): p(a){}
	std::array<Vector3d, 1> dividePatch() const {
		std::array<Vector3d, 1> pts={p};
		return pts;
	}
};
class Face: public TriLinearBezier{
public:
	Face():TriLinearBezier(){}
	Face(const std::array<Vector3d, 3>& p): TriLinearBezier(p) {}
};

template <typename PrimType1, typename PrimType2>
static double primitiveCheck(const std::array<Vector3d,PrimType1::cntCp> &ptPos1, 
						const std::array<Vector3d,PrimType1::cntCp> &ptVel1, 
						const std::array<Vector3d,PrimType2::cntCp> &ptPos2, 
						const std::array<Vector3d,PrimType2::cntCp> &ptVel2,
						const double upperTime = DeltaT){
	// std::cout<<"done!\n";
	std::array<Vector3d,3> axes = {Vector3d::Unit(0), Vector3d::Unit(1), Vector3d::Unit(2)};
	// axes.clear();
	// if(bbType==BoundingBoxType::AABB){
	// 	axes = {Vector3d::Unit(0), Vector3d::Unit(1), Vector3d::Unit(2)};
	// }
	// else if(bbType==BoundingBoxType::DOP14){
	// 	axes = {Vector3d::Unit(0), Vector3d::Unit(1), Vector3d::Unit(2), 
	// 			Vector3d(1,1,1).normalized(), Vector3d(-1,1,1).normalized(), Vector3d(-1,-1,1).normalized()};
	// }
	// else if(bbType==BoundingBoxType::OBB){
	// 	Vector3d lu1 = (ptPos1[Edge::cornerId(2)]-ptPos1[Edge::cornerId(0)]).normalized();//u延展的方向
	// 	Vector3d lv1 = (ptPos1[Edge::cornerId(1)]-ptPos1[Edge::cornerId(0)]);//v延展的方向
	// 	lv1 = (lv1-lv1.dot(lu1)*lu1).eval();
	// 	Vector3d ln1 = lu1.cross(lv1);
	// 	Vector3d lu2 = (ptPos2[Edge::cornerId(2)]-ptPos2[Edge::cornerId(0)]).normalized();//u延展的方向
	// 	Vector3d lv2 = (ptPos2[Edge::cornerId(1)]-ptPos2[Edge::cornerId(0)]);//v延展的方向
	// 	lv2 = (lv2-lv2.dot(lu2)*lu2).eval();
	// 	Vector3d ln2 = lu2.cross(lv2);
	// 	axes = {lu1,lv1,ln1,lu2,lv2,ln2, 
	// 		lu1.cross(lu2), lu1.cross(lv2), lu1.cross(ln2), 
	// 		lv1.cross(lu2), lv1.cross(lv2), lv1.cross(ln2), 
	// 		ln1.cross(lu2), ln1.cross(lv2), ln1.cross(ln2)};
	// }
	std::vector<Array2d> feasibleIntvs;
	feasibleIntvs.clear();

	for(auto& axis:axes){
		// AxisCheck(ptPos1, ptVel1, ptPos2, ptVel2, axis);
		std::vector<Line> lines1, lines2;
		std::vector<Line> ch1, ch2;
		std::vector<double> pts1, pts2;
		lines1.clear(); lines2.clear();
		for(int i = 0; i < PrimType1::cntCp; i++) lines1.emplace_back(ptVel1[i].dot(axis), ptPos1[i].dot(axis));
		for(int i = 0; i < PrimType2::cntCp; i++) lines2.emplace_back(ptVel2[i].dot(axis), ptPos2[i].dot(axis));
		
		std::sort(lines1.begin(), lines1.end());
		std::sort(lines2.begin(), lines2.end()); 
		// auto CHCheck=[&](std::vector<Line> lineSet1, std::vector<Line> lineSet2)
		{
			ch1.clear(); ch2.clear();
			pts1.clear(); pts2.clear();
			getCH(lines1, ch1, pts1, true, upperTime);
			getCH(lines2, ch2, pts2, false, upperTime);//1/20
			const auto intvT = linearCHIntersect(ch1, ch2, pts1, pts2, upperTime);
			if(intvT[0]!=-1)feasibleIntvs.push_back(intvT);
		}
		{
			ch1.clear(); ch2.clear();
			pts1.clear(); pts2.clear();
			getCH(lines2, ch1, pts1, true, upperTime);
			getCH(lines1, ch2, pts2, false, upperTime);//1/20
			const auto intvT = linearCHIntersect(ch1, ch2, pts1, pts2, upperTime);
			if(intvT[0]!=-1)feasibleIntvs.push_back(intvT);
		}
	}

	if (feasibleIntvs.size()==0) return 0; //这意味着整段时间都有碰撞
	//无碰撞发生的并，剩下的就是有碰撞发生的
	std::sort(feasibleIntvs.begin(), feasibleIntvs.end(), 
		[](const Array2d& intv1, const Array2d& intv2){
			return (intv1(0)<intv2(0));
		});
	if(feasibleIntvs[0](0)>0) return 0;
	double minT = feasibleIntvs[0](1);
	for(int i=1;i<feasibleIntvs.size();i++)
		if(feasibleIntvs[i](0)<minT) //不能加等，因为无碰撞给的是开区间，如果有),(的情况加等号会把这个情况漏掉
			minT=std::max(minT, feasibleIntvs[i](1));
		else break;
	
	if(minT<upperTime)return minT;
	else return -1;
}

void getEECH(Line lines[2], std::vector<Line>& ch, std::vector<double>& pts,
			 const bool getUpperCH = true, const double& upperT = DeltaT) {
	if(!getUpperCH)std::swap(lines[0],lines[1]);
	ch.clear();
	pts.clear();

	if(lines[0]==lines[1]){
		pts.push_back(0);
		pts.push_back(upperT);
		ch.push_back(lines[0]);
	}
	else{
		double intsctX = lineIntersect_x(lines[1], lines[0]);
		if(intsctX<=0){
			pts.push_back(0);
			pts.push_back(upperT);
			ch.push_back(lines[1]);
		}
		else if(intsctX>=upperT){
			pts.push_back(0);
			pts.push_back(upperT);
			ch.push_back(lines[0]);
		}
		else{
			pts.push_back(0);
			pts.push_back(intsctX);
			pts.push_back(upperT);
			ch.push_back(lines[0]);
			ch.push_back(lines[1]);
		}
	}
}
static double primitiveEECheck(const Vector3d &ptPos10, const Vector3d &ptPos11,
						const Vector3d &ptVel10, const Vector3d &ptVel11,
						const Vector3d &ptPos20, const Vector3d &ptPos21,
						const Vector3d &ptVel20, const Vector3d &ptVel21,
						const double upperTime = DeltaT){
	// Vector3d axes[3] = {Vector3d::Unit(0), Vector3d::Unit(1), Vector3d::Unit(2)};
	Vector3d axes[3] = {ptPos11-ptPos10, ptPos21-ptPos20, (ptPos11-ptPos10).cross(ptPos21-ptPos20)};
	if(axes[2].norm()<1e-12){
		std::cerr<<"parallel!\n"<<axes[0]<<"\n\n"<<axes[1];
		exit(-1);
	}
	for(auto& axis:axes)axis.normalize();

	std::vector<Array2d> feasibleIntvs;
	feasibleIntvs.clear();
	for(auto& axis:axes){
		Line lines1[2]={Line(ptVel10.dot(axis), ptPos10.dot(axis)), 
						Line(ptVel11.dot(axis), ptPos11.dot(axis))}, 
				lines2[2] = {Line(ptVel20.dot(axis), ptPos20.dot(axis)), 
						Line(ptVel21.dot(axis), ptPos21.dot(axis))};
		std::vector<Line> ch1, ch2;
		std::vector<double> pts1, pts2;
		std::sort(lines1, lines1+2);
		std::sort(lines2, lines2+2); 
		{
			auto l1=lines1, l2=lines2;
			ch1.clear(); ch2.clear();
			pts1.clear(); pts2.clear();
			getEECH(l1, ch1, pts1, true, upperTime);
			getEECH(l2, ch2, pts2, false, upperTime);//1/20
			const auto intvT = linearCHIntersect(ch1, ch2, pts1, pts2, upperTime);
			if(intvT[0]!=-1)feasibleIntvs.push_back(intvT);
		}
		{
			auto l1=lines1, l2=lines2;
			ch1.clear(); ch2.clear();
			pts1.clear(); pts2.clear();
			getEECH(l2, ch1, pts1, true, upperTime);
			getEECH(l1, ch2, pts2, false, upperTime);//1/20
			const auto intvT = linearCHIntersect(ch1, ch2, pts1, pts2, upperTime);
			if(intvT[0]!=-1)feasibleIntvs.push_back(intvT);
		}
	}

	if (feasibleIntvs.size()==0) return 0; //这意味着整段时间都有碰撞
	//无碰撞发生的并，剩下的就是有碰撞发生的
	std::sort(feasibleIntvs.begin(), feasibleIntvs.end(), 
		[](const Array2d& intv1, const Array2d& intv2){
			return (intv1(0)<intv2(0));
		});
	if(feasibleIntvs[0](0)>0) return 0;
	double minT = feasibleIntvs[0](1);
	for(int i=1;i<feasibleIntvs.size();i++)
		if(feasibleIntvs[i](0)<minT) //不能加等，因为无碰撞给的是开区间，如果有),(的情况加等号会把这个情况漏掉
			minT=std::max(minT, feasibleIntvs[i](1));
		else break;
	
	if(minT<upperTime)return minT;
	else return -1;
}

static double VFTest(const Point &CpPos1, const Point &CpVel1, 
						const Face &CpPos2, const Face &CpVel2,
						Array2d& uv2, 
						const double upperTime = DeltaT) {
	struct PatchPair{
		TriParamBound pb2;
		double tLower;
		PatchPair(const TriParamBound& c2, 
				double t = std::numeric_limits<double>::infinity()): pb2(c2), tLower(t) {}
		bool operator<(PatchPair const &o) const { return tLower > o.tLower; }
		double calcL1Dist(const Point &CpPos1, const Point &CpVel1, 
						const Face &CpPos2, const Face &CpVel2) const{
			auto const ptPos2 = CpPos2.divideBezierPatch(pb2);
			double d=0;
			for(int axis=0;axis<3;axis++){
				double maxv = ptPos2[0][axis], minv=maxv;
				for(int i = 1; i < Face::cntCp; i++) {
					maxv=std::max(maxv, ptPos2[i][axis]);
					minv=std::min(minv, ptPos2[i][axis]);
				}
				d=std::max(d,maxv-minv);
			}
			return d;
		}
	};

	using steady_clock = std::chrono::steady_clock;
	using duration = std::chrono::duration<double>;
	const auto initialTime = steady_clock::now();

	std::priority_queue<PatchPair> heap;
	TriParamBound initParam2;
	auto const ptPos1 = CpPos1.dividePatch();
	auto const ptVel1 = CpVel1.dividePatch();
	auto const ptPos2 = CpPos2.divideBezierPatch(initParam2);
	auto const ptVel2 = CpVel2.divideBezierPatch(initParam2);
	double colTime = primitiveCheck<Point,Face>(ptPos1, ptVel1, ptPos2, ptVel2, upperTime);
	if (colTime>=0 && colTime<=upperTime)
		heap.emplace(initParam2, colTime);
	// std::cout<<"done!\n";

	while (!heap.empty()) {
		auto const cur = heap.top();
		heap.pop();
		cnt++;
		if (cur.calcL1Dist(CpPos1, CpVel1, CpPos2, CpVel2) < MinL1Dist) {
			uv2 = cur.pb2.centerParam();
			const auto endTime = steady_clock::now();
			return cur.tLower;
		}

		// Divide the current patch into four-to-four pieces
		for (int j = 0; j < 4; j++) {
			TriParamBound divUvB2(cur.pb2.interpSubpatchParam(j));
			auto const ptPos1 = CpPos1.dividePatch();
			auto const ptVel1 = CpVel1.dividePatch();
			auto const ptPos2 = CpPos2.divideBezierPatch(divUvB2);
			auto const ptVel2 = CpVel2.divideBezierPatch(divUvB2);
			double colTime = primitiveCheck<Point,Face>(ptPos1, ptVel1, ptPos2, ptVel2, upperTime);
			if (colTime>=0 && colTime<=upperTime){
				heap.emplace(divUvB2, colTime);
			}
		}
	}
	const auto endTime = steady_clock::now();
	return -1;
}

static double modifiedEETest(const Vector3d &CpPos10, const Vector3d &CpPos11,
						const Vector3d &CpVel10, const Vector3d &CpVel11,
						const Vector3d &CpPos20, const Vector3d &CpPos21,
						const Vector3d &CpVel20, const Vector3d &CpVel21,
						double& u1, double& u2, 
						const double upperTime = DeltaT) {
	struct PatchPair{
		Array2d pb1;
		Array2d pb2;
		double tLower;
		PatchPair(const Array2d& c1, const Array2d& c2, 
				double t = std::numeric_limits<double>::infinity()): pb1(c1), pb2(c2), tLower(t) {}
		bool operator<(PatchPair const &o) const { return tLower > o.tLower; }
		double calcL1Dist(const Vector3d &ptPos10, const Vector3d &ptPos11,
						const Vector3d &ptPos20, const Vector3d &ptPos21) const{
			double d1=std::abs((pb1[1]-pb1[0]))*((ptPos10-ptPos11).cwiseAbs().maxCoeff());
			double d2=std::abs((pb2[1]-pb2[0]))*((ptPos20-ptPos21).cwiseAbs().maxCoeff());
			return std::max(d1,d2);
		}
		double calcDist(const Vector3d &ptPos10, const Vector3d &ptPos11,
						const Vector3d &ptPos20, const Vector3d &ptPos21) const{
			double d1=std::abs((pb1[1]-pb1[0]))*((ptPos10-ptPos11).norm());
			double d2=std::abs((pb2[1]-pb2[0]))*((ptPos20-ptPos21).norm());
			return std::max(d1,d2);
		}
	};

	std::priority_queue<PatchPair> heap;
	Array2d initParam(0,1);
	double colTime = primitiveEECheck(CpPos10, CpPos11, CpVel10, CpVel11, 
										CpPos20, CpPos21, CpVel20, CpVel21, upperTime);
	if (colTime>=0 && colTime<=upperTime)
		heap.emplace(initParam, initParam, colTime);

	while (!heap.empty()) {
		auto const cur = heap.top();
		heap.pop();
		cnt++;
		double mid1 = (cur.pb1[0]+cur.pb1[1])*0.5, mid2 = (cur.pb2[0]+cur.pb2[1])*0.5;
		if (cur.calcDist(CpPos10, CpPos11, CpPos20, CpPos21) < MinDist) {
			u1=mid1, u2=mid2;
			return cur.tLower;
		}

		for (int i = 0; i < 2; i++) {
			Array2d divUvB1[2]={Array2d(cur.pb1[0],mid1), Array2d(mid1, cur.pb1[1])};
			for (int j = 0; j < 2; j++) {
				Array2d divUvB2[2]={Array2d(cur.pb2[0],mid2), Array2d(mid2, cur.pb2[1])};
				auto const ptVel10 = (1-divUvB1[i][0])*CpVel10+divUvB1[i][0]*CpVel11;
				auto const ptVel11 = (1-divUvB1[i][1])*CpVel10+divUvB1[i][1]*CpVel11;
				auto const ptVel20 = (1-divUvB2[j][0])*CpVel20+divUvB2[j][0]*CpVel21;
				auto const ptVel21 = (1-divUvB2[j][1])*CpVel20+divUvB2[j][1]*CpVel21;
				auto const ptPos10 = (1-divUvB1[i][0])*CpPos10+divUvB1[i][0]*CpPos11+cur.tLower*ptVel10;
				auto const ptPos11 = (1-divUvB1[i][1])*CpPos10+divUvB1[i][1]*CpPos11+cur.tLower*ptVel11;
				auto const ptPos20 = (1-divUvB2[j][0])*CpPos20+divUvB2[j][0]*CpPos21+cur.tLower*ptVel20;
				auto const ptPos21 = (1-divUvB2[j][1])*CpPos20+divUvB2[j][1]*CpPos21+cur.tLower*ptVel21;
				colTime = cur.tLower + primitiveEECheck(ptPos10, ptPos11, ptVel10, ptVel11, 
										ptPos20, ptPos21, ptVel20, ptVel21, upperTime-cur.tLower);
				if (colTime>=0 && colTime<=upperTime){
					heap.emplace(divUvB1[i], divUvB2[j], colTime);
				}
			}
		}
	}

	return -1;
}
