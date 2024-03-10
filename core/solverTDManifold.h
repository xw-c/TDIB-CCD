# pragma once
#include "mathOps.h"
template<typename ParamObj1, typename ParamObj2, typename ParamBound1, typename ParamBound2>
class SolverTDManifold{
	std::array<Vector3d, ParamObj1::cntCp> posStart1, posEnd1;
	std::array<Vector3d, ParamObj2::cntCp> posStart2, posEnd2;
	std::array<Vector3d, 2> aabb1, aabb2;

	void calcPatches(const ParamObj1 &CpPos1, const ParamObj1 &CpVel1, 
							const ParamObj2 &CpPos2, const ParamObj2 &CpVel2,
							const ParamBound1 &divUvB1, const ParamBound2 &divUvB2,
							const Array2d divTime = Array2d(0,DeltaT)) {
		// init patches
		posStart1 = CpPos1.divideBezierPatch(divUvB1), posEnd1 = posStart1;
		posStart2 = CpPos2.divideBezierPatch(divUvB2), posEnd2 = posStart2;
		auto ptVel1 = CpVel1.divideBezierPatch(divUvB1);
		auto ptVel2 = CpVel2.divideBezierPatch(divUvB2);
		for(int i=0;i<ParamObj1::cntCp;i++){
			posStart1[i]+=ptVel1[i]*divTime[0],
			posEnd1[i]+=ptVel1[i]*divTime[1];
		}
		for(int i=0;i<ParamObj2::cntCp;i++){
			posStart2[i]+=ptVel2[i]*divTime[0],
			posEnd2[i]+=ptVel2[i]*divTime[1];
		}
	}
	void calcAABBs(const ParamObj1 &CpPos1, const ParamObj1 &CpVel1, 
							const ParamObj2 &CpPos2, const ParamObj2 &CpVel2,
							const ParamBound1 &divUvB1, const ParamBound2 &divUvB2,
							const Array2d divTime = Array2d(0,DeltaT)){
		calcPatches(CpPos1, CpVel1,CpPos2, CpVel2, divUvB1, divUvB2, divTime);
		
		aabb1[0] = aabb2[0] = Vector3d::Constant(INFT),
		aabb1[1] = aabb2[1] = Vector3d::Constant(-INFT); 
		for(const auto& pos:posStart1){
			aabb1[0] = pos.cwiseMin(aabb1[0]);
			aabb1[1] = pos.cwiseMax(aabb1[1]);
		}
		for(const auto& pos:posStart2){
			aabb2[0] = pos.cwiseMin(aabb2[0]);
			aabb2[1] = pos.cwiseMax(aabb2[1]);
		}
	}
	double primitiveMaxDist(const CCDRoot& r){
		Vector3d aaExtent1 = (r.aabb1[1]-aabb1[0]).cwiseMax(aabb1[1]-r.aabb1[0]),
		aaExtent2 = (r.aabb2[1]-aabb2[0]).cwiseMax(aabb2[1]-r.aabb2[0]);
		return std::max(aaExtent1.maxCoeff(), aaExtent2.maxCoeff());
	}

public:
	bool primitiveCheck(const ParamObj1 &CpPos1, const ParamObj1 &CpVel1, 
						const ParamObj2 &CpPos2, const ParamObj2 &CpVel2,
						const ParamBound1 &divUvB1, const ParamBound2 &divUvB2,
						Array2d& colTime,
						const Array2d& timeIntv = Array2d(0,DeltaT)) {
		auto ptPos1 = CpPos1.divideBezierPatch(divUvB1);
		auto ptVel1 = CpVel1.divideBezierPatch(divUvB1);
		auto ptPos2 = CpPos2.divideBezierPatch(divUvB2);
		auto ptVel2 = CpVel2.divideBezierPatch(divUvB2);
		for(int i=0;i<ParamObj1::cntCp;i++){
			ptPos1[i]+=ptVel1[i]*timeIntv[0];
		}
		for(int i=0;i<ParamObj2::cntCp;i++){
			ptPos2[i]+=ptVel2[i]*timeIntv[0];
		}
		const double upperTime = timeIntv[1] - timeIntv[0];
		std::vector<Vector3d> axes;
		setAxes<ParamObj1, ParamObj2>(ptPos1, ptPos2, axes);
		

		std::vector<Array2d> feasibleIntvs;
		feasibleIntvs.clear();

		auto AxisCheck=[&](std::vector<Line> lines1, std::vector<Line> lines2){
			std::vector<Line> ch1, ch2;
			std::vector<double> pts1, pts2;
			ch1.clear(); ch2.clear();
			pts1.clear(); pts2.clear();

			getCH(lines1, ch1, pts1, true, upperTime);
			getCH(lines2, ch2, pts2, false, upperTime);
			const auto intvT = linearCHIntersect(ch1, ch2, pts1, pts2, upperTime);
			if(intvT[0]!=-1)feasibleIntvs.push_back(intvT);
		};

		for(const auto& axis:axes){
			std::vector<Line> ptLines1, ptLines2;
			ptLines1.clear(); ptLines2.clear();
			for(int i = 0; i < ParamObj1::cntCp; i++) ptLines1.emplace_back(ptVel1[i].dot(axis), ptPos1[i].dot(axis));
			for(int i = 0; i < ParamObj2::cntCp; i++) ptLines2.emplace_back(ptVel2[i].dot(axis), ptPos2[i].dot(axis));
			std::sort(ptLines1.begin(), ptLines1.end());
			std::sort(ptLines2.begin(), ptLines2.end());
			AxisCheck(ptLines1, ptLines2);
			AxisCheck(ptLines2, ptLines1);
		}
		// if (feasibleIntvs.size()==0) return 0; //这意味着整段时间都有碰撞

		// if(feasibleIntvs[0](0)>0) return timeIntv[0];
		if (feasibleIntvs.size()==0) {
			//这意味着整段时间都有碰撞
			colTime = timeIntv;
			return true; 
		}

		//无碰撞发生的并，剩下的就是有碰撞发生的
		double minT = 0, maxT = upperTime;
		std::sort(feasibleIntvs.begin(), feasibleIntvs.end(), 
			[](const Array2d& intv1, const Array2d& intv2){
				return (intv1(0)<intv2(0));
			});
		// for(const auto&l:feasibleIntvs)std::cout<<"intv:"<<l.transpose()<<"\n";
		if(feasibleIntvs[0](0)<=0){
			minT = feasibleIntvs[0](1);
			for(int i=1;i<feasibleIntvs.size();i++)
				if(feasibleIntvs[i](0)<minT) //不能加等，因为无碰撞给的是开区间，如果有),(的情况加等号会把这个情况漏掉
					minT=std::max(minT, feasibleIntvs[i](1));
				else break;
		}
		if(minT >= upperTime){ colTime = Array2d(-1,-1); return false; }
		
		std::sort(feasibleIntvs.begin(), feasibleIntvs.end(), 
			[](const Array2d& intv1, const Array2d& intv2){
				return (intv1(1)>intv2(1));
			});
		// for(const auto&l:feasibleIntvs)std::cout<<"intv:"<<l.transpose()<<"\n";
		if(feasibleIntvs[0](1)>=upperTime){
			maxT = feasibleIntvs[0](0);
			for(int i=1;i<feasibleIntvs.size();i++)
				if(feasibleIntvs[i](1)>maxT) //不能加等，因为无碰撞给的是开区间，如果有),(的情况加等号会把这个情况漏掉
					maxT=std::min(maxT, feasibleIntvs[i](0));
				else break;
		}
		if(minT >= upperTime){ colTime = Array2d(-1,-1); return false; }

		// std::cout<<minT<<"\n";
		colTime = Array2d(minT + timeIntv[0], maxT + timeIntv[0]); 
		return true;
	}
	double solveCCD(const ParamObj1 &CpPos1, const ParamObj1 &CpVel1, 
						const ParamObj2 &CpPos2, const ParamObj2 &CpVel2,
						std::multiset<CCDRoot> & solutSet,
						const double upperTime = DeltaT,
						const double deltaDist = MinL1Dist) {
		struct PatchPair{
			ParamBound1 pb1;
			ParamBound2 pb2;
			Array2d tIntv;
			PatchPair(const ParamBound1& c1, const ParamBound2& c2, 
					Array2d t = Array2d(0,DeltaT)): pb1(c1), pb2(c2), tIntv(t) {}
			bool operator<(PatchPair const &o) const { return tIntv[1] > o.tIntv[1]; }
			double calcL1Dist(const std::array<Vector3d, 2> &aabb1, 
							const std::array<Vector3d, 2> &aabb2) const{
				double d1=(aabb1[1]-aabb1[0]).maxCoeff();
				double d2=(aabb2[1]-aabb2[0]).maxCoeff();
				return std::max(d1, d2);
			}
		};

		using steady_clock = std::chrono::steady_clock;
		using duration = std::chrono::duration<double>;
		const auto initialTime = steady_clock::now();

		std::priority_queue<PatchPair> heap;
		ParamBound1 initParam1;
		ParamBound2 initParam2;
		Array2d initTimeIntv(0,upperTime), colTime;
		if (primitiveCheck(CpPos1, CpVel1, CpPos2, CpVel2, initParam1, initParam2, colTime, initTimeIntv))
			heap.emplace(initParam1, initParam2, colTime);
		// if (primitiveCheck(CpPos1, CpVel1, CpPos2, CpVel2, initParam1, initParam2, colTime, 0, upperTime))
		// 	heap.emplace(initParam1, initParam2, colTime);

		double leastUb = upperTime;
		while (!heap.empty()) {
			auto const cur = heap.top();
			heap.pop();
			// cnt++;
			// if(SHOWANS) std::cout<<cnt<<"\n";
			if (cur.tIntv[0] > leastUb + MeantimeEpsilon)
				continue;

			calcAABBs(CpPos1, CpVel1, CpPos2, CpVel2, cur.pb1, cur.pb2, cur.tIntv);
			if(cur.tIntv[1] > leastUb - MeantimeEpsilon && !solutSet.empty()){
				bool discard = false;
				for(const auto& r:solutSet)
					if (primitiveMaxDist(r) < SeparationDist){
						discard = true;
						break;
					}
				if (discard) continue;
			}

			if (cur.calcL1Dist(aabb1, aabb2) < deltaDist) {
					Array2d uv1 = cur.pb1.centerParam();
					Array2d uv2 = cur.pb2.centerParam();
					leastUb = std::min(leastUb, cur.tIntv[0]);
					while(!solutSet.empty() && solutSet.begin()->t > leastUb + MeantimeEpsilon)
						solutSet.erase(solutSet.begin());
					solutSet.insert(CCDRoot(uv1, uv2, aabb1, aabb2, cur.tIntv[0]));
				continue;
			}

			// Divide the current patch into two sets of four-to-four pieces
			for (int i = 0; i < 4; i++) {
				ParamBound1 divUvB1(cur.pb1.interpSubpatchParam(i));
				for (int j = 0; j < 4; j++) {
					ParamBound2 divUvB2(cur.pb2.interpSubpatchParam(j));
					if (cur.tIntv[0] < leastUb + MeantimeEpsilon && primitiveCheck(CpPos1, CpVel1, CpPos2, CpVel2, divUvB1, divUvB2, colTime, cur.tIntv)){
						heap.emplace(divUvB1, divUvB2, colTime);
					}
				}
			}
		}

		const auto endTime = steady_clock::now();
		if(SHOWANS)
			std::cout << "used seconds: " <<
				duration(endTime - initialTime).count()
				<< std::endl;
		if(solutSet.empty()) return -1;
		else return solutSet.rbegin()->t;
	}
};