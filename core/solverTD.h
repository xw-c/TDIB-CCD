# pragma once
#include "mathOps.h"
#include "triBezier.h"
#include "recBezier.h"
#include "recRationalBezier.h"
template<typename ParamObj1, typename ParamObj2, typename ParamBound1, typename ParamBound2>
class SolverTD{
public:
	static bool primitiveCheck(const ParamObj1 &CpPos1, const ParamObj1 &CpVel1, 
						const ParamObj2 &CpPos2, const ParamObj2 &CpVel2,
						const ParamBound1 &divUvB1, const ParamBound2 &divUvB2,
						double& colTime,
						const double lowerTime = 0, double upperTime = DeltaT) {
		auto ptPos1 = CpPos1.divideBezierPatch(divUvB1);
		auto ptVel1 = CpVel1.divideBezierPatch(divUvB1);
		auto ptPos2 = CpPos2.divideBezierPatch(divUvB2);
		auto ptVel2 = CpVel2.divideBezierPatch(divUvB2);
		for(int i=0;i<ParamObj1::cntCp;i++){
			ptPos1[i]+=ptVel1[i]*lowerTime;
		}
		for(int i=0;i<ParamObj2::cntCp;i++){
			ptPos2[i]+=ptVel2[i]*lowerTime;
		}
		upperTime-=lowerTime;
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

		//无碰撞发生的并，剩下的就是有碰撞发生的
		std::sort(feasibleIntvs.begin(), feasibleIntvs.end(), 
			[](const Array2d& intv1, const Array2d& intv2){
				return (intv1(0)<intv2(0));
			});
		// if(feasibleIntvs[0](0)>0) return lowerTime;
		if (feasibleIntvs.size()==0||feasibleIntvs[0](0)>0) {
			//这意味着整段时间都有碰撞
			colTime = lowerTime;
			return true; 
		}
		// for(const auto&l:feasibleIntvs)std::cout<<"intv:"<<l.transpose()<<"\n";
		double minT = feasibleIntvs[0](1);
		for(int i=1;i<feasibleIntvs.size();i++)
			if(feasibleIntvs[i](0)<minT) //不能加等，因为无碰撞给的是开区间，如果有),(的情况加等号会把这个情况漏掉
				minT=std::max(minT, feasibleIntvs[i](1));
			else break;
		// std::cout<<minT<<"\n";
		if(minT<upperTime){colTime= minT+lowerTime;return true;}
		else {colTime = -1;return false;}
	}
						
	static double solveCCD(const ParamObj1 &CpPos1, const ParamObj1 &CpVel1, 
						const ParamObj2 &CpPos2, const ParamObj2 &CpVel2,
						Array2d& uv1, Array2d& uv2, 
						const double upperTime = DeltaT) {
		struct PatchPair{
			ParamBound1 pb1;
			ParamBound2 pb2;
			double tLower;
			PatchPair(const ParamBound1& c1, const ParamBound2& c2, 
					double t = std::numeric_limits<double>::infinity()): pb1(c1), pb2(c2), tLower(t) {}
			bool operator<(PatchPair const &o) const { return tLower > o.tLower; }
			double calcL1Dist(const ParamObj1 &CpPos1, const ParamObj1 &CpVel1, 
							const ParamObj2 &CpPos2, const ParamObj2 &CpVel2) const{
				auto const ptPos1 = CpPos1.divideBezierPatch(pb1);
				auto const ptPos2 = CpPos2.divideBezierPatch(pb2);
				double d1=calcAAExtent<ParamObj1>(ptPos1);
				double d2=calcAAExtent<ParamObj1>(ptPos2);
				return std::max(d1, d2);
			}
		};

		using steady_clock = std::chrono::steady_clock;
		using duration = std::chrono::duration<double>;
		const auto initialTime = steady_clock::now();

		std::priority_queue<PatchPair> heap;
		ParamBound1 initParam1;
		ParamBound2 initParam2;
		double colTime;
		if (primitiveCheck(CpPos1, CpVel1, CpPos2, CpVel2, initParam1, initParam2, colTime, 0, upperTime))
			heap.emplace(initParam1, initParam2, colTime);
		// cnt=1;
		while (!heap.empty()) {
			auto const cur = heap.top();
			heap.pop();
			// cnt++;
			// if(SHOWANS) std::cout<<cnt<<"\n";

			// Decide whether the algorithm converges
			if (cur.calcL1Dist(CpPos1, CpVel1, CpPos2, CpVel2) < MinL1Dist) {
				uv1 = cur.pb1.centerParam();
				uv2 = cur.pb2.centerParam();
				const auto endTime = steady_clock::now();
				if(SHOWANS)
					std::cout << "min time: "<<  cur.tLower 
						<< "\nused seconds: " << duration(endTime - initialTime).count()
						<< std::endl;
				return cur.tLower;
			}

			// Divide the current patch into four-to-four pieces
			for (int i = 0; i < 4; i++) {
				ParamBound1 divUvB1(cur.pb1.interpSubpatchParam(i));
				for (int j = 0; j < 4; j++) {
					ParamBound2 divUvB2(cur.pb2.interpSubpatchParam(j));
					if (primitiveCheck(CpPos1, CpVel1, CpPos2, CpVel2, divUvB1, divUvB2, colTime, cur.tLower, upperTime)){
						heap.emplace(divUvB1, divUvB2, colTime);
					}
				}
			}
		}

		const auto endTime = steady_clock::now();
		if(SHOWANS)
			std::cout << "used seconds: " << duration(endTime - initialTime).count()
				<< std::endl;
		return -1;
	}
};
