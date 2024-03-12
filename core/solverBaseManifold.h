# pragma once
#include "mathOps.h"
template<typename ParamObj1, typename ParamObj2, typename ParamBound1, typename ParamBound2>
class SolverBaseManifold{
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
		for(const auto& pos:posEnd1){
			aabb1[0] = pos.cwiseMin(aabb1[0]);
			aabb1[1] = pos.cwiseMax(aabb1[1]);
		}
		for(const auto& pos:posEnd2){
			aabb2[0] = pos.cwiseMin(aabb2[0]);
			aabb2[1] = pos.cwiseMax(aabb2[1]);
		}
	}
	bool separationCheck(const CCDIntv<ParamBound1, ParamBound2> & r){
		Vector3d aaExtent1 = (r.aabb1[1]-aabb1[0]).cwiseMax(aabb1[1]-r.aabb1[0]),
		aaExtent2 = (r.aabb2[1]-aabb2[0]).cwiseMax(aabb2[1]-r.aabb2[0]);
		if(std::max(aaExtent1.norm(), aaExtent2.norm()) < SeparationEucDist) return false;
		return true;
	}
	bool primitiveCheck(const ParamObj1 &CpPos1, const ParamObj1 &CpVel1, 
							const ParamObj2 &CpPos2, const ParamObj2 &CpVel2,
							const ParamBound1 &divUvB1, const ParamBound2 &divUvB2,
							const Array2d divTime = Array2d(0,DeltaT)) {
		calcPatches(CpPos1, CpVel1,CpPos2, CpVel2, divUvB1, divUvB2, divTime);

		std::vector<Vector3d> axes;
		axes.clear();
		setAxes<ParamObj1, ParamObj2>(posStart1, posStart2, axes);

		for(auto& axis:axes){
			double maxProj1 = -INFT, minProj1 = INFT;
			for(const auto&p:posStart1){
				maxProj1 = std::max(maxProj1, p.dot(axis));
				minProj1 = std::min(minProj1, p.dot(axis));
			}
			for(const auto&p:posEnd1){
				maxProj1 = std::max(maxProj1, p.dot(axis));
				minProj1 = std::min(minProj1, p.dot(axis));
			}
			double maxProj2 = -INFT, minProj2 = INFT;
			for(const auto&p:posStart2){
				maxProj2 = std::max(maxProj2, p.dot(axis));
				minProj2 = std::min(minProj2, p.dot(axis));
			}
			for(const auto&p:posEnd2){
				maxProj2 = std::max(maxProj2, p.dot(axis));
				minProj2 = std::min(minProj2, p.dot(axis));
			}
			if(maxProj2<minProj1 || maxProj1<minProj2) return false;
		}
		return true;
	}
	
public:
	double solveCCD(const ParamObj1 &CpPos1, const ParamObj1 &CpVel1, 
						const ParamObj2 &CpPos2, const ParamObj2 &CpVel2,
						std::multiset<CCDIntv<ParamBound1, ParamBound2> > & solutSet,
						const double upperTime = DeltaT,
						const double deltaDist = MinL1Dist) {
		struct PatchPair{
			ParamBound1 pb1;
			ParamBound2 pb2;
			Array2d tIntv;
			PatchPair(const ParamBound1& c1, const ParamBound2& c2, 
					Array2d t = Array2d(0,DeltaT)): pb1(c1), pb2(c2), tIntv(t) {}
			// 按照lb升序排列-->按照ub升序排列
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
		Array2d initTimeIntv(0,upperTime);
		if (primitiveCheck(CpPos1, CpVel1, CpPos2, CpVel2, initParam1, initParam2, initTimeIntv))
			heap.emplace(initParam1, initParam2, initTimeIntv);

		// cnt=1;
		// u得作为上限，不能用下限，会慢一个数量级
		// 我发现所有patch pair测试都是deltaT的上界才能显示出一些框架的bug
		double leastUB = upperTime;//solutSet.empty() ? upperTime : solutSet.rbegin()->t;
		while (!heap.empty()) {
			auto const cur = heap.top();
			heap.pop();
			// cnt++;
			// if(SHOWANS) std::cout<<cnt<<"\n";
			if (cur.tIntv[0] > leastUB + MeantimeEpsilon)
				continue;

			calcAABBs(CpPos1, CpVel1, CpPos2, CpVel2, cur.pb1, cur.pb2, cur.tIntv);
			if(/*cur.tIntv[0] > leastUB - MeantimeEpsilon && */!solutSet.empty()){
				bool discard = false;
				for(const auto& r:solutSet)
					if (!separationCheck(r)){
						discard = true;
						break;
					}
				if (discard) continue;
			}

			if (cur.calcL1Dist(aabb1, aabb2) < MinL1Dist) {
				// if(cur.tIntv[0]<leastUB + MeantimeEpsilon){
					// std::cout<<cur.calcL1Dist(aabb1, aabb2)<<"\n";
					// bool discard = false;
					// if(!solutSet.empty()){
					// 	for(const auto& r:solutSet)
					// 		if (primitiveMaxDist(r) < SeparationEucDist){
					// 			discard = true;
					// 			break;
					// 		}
					// }
					// if(!discard){
					Array2d uv1 = cur.pb1.centerParam();
					Array2d uv2 = cur.pb2.centerParam();
					leastUB = std::min(leastUB, cur.tIntv[1]);
					while(!solutSet.empty() && solutSet.begin()->tIntv[0] > leastUB + MeantimeEpsilon)
						solutSet.erase(solutSet.begin());
					solutSet.insert(CCDIntv<ParamBound1, ParamBound2>(cur.pb1, cur.pb2, cur.tIntv, aabb1, aabb2));

					// bool discard = false;
					// for(const auto& r:solutSet){
					// 	if(DEBUG){
					// 	// std::cout<<primitiveMaxDist(r)<<"\n";
					// 	// std::cout<<aabb1[0].transpose()<<"\n"<<aabb1[1].transpose()<<"\n";
					// 	// std::cout<<r.aabb1[0].transpose()<<"\n"<<r.aabb1[1].transpose()<<"\n";
					// 	// std::cout<<aabb2[0].transpose()<<"\n"<<aabb2[1].transpose()<<"\n";
					// 	// std::cout<<r.aabb2[0].transpose()<<"\n"<<r.aabb2[1].transpose()<<"\n";
					// 	}
					// 	if (primitiveMaxDist(r) < SeparationEucDist){
					// 		discard = true;
					// 		break;
					// 	}
					// }
					// // if(DEBUG){
					// // std::cout<<discard<<'\n';
					// // std::cin.get();}
					// if (discard) continue;
					// else 
					// solutSet.insert(CCDRoot(uv1, uv2, cur.tIntv[0], aabb1, aabb2));
					// }
				// }
				continue;
			}

			// Divide the current patch into two sets of four-to-four pieces
			double tMid = (cur.tIntv[0]+cur.tIntv[1])*0.5;
			Array2d divTime1(cur.tIntv[0],tMid), divTime2(tMid, cur.tIntv[1]);
			for (int i = 0; i < 4; i++) {
				ParamBound1 divUvB1(cur.pb1.interpSubpatchParam(i));
				for (int j = 0; j < 4; j++) {
					ParamBound2 divUvB2(cur.pb2.interpSubpatchParam(j));
					calcPatches(CpPos1, CpVel1, CpPos2, CpVel2, divUvB1, divUvB2, divTime1);
					if (divTime1[0]<leastUB + MeantimeEpsilon && primitiveCheck(CpPos1, CpVel1, CpPos2, CpVel2, divUvB1, divUvB2, divTime1)){
						heap.emplace(divUvB1, divUvB2, divTime1);
					}
					if (divTime2[0]<leastUB + MeantimeEpsilon && primitiveCheck(CpPos1, CpVel1, CpPos2, CpVel2, divUvB1, divUvB2, divTime2)){
						heap.emplace(divUvB1, divUvB2, divTime2);
					}
				}
			}
		}

		const auto endTime = steady_clock::now();
		if(SHOWANS)
			std::cout << "used seconds: " <<
				duration(endTime - initialTime).count()
				<< std::endl;
		std::cout<<leastUB<<"\n";
		if(solutSet.empty()) return -1;
		else return leastUB;
	}
};