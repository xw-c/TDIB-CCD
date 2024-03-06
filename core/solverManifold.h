# pragma once
#include "mathOps.h"
template<typename ParamObj1, typename ParamObj2, typename ParamBound1, typename ParamBound2>
class SolverManifold{
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
	// static bool primitiveCheck(const ParamObj1 &CpPos1, const ParamObj1 &CpVel1, 
	// 						const ParamObj2 &CpPos2, const ParamObj2 &CpVel2,
	// 						const ParamBound1 &divUvB1, const ParamBound2 &divUvB2,
	// 						const Array2d divTime = Array2d(0,DeltaT)) {
	// 	initInclusions(CpPos1, CpVel1,CpPos2, CpVel2, divUvB1, divUvB2, divTime);
	// 	if(bbType==BoundingBoxType::AABB){
	// 		for(int i=0;i<3;i++){
	// 			if(bb2[1][i] < bb1[0][i] || bb1[1][i] < bb2[0][i]) return false;
	// 		}
	// 		return true;
	// 	}
	// 	else if(bbType==BoundingBoxType::OBB){
	// 		for(int i=0;i<3;i++){
	// 			double proj[2] = {bb1[0].dot(axes2[i]), bb1[1].dot(axes2[i])};
	// 			if(std::max(proj[0], proj[1]) < bb2[0][i] || std::min(proj[0], proj[1]) < bb2[1][i])
	// 				return false;
	// 		}
	// 		for(int i=0;i<3;i++){
	// 			double proj[2] = {bb2[0].dot(axes1[i]), bb2[1].dot(axes1[i])};
	// 			if(std::max(proj[0], proj[1]) < bb1[0][i] || std::min(proj[0], proj[1]) < bb1[1][i])
	// 				return false;
	// 		}
	// 		for(int i=0;i<3;i++)
	// 			for(int j=0;j<3;j++){
	// 				Vector3d axis = axes1[i].cross(axes2[i]);
	// 				double proj1[2] = {bb1[0].dot(axis), bb1[1].dot(axis)},
	// 				proj2[2] = {bb2[0].dot(axis), bb2[1].dot(axis)};
	// 				double maxProj1 = std::max(proj1[0], proj1[1]),
	// 				minProj1 = std::min(proj1[0], proj1[1]),
	// 				maxProj2 = std::max(proj2[0], proj2[1]),
	// 				minProj2 = std::min(proj2[0], proj2[1]);
	// 				if(maxProj2<minProj1 || maxProj1<minProj2) return false;
	// 			}
	// 		return true;
	// 	}
	// }
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
		if (primitiveCheck(CpPos1, CpVel1, CpPos2, CpVel2, initParam1, initParam2))
			heap.emplace(initParam1, initParam2);

		// cnt=1;
		double leastUb = upperTime;//solutSet.empty() ? upperTime : solutSet.rbegin()->t;
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
				if(cur.tIntv[1]<leastUb + MeantimeEpsilon){
					Array2d uv1 = cur.pb1.centerParam();
					Array2d uv2 = cur.pb2.centerParam();
					leastUb = std::min(leastUb, cur.tIntv[1]);
					while(!solutSet.empty() && solutSet.begin()->t > leastUb + MeantimeEpsilon)
						solutSet.erase(solutSet.begin());

					// bool discard = false;
					// for(const auto& r:solutSet){
					// 	if(DEBUG){
					// 	// std::cout<<primitiveMaxDist(r)<<"\n";
					// 	// std::cout<<aabb1[0].transpose()<<"\n"<<aabb1[1].transpose()<<"\n";
					// 	// std::cout<<r.aabb1[0].transpose()<<"\n"<<r.aabb1[1].transpose()<<"\n";
					// 	// std::cout<<aabb2[0].transpose()<<"\n"<<aabb2[1].transpose()<<"\n";
					// 	// std::cout<<r.aabb2[0].transpose()<<"\n"<<r.aabb2[1].transpose()<<"\n";
					// 	}
					// 	if (primitiveMaxDist(r) < SeparationDist){
					// 		discard = true;
					// 		break;
					// 	}
					// }
					// // if(DEBUG){
					// // std::cout<<discard<<'\n';
					// // std::cin.get();}
					// if (discard) continue;
					// else 
					solutSet.insert(CCDRoot(uv1, uv2, aabb1, aabb2, cur.tIntv[0]));
				}
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
					if (divTime1[0]<leastUb + MeantimeEpsilon && primitiveCheck(CpPos1, CpVel1, CpPos2, CpVel2, divUvB1, divUvB2, divTime1)){
						heap.emplace(divUvB1, divUvB2, divTime1);
					}
					if (divTime2[0]<leastUb + MeantimeEpsilon && primitiveCheck(CpPos1, CpVel1, CpPos2, CpVel2, divUvB1, divUvB2, divTime2)){
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
		if(solutSet.empty()) return -1;
		else return solutSet.rbegin()->t;
	}
};