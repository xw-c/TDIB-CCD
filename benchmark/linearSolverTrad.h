# pragma once
#include "linearGeom.h"

class LinearSolverTrad{
public:
	static bool primitiveEECheck(const std::array<Vector3d,Edge::cntCp> &ptPos1, 
						const std::array<Vector3d,Edge::cntCp> &ptVel1, 
						const std::array<Vector3d,Edge::cntCp> &ptPos2, 
						const std::array<Vector3d,Edge::cntCp> &ptVel2,
						const BoundingBoxType& bb,
						const Array2d& initTimeIntv) {
		
		std::array<Vector3d,Edge::cntCp> posStart1, posEnd1, posStart2, posEnd2;
		for(int i=0;i<Edge::cntCp;i++){
			posStart1[i]=ptPos1[i]+ptVel1[i]*initTimeIntv[0],
			posEnd1[i]=ptPos1[i]+ptVel1[i]*initTimeIntv[1];
		}
		for(int i=0;i<Edge::cntCp;i++){
			posStart2[i]=ptPos2[i]+ptVel2[i]*initTimeIntv[0],
			posEnd2[i]=ptPos2[i]+ptVel2[i]*initTimeIntv[1];
		}

		std::vector<Vector3d> axes;
		axes.clear();
		if(bb==BoundingBoxType::AABB){
			axes = {Vector3d::Unit(0), Vector3d::Unit(1), Vector3d::Unit(2)};
		}
		else if(bb==BoundingBoxType::OBB){
			Vector3d lu = Edge::direction(ptPos1) + initTimeIntv[0]*Edge::direction(ptVel1);
			Vector3d lvtmp = Edge::direction(ptPos2) + initTimeIntv[0]*Edge::direction(ptVel2);
			Vector3d randu = Vector3d::Random()*lu.norm()*0.01;
			Vector3d randv = Vector3d::Random()*lvtmp.norm()*0.01;
			lu+=randu;
			lvtmp+=randv;
			Vector3d ln = lu.cross(lvtmp);
			Vector3d lv = ln.cross(lu);
			axes = {lu, lv, ln};
		}
		for(auto& axis:axes){
			double maxProj1 = -std::numeric_limits<double>::infinity(), minProj1 = std::numeric_limits<double>::infinity();
			for(const auto&p:posStart1){
				maxProj1 = std::max(maxProj1, p.dot(axis));
				minProj1 = std::min(minProj1, p.dot(axis));
			}
			for(const auto&p:posEnd1){
				maxProj1 = std::max(maxProj1, p.dot(axis));
				minProj1 = std::min(minProj1, p.dot(axis));
			}
			double maxProj2 = -std::numeric_limits<double>::infinity(), minProj2 = std::numeric_limits<double>::infinity();
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

	static double solveEETest(const Edge &CpPos1, const Edge &CpVel1, 
						const Edge &CpPos2, const Edge &CpVel2,
						double& u1, double& u2, 
						const BoundingBoxType& bb,
						const double upperTime,
						const double deltaDist) {
		std::priority_queue<EEPair> heap;
		Array2d initParam(0,1), initTimeIntv(0,upperTime);
		if (primitiveEECheck(CpPos1.ctrlp, CpVel1.ctrlp, CpPos2.ctrlp, CpVel2.ctrlp, bb, initTimeIntv))
			heap.emplace(initParam, initParam, initTimeIntv);
		
		while (!heap.empty()) {
			auto const cur = heap.top();
			heap.pop();

			// Decide whether the algorithm converges
			double mid1 = (cur.pb1[0]+cur.pb1[1])*0.5, mid2 = (cur.pb2[0]+cur.pb2[1])*0.5;
			if (cur.calc4dWidth() < deltaDist) {
				u1=mid1, u2=mid2;
				return cur.tIntv[0];
			}

			// Divide the current patch into four-to-four pieces
			double tMid = (cur.tIntv[0]+cur.tIntv[1])*0.5;
			Array2d divTime1(cur.tIntv[0],tMid), divTime2(tMid, cur.tIntv[1]);
			for (int i = 0; i < 2; i++) {
				Array2d divUvB1[2]={Array2d(cur.pb1[0],mid1), Array2d(mid1, cur.pb1[1])};
				for (int j = 0; j < 2; j++) {
					Array2d divUvB2[2]={Array2d(cur.pb2[0],mid2), Array2d(mid2, cur.pb2[1])};
					auto ptPos1 = CpPos1.divideBezierPatch(divUvB1[i]);
					auto ptVel1 = CpVel1.divideBezierPatch(divUvB1[i]);
					auto ptPos2 = CpPos2.divideBezierPatch(divUvB2[j]);
					auto ptVel2 = CpVel2.divideBezierPatch(divUvB2[j]);
					if (primitiveEECheck(ptPos1, ptVel1, ptPos2, ptVel2, bb, divTime1)){
						heap.emplace(divUvB1[i], divUvB2[j], divTime1);
					}
					if (primitiveEECheck(ptPos1, ptVel1, ptPos2, ptVel2, bb, divTime2)){
						heap.emplace(divUvB1[i], divUvB2[j], divTime2);
					}
				}
			}
		}
		return -1;
	}

	static bool primitiveVFCheck(const Vector3d &ptPos1, const Vector3d &ptVel1, 
						const std::array<Vector3d,Face::cntCp> &ptPos2, 
						const std::array<Vector3d,Face::cntCp> &ptVel2,
						const BoundingBoxType& bb,
						const Array2d& initTimeIntv) {
		auto posStart1 = ptPos1+initTimeIntv[0]*ptVel1;
		auto posEnd1 = ptPos1+initTimeIntv[1]*ptVel1;
		std::array<Vector3d,Face::cntCp> posStart2, posEnd2;
		for(int i=0;i<Face::cntCp;i++){
			posStart2[i]=ptPos2[i]+ptVel2[i]*initTimeIntv[0],
			posEnd2[i]=ptPos2[i]+ptVel2[i]*initTimeIntv[1];
		}

		std::vector<Vector3d> axes;
		if(bb==BoundingBoxType::AABB){
			axes = {Vector3d::Unit(0), Vector3d::Unit(1), Vector3d::Unit(2)};
		}
		else if(bb==BoundingBoxType::OBB){
			Vector3d lu = Face::axisU(ptPos2) + initTimeIntv[0]*Face::axisU(ptVel2);
			Vector3d lvtmp = Face::axisV(ptPos2) + initTimeIntv[0]*Face::axisV(ptVel2);
			lu[0]*=1.01;
			Vector3d ln = lu.cross(lvtmp);
			Vector3d lv = ln.cross(lu);
			axes = {lu, lv, ln};
		}

		for(auto& axis:axes){
			double maxProj1 = posStart1.dot(axis), minProj1 = posStart1.dot(axis);
			maxProj1 = std::max(maxProj1,posEnd1.dot(axis)), minProj1 = std::min(minProj1,posEnd1.dot(axis));
			double maxProj2 = -std::numeric_limits<double>::infinity(), minProj2 = std::numeric_limits<double>::infinity();
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

	static double solveVFTest(const Vector3d &CpPos1, const Vector3d &CpVel1, 
						const Face &CpPos2, const Face &CpVel2,
						Array2d& uv,
						const BoundingBoxType& bb,
						const double upperTime,
						const double deltaDist) {
		std::priority_queue<VFPair> heap;
		TriParamBound initParam;
		Array2d initTimeIntv(0,upperTime);
		if (primitiveVFCheck(CpPos1, CpVel1, CpPos2.ctrlp, CpVel2.ctrlp, bb, initTimeIntv))
			heap.emplace(initParam, initTimeIntv);
		
		while (!heap.empty()) {
			auto const cur = heap.top();
			heap.pop();

			// Decide whether the algorithm converges
			if (cur.calc4dWidth() < deltaDist) {
				uv = cur.pb.centerParam();
				return cur.tIntv[0];
			}

			// Divide the current patch into four-to-four pieces
			double tMid = (cur.tIntv[0]+cur.tIntv[1])*0.5;
			Array2d divTime1(cur.tIntv[0],tMid), divTime2(tMid, cur.tIntv[1]);
			for (int j = 0; j < 4; j++) {
				TriParamBound divUvB2(cur.pb.interpSubpatchParam(j));
				auto ptPos2 = CpPos2.divideBezierPatch(divUvB2);
				auto ptVel2 = CpVel2.divideBezierPatch(divUvB2);
				if (primitiveVFCheck(CpPos1, CpVel1, ptPos2, ptVel2, bb, divTime1)){
					heap.emplace(divUvB2, divTime1);
				}
				if (primitiveVFCheck(CpPos1, CpVel1, ptPos2, ptVel2, bb, divTime2)){
					heap.emplace(divUvB2, divTime2);
				}
			}
		}
		return -1;
	}
};
