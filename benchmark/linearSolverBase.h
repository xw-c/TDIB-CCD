# pragma once
#include "linearGeom.h"

class LinearSolverBase{
private:
	static Array2d axisCheck(std::vector<Line> lines1, std::vector<Line> lines2,
						const Array2d& timeIntv){
		std::vector<Line> ch1, ch2;
		ch1.clear(); ch2.clear();

			// for(const auto& l:lines1)std::cout<<"lines1:  "<<l.k<<" "<<l.b<<"\n";
			// for(const auto& l:lines2)std::cout<<"lines2:  "<<l.k<<" "<<l.b<<"\n";
		robustCH(lines1, ch1, true, timeIntv);
		robustCH(lines2, ch2, false, timeIntv);
			// for(const auto& l:ch1)std::cout<<"ch1:  "<<l.k<<" "<<l.b<<"\n";
			// for(const auto& l:ch2)std::cout<<"ch2:  "<<l.k<<" "<<l.b<<"\n";
		return robustHullIntersect(ch1, ch2, timeIntv);
	};
	static void robustCH(std::vector<Line>& lines, std::vector<Line>& ch, 
				const bool getMaxCH, const Array2d& tIntv) {
		if(!getMaxCH)std::reverse(lines.begin(),lines.end());
		lines.erase(std::unique(lines.begin(), lines.end()), lines.end()); // 去重
		ch.clear();
		ch.push_back(lines[0]);
		int alpha = 1;
		while(alpha < lines.size()){
			int beta = ch.size()-1;
			while(beta > 0){
				double chfp = (ch[beta].k-ch[beta-1].k)*(lines[alpha].b-ch[beta-1].b)
							-(lines[alpha].k-ch[beta-1].k)*(ch[beta].b-ch[beta-1].b);
				if(chfp>=0){
					ch.pop_back();
					beta--;
				}
				else break;
				
			}
			if(beta==0){
				double chStart = tIntv[0]*(lines[alpha].k-ch[0].k)+(lines[alpha].b-ch[0].b);
				if((getMaxCH&&chStart>=0)||(!getMaxCH&&chStart<=0))
					ch.pop_back();
			}
			if(ch.empty())ch.push_back(lines[alpha]);
			else{
				double chEnd = tIntv[1]*(lines[alpha].k-ch[beta].k)+(lines[alpha].b-ch[beta].b);
				if((getMaxCH&&chEnd>0)||(!getMaxCH&&chEnd<0))
					ch.push_back(lines[alpha]);
			}
			alpha++;
		}
		if(ch.empty()){
			std::cerr<<"empty CH!\n";
			exit(-1);
		}
	}
	static Array2d robustHullIntersect(const std::vector<Line>& ch1, const std::vector<Line>& ch2, 
							const Array2d& tIntv) {
		int id1=0, id2=0;
		double intvL=-1, intvR=-1;
		if(ch1[0].k*tIntv[0]+ch1[0].b<ch2[0].k*tIntv[0]+ch2[0].b)intvL=tIntv[0];
		else{
			while(id1<ch1.size()&&id2<ch2.size()){
				if(ch1[id1].k>=ch2[id2].k){
					break;
				}
				double hifp1, hifp2;
				if(id1<ch1.size()-1)
					hifp1=(ch1[id1+1].k-ch2[id2].k)*(ch1[id1].b-ch2[id2].b)
							-(ch1[id1].k-ch2[id2].k)*(ch1[id1+1].b-ch2[id2].b);
				else 
					hifp1=tIntv[1]*(ch1[id1].k-ch2[id2].k)+(ch1[id1].b-ch2[id2].b);
				if(id2<ch2.size()-1)
					hifp2=(ch1[id1].k-ch2[id2+1].k)*(ch1[id1].b-ch2[id2].b)
							-(ch1[id1].k-ch2[id2].k)*(ch1[id1].b-ch2[id2+1].b);
				else
					hifp2=tIntv[1]*(ch1[id1].k-ch2[id2].k)+(ch1[id1].b-ch2[id2].b);
				if(hifp1<0){
					if(hifp2<0){
						intvL = -(ch1[id1].b-ch2[id2].b)/(ch1[id1].k-ch2[id2].k);
						break;
					}
					else id2++;
				}
				else {
					id1++;
					if(hifp2<0);
					else id2++;
				}
			}
			if(intvL==-1||intvL>=tIntv[1])return Array2d(-1,-1);
		}

		id1 = ch1.size()-1, id2 = ch2.size()-1;
		if((ch1[id1].k-ch2[id2].k)*tIntv[1]+(ch1[id1].b-ch2[id2].b)<0)intvR=tIntv[1];
		else{
			while(id1>=0&&id2>=0){
				if(ch1[id1].k<=ch2[id2].k){
					std::cerr<<"end at strange slopes?\n";
					exit(-1);
				}
				double hifp1, hifp2;
				if(id1>0)
					hifp1=(ch1[id1].k-ch2[id2].k)*(ch1[id1-1].b-ch2[id2].b)
							-(ch1[id1-1].k-ch2[id2].k)*(ch1[id1].b-ch2[id2].b);
				else 
					hifp1=tIntv[0]*(ch1[id1].k-ch2[id2].k)+(ch1[id1].b-ch2[id2].b);
				if(id2>0)
					hifp2=(ch1[id1].k-ch2[id2].k)*(ch1[id1].b-ch2[id2-1].b)
							-(ch1[id1].k-ch2[id2-1].k)*(ch1[id1].b-ch2[id2].b);
				else
					hifp2=tIntv[0]*(ch1[id1].k-ch2[id2].k)+(ch1[id1].b-ch2[id2].b);
				if(hifp1<0){
					if(hifp2<0){
						intvR = -(ch1[id1].b-ch2[id2].b)/(ch1[id1].k-ch2[id2].k);
						break;
					}
					else id2--;
				}
				else {
					id1--;
					if(hifp2<0);
					else id2--;
				}
			}
			if(intvR==-1){
				std::cerr<<"intvL done but no intvR?\n";
				exit(-1);
			}
			if(intvR<=intvL)
				return Array2d(-1,-1);
		}

		intvL = std::max(intvL, tIntv[0]);
		intvR = std::min(intvR, tIntv[1]);

		if(intvL>intvR||intvL<tIntv[0]||intvR>tIntv[1]){
			std::cerr<<"error intersection!\n";
			std::cerr<<intvL<<" "<<intvR<<", in range"<<tIntv[0]<<" "<<tIntv[0]<<"\n";
			exit(-1);
		}
		else return Array2d(intvL,intvR);
	}			
	static bool intvMerge(std::vector<Array2d>& feasibleIntvs, Array2d& colTime, const Array2d& initTimeIntv){
		if (feasibleIntvs.size()==0) {
			colTime = initTimeIntv;
			return true; 
		}
		double minT = initTimeIntv[0], maxT = initTimeIntv[1];
		std::sort(feasibleIntvs.begin(), feasibleIntvs.end(), 
			[](const Array2d& intv1, const Array2d& intv2){
				return (intv1(0)<intv2(0));
			});
		if(feasibleIntvs[0](0)<=initTimeIntv[0]){
			minT = feasibleIntvs[0](1);
			for(int i=1;i<feasibleIntvs.size();i++)
				if(feasibleIntvs[i](0)<minT) //不能加等，因为无碰撞给的是开区间，如果有),(的情况加等号会把这个情况漏掉
					minT=std::max(minT, feasibleIntvs[i](1));
				else break;
		}
		if(minT > initTimeIntv[1]) { colTime = Array2d(-1,-1); return false; }
		
		std::sort(feasibleIntvs.begin(), feasibleIntvs.end(), 
			[](const Array2d& intv1, const Array2d& intv2){
				return (intv1(1)>intv2(1));
			});
		if(feasibleIntvs[0](1)>=initTimeIntv[1]){
			maxT = feasibleIntvs[0](0);
			for(int i=1;i<feasibleIntvs.size();i++)
				if(feasibleIntvs[i](1)>maxT) //不能加等，因为无碰撞给的是开区间，如果有),(的情况加等号会把这个情况漏掉
					maxT=std::min(maxT, feasibleIntvs[i](0));
				else break;
		}
		if(initTimeIntv[0] > maxT) { colTime = Array2d(-1,-1); return false; }
		// maxT=std::max(minT, maxT);
		else if(minT > maxT) { colTime = Array2d(minT, initTimeIntv[1]);}//{ colTime = Array2d(-1,-1); return false; }
		colTime = Array2d(minT, maxT); 
		return true;
	}
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
			//axis变成0即为退化情况，但看起来不用特殊处理
			Vector3d lu = Edge::direction(ptPos1) + initTimeIntv[0]*Edge::direction(ptVel1);//u延展的方向
			Vector3d lvtmp = Edge::direction(ptPos2) + initTimeIntv[0]*Edge::direction(ptVel2);//u延展的方向
			// lu[0]*=1.01;
			Vector3d randu = Vector3d::Random()*lu.norm()*0.01;
			Vector3d randv = Vector3d::Random()*lvtmp.norm()*0.01;
			lu+=randu;
			lvtmp+=randv;
			// lvtmp[0]*=1.5;
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
			// if (cur.calcL1Dist(CpPos1, CpVel1, CpPos2, CpVel2) < deltaDist) {
			if (cur.calcWidth() < deltaDist) {
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
			Vector3d lu = Face::axisU(ptPos2) + initTimeIntv[0]*Face::axisU(ptVel2);//u延展的方向
			Vector3d lvtmp = Face::axisV(ptPos2) + initTimeIntv[0]*Face::axisV(ptVel2);//u延展的方向
			lu[0]*=1.01;
			// lvtmp[0]*=1.01;
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
			// if (cur.calcL1Dist(CpPos1, CpVel1, CpPos2, CpVel2) < deltaDist) {
			if (cur.calcWidth() < deltaDist) {
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
					// if(colTime[1]-colTime[0]<1e-12)return colTime[0];
					heap.emplace(divUvB2, divTime1);
				}
				if (primitiveVFCheck(CpPos1, CpVel1, ptPos2, ptVel2, bb, divTime2)){
					// if(colTime[1]-colTime[0]<1e-12)return colTime[0];
					heap.emplace(divUvB2, divTime2);
				}
			}
		}
		return -1;
	}
};
