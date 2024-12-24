# pragma once
#include "linearGeom.h"

class LinearSolverTD{
private:
	static Vector3d clampAxis(const Vector3d& axis){
		Vector3d clp = axis;
		const double maxcoeff=axis.cwiseAbs().maxCoeff();
		for(int i=0;i<3;i++)if(std::abs(clp[i])<=maxcoeff*1e-10)clp[i]=0;
		return clp;
	}
	static Array2d axisCheck(std::vector<Line> lines1, std::vector<Line> lines2,
						const Array2d& timeIntv){
		std::vector<Line> ch1, ch2;
		ch1.clear(); ch2.clear();

		for( auto& l:lines1)l.b+=1e-12*abs(l.b);
		for( auto& l:lines2)l.b-=1e-12*abs(l.b);
		calcBoundaries(lines1, ch1, true, timeIntv);
		calcBoundaries(lines2, ch2, false, timeIntv);
		return boundaryIntersect(ch1, ch2, timeIntv);
	};
	static void calcBoundaries(std::vector<Line>& lines, std::vector<Line>& ch, 
				const bool getMaxCH, const Array2d& tIntv) {
		if(!getMaxCH)std::reverse(lines.begin(),lines.end());
		lines.erase(std::unique(lines.begin(), lines.end()), lines.end());
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
	static Array2d boundaryIntersect(const std::vector<Line>& ch1, const std::vector<Line>& ch2, 
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
					return Array2d(-1,-1);
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
			return Array2d(-1,-1);
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
		if(feasibleIntvs[0](0)<initTimeIntv[0]){
			minT = std::max(minT,feasibleIntvs[0](1));
			for(int i=1;i<feasibleIntvs.size();i++)
				if(feasibleIntvs[i](0)<minT)
					minT=std::max(minT, feasibleIntvs[i](1));
				else break;
		}
		if(minT > initTimeIntv[1]) { colTime = Array2d(-1,-1); return false; }
		
		std::sort(feasibleIntvs.begin(), feasibleIntvs.end(), 
			[](const Array2d& intv1, const Array2d& intv2){
				return (intv1(1)>intv2(1));
			});
		if(feasibleIntvs[0](1)>initTimeIntv[1]){
			maxT = std::min(maxT, feasibleIntvs[0](0));
			for(int i=1;i<feasibleIntvs.size();i++)
				if(feasibleIntvs[i](1)>maxT)
					maxT=std::min(maxT, feasibleIntvs[i](0));
				else break;
		}
		if(initTimeIntv[0] > maxT) { colTime = Array2d(-1,-1); return false; }
		if(minT >= maxT) { colTime = Array2d(maxT, minT);}
		else colTime = Array2d(minT, maxT); 
		return true;
	}
public:
	static bool primitiveEECheck(const std::array<Vector3d,Edge::cntCp> &ptPos1, 
						const std::array<Vector3d,Edge::cntCp> &ptVel1, 
						const std::array<Vector3d,Edge::cntCp> &ptPos2, 
						const std::array<Vector3d,Edge::cntCp> &ptVel2,
						Array2d& colTime,
						const BoundingBoxType& bb,
						const Array2d& initTimeIntv) {
		Array2d timeIntv(initTimeIntv[0]-1e-10, initTimeIntv[1]+1e-10);
		std::vector<Vector3d> axes;
		if(bb==BoundingBoxType::AABB){
			axes = {Vector3d::Unit(0), Vector3d::Unit(1), Vector3d::Unit(2)};
		}
		else if(bb==BoundingBoxType::OBB){
			Vector3d lu = Edge::direction(ptPos1) + initTimeIntv[0]*Edge::direction(ptVel1);
			Vector3d lvtmp = Edge::direction(ptPos2) + initTimeIntv[0]*Edge::direction(ptVel2);
			Vector3d ln = lu.cross(lvtmp);
			Vector3d lv = ln.cross(lu);
			axes = {lu, lv, ln};
		}

		std::vector<Array2d> feasibleIntvs;
		feasibleIntvs.clear();
		
		for( auto& axis:axes){
			std::vector<Line> ptLines1, ptLines2;
			ptLines1.clear(); ptLines2.clear();
			for(int i = 0; i < Edge::cntCp; i++) ptLines1.emplace_back(ptVel1[i].dot(axis), ptPos1[i].dot(axis));
			for(int i = 0; i < Edge::cntCp; i++) ptLines2.emplace_back(ptVel2[i].dot(axis), ptPos2[i].dot(axis));
			
			std::sort(ptLines1.begin(), ptLines1.end());
			std::sort(ptLines2.begin(), ptLines2.end());
			auto intvT = axisCheck(ptLines1, ptLines2, timeIntv);
			if(intvT[0]!=-1)feasibleIntvs.push_back(intvT);
			intvT = axisCheck(ptLines2, ptLines1, timeIntv);
			if(intvT[0]!=-1)feasibleIntvs.push_back(intvT);
		}
		auto b=intvMerge(feasibleIntvs, colTime, initTimeIntv);
		return b;
	}

	static double solveEETest(const Edge &CpPos1, const Edge &CpVel1, 
						const Edge &CpPos2, const Edge &CpVel2,
						double& u1, double& u2, 
						const BoundingBoxType& bb,
						const double upperTime,
						const double deltaDist) {
		std::priority_queue<EEPair> heap;
		Array2d initParam(0,1), initTimeIntv(0,upperTime), colTime;
		if (primitiveEECheck(CpPos1.ctrlp, CpVel1.ctrlp, CpPos2.ctrlp, CpVel2.ctrlp, colTime, bb, initTimeIntv))
			heap.emplace(initParam, initParam, colTime);
		
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
			for (int i = 0; i < 2; i++) {
				Array2d divUvB1[2]={Array2d(cur.pb1[0],mid1), Array2d(mid1, cur.pb1[1])};
				for (int j = 0; j < 2; j++) {
					Array2d divUvB2[2]={Array2d(cur.pb2[0],mid2), Array2d(mid2, cur.pb2[1])};
					auto ptPos1 = CpPos1.divideBezierPatch(divUvB1[i]);
					auto ptVel1 = CpVel1.divideBezierPatch(divUvB1[i]);
					auto ptPos2 = CpPos2.divideBezierPatch(divUvB2[j]);
					auto ptVel2 = CpVel2.divideBezierPatch(divUvB2[j]);
					if (primitiveEECheck(ptPos1, ptVel1, ptPos2, ptVel2, colTime, bb, cur.tIntv)){
						heap.emplace(divUvB1[i], divUvB2[j], colTime);
					}
				}
			}
		}
		return -1;
	}

	static bool primitiveVFCheck(const Vector3d &ptPos1, const Vector3d &ptVel1, 
						const std::array<Vector3d,Face::cntCp> &ptPos2, 
						const std::array<Vector3d,Face::cntCp> &ptVel2,
						Array2d& colTime,
						const BoundingBoxType& bb,
						const Array2d& initTimeIntv) {
		Array2d timeIntv(initTimeIntv[0]-1e-10, initTimeIntv[1]+1e-10);
		std::vector<Vector3d> axes;
		if(bb==BoundingBoxType::AABB){
			axes = {Vector3d::Unit(0), Vector3d::Unit(1), Vector3d::Unit(2)};
		}
		else if(bb==BoundingBoxType::OBB){
			Vector3d lu = Face::axisU(ptPos2) + initTimeIntv[0]*Face::axisU(ptVel2);
			Vector3d lvtmp = Face::axisV(ptPos2) + initTimeIntv[0]*Face::axisV(ptVel2);
			Vector3d ln = lu.cross(lvtmp);
			Vector3d lv = ln.cross(lu);
			axes = {lu, lv, ln};
		}

		std::vector<Array2d> feasibleIntvs;
		feasibleIntvs.clear();
		
		for(const auto& axis:axes){
			std::vector<Line> ptLines1, ptLines2;
			ptLines1.clear(); ptLines2.clear();
			ptLines1.emplace_back(ptVel1.dot(axis), ptPos1.dot(axis));
			for(int i = 0; i < Face::cntCp; i++) ptLines2.emplace_back(ptVel2[i].dot(axis), ptPos2[i].dot(axis));
			std::sort(ptLines1.begin(), ptLines1.end());
			std::sort(ptLines2.begin(), ptLines2.end());
			auto intvT = axisCheck(ptLines1, ptLines2, timeIntv);
			if(intvT[0]!=-1)feasibleIntvs.push_back(intvT);
			intvT = axisCheck(ptLines2, ptLines1, timeIntv);
			if(intvT[0]!=-1)feasibleIntvs.push_back(intvT);
		}
		return intvMerge(feasibleIntvs, colTime, initTimeIntv);
	}

	static double solveVFTest(const Vector3d &CpPos1, const Vector3d &CpVel1, 
						const Face &CpPos2, const Face &CpVel2,
						Array2d& uv,
						const BoundingBoxType& bb,
						const double upperTime,
						const double deltaDist) {
		std::priority_queue<VFPair> heap;
		TriParamBound initParam;
		Array2d initTimeIntv(0,upperTime), colTime;
		if (primitiveVFCheck(CpPos1, CpVel1, CpPos2.ctrlp, CpVel2.ctrlp, colTime, bb, initTimeIntv))
			heap.emplace(initParam, colTime);
		
		while (!heap.empty()) {
			auto const cur = heap.top();
			heap.pop();

			if (cur.calc4dWidth() < deltaDist) {
				uv = cur.pb.centerParam();
				return cur.tIntv[0];
			}

			// Divide the current patch into four-to-four pieces
			for (int j = 0; j < 4; j++) {
				TriParamBound divUvB2(cur.pb.interpSubpatchParam(j));
				auto ptPos2 = CpPos2.divideBezierPatch(divUvB2);
				auto ptVel2 = CpVel2.divideBezierPatch(divUvB2);
				if (primitiveVFCheck(CpPos1, CpVel1, ptPos2, ptVel2, colTime, bb, cur.tIntv)){
					heap.emplace(divUvB2, colTime);
				}
			}
		}
		return -1;
	}
};
