# pragma once
#include "utils.h"
template<typename ParamObj1, typename ParamObj2, typename ParamBound1, typename ParamBound2>
class SolverTD{
public:
	struct Line
	{
		double k, b;
		Line():k(0),b(0){}
		void set(const double& _k, const double &_b){k=_k,b=_b;}
		Line(const double& k,const double& b): k(k), b(b) {}
		bool operator<(const Line &l) const { 
			return k < l.k || (k == l.k && b > l.b); // 相同斜率的直线中只有截距最大的被留下来
		}
		bool operator==(const Line &l) const {return k == l.k;}
	};
	
	static bool primitiveCheck(const ParamObj1 &CpPos1, const ParamObj1 &CpVel1, 
						const ParamObj2 &CpPos2, const ParamObj2 &CpVel2,
						const ParamBound1 &divUvB1, const ParamBound2 &divUvB2,
						Array2d& colTime,
						const BoundingBoxType& bb,
						const Array2d& initTimeIntv = Array2d(0,DeltaT)) {
		auto ptPos1 = CpPos1.divideBezierPatch(divUvB1);
		auto ptVel1 = CpVel1.divideBezierPatch(divUvB1);
		auto ptPos2 = CpPos2.divideBezierPatch(divUvB2);
		auto ptVel2 = CpVel2.divideBezierPatch(divUvB2);
		// Enlarge the time interval by a small margin so that the end points are correctly treated 
		Array2d timeIntv(initTimeIntv[0]-1e-6, initTimeIntv[1]+1e-6);

		std::vector<Vector3d> axes;
		setAxes<ParamObj1, ParamObj2>(ptPos1, ptVel1, ptPos2, ptVel2, axes, bb, initTimeIntv[0]);

		std::vector<Array2d> feasibleIntvs;
		feasibleIntvs.clear();
		
		auto AxisCheck=[&](std::vector<Line> lines1, std::vector<Line> lines2){
			std::vector<Line> ch1, ch2;
			ch1.clear(); ch2.clear();
			for( auto& l:lines1)l.b += 1e-6 * abs(l.b);
			for( auto& l:lines2)l.b -= 1e-6 * abs(l.b);
			calcBoundaries(lines1, ch1, true, timeIntv);
			calcBoundaries(lines2, ch2, false, timeIntv);
			const auto intvT = boundaryIntersect(ch1, ch2, timeIntv);
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
		
		if (feasibleIntvs.size()==0) {
			// Meaning that collision is possible in the entire time interval
			colTime = initTimeIntv;
			return true; 
		}
		double minT = initTimeIntv[0], maxT = initTimeIntv[1];
		// Find the lower bound
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
		if(minT > maxT){ colTime = Array2d(-1,-1); return false; }
		
		// Find the upper bound
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

		if(minT >= maxT) { colTime = Array2d(maxT, minT);}
		else colTime = Array2d(minT, maxT); 
		return true;
	}

	static void calcBoundaries(std::vector<Line>& lines, std::vector<Line>& ch, 
				const bool getMaxCH, const Array2d& tIntv) {
		// True for calculating the max func, false for calculating the min func
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
			std::cerr<<"WARNING: you got an empty set of boundary! This output is not as expected and may be caused by incorrect usage or floating-point errors. It is recommended to pause the program for inspection.\n";
			exit(-1);
		}
	}

	static Array2d boundaryIntersect(const std::vector<Line>& ch1, const std::vector<Line>& ch2, 
							const Array2d& tIntv) {
		int id1=0, id2=0;

		// Find the left end of the intersection of the boundaries
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

		// Find the right end of the intersection of the boundaries
		id1 = ch1.size()-1, id2 = ch2.size()-1;
		if((ch1[id1].k-ch2[id2].k)*tIntv[1]+(ch1[id1].b-ch2[id2].b)<0)intvR=tIntv[1];
		else{
			while(id1>=0&&id2>=0){
				if(ch1[id1].k<=ch2[id2].k){
					// std::cerr<<"WARNING: calculation of boundary intersection ends at strange slopes! This output is not as expected and may be caused by incorrect usage or floating-point errors. It is recommended to pause the program for inspection.\n";
					// exit(-1);
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
				std::cerr<<"WARNING: find the left end of boundary intersection but no right end! This output is not as expected and may be caused by incorrect usage or floating-point errors. It is recommended to pause the program for inspection.\n";
				std::cerr<<"("<<intvL<<", "<<intvR<<"), in candidate ("<<tIntv[0]<<", "<<tIntv[1]<<")\n";
				exit(-1);
			}
			if(intvR<=intvL)
				return Array2d(-1,-1);
		}

		if(intvL>intvR||intvL<tIntv[0]||intvR>tIntv[1]){
			// std::cerr<<"WARNING: you got an in-valid sub-time interval! This output is not as expected and may be caused by incorrect usage or floating-point errors. It is recommended to pause the program for inspection.\n";
			// std::cerr<<"("<<intvL<<", "<<intvR<<"), in candidate ("<<tIntv[0]<<", "<<tIntv[1]<<")\n";
			// exit(-1);
			return Array2d(-1,-1);
		}
		else return Array2d(intvL,intvR);
	}			
public:
	static double solveCCD(const ParamObj1 &CpPos1, const ParamObj1 &CpVel1, 
						const ParamObj2 &CpPos2, const ParamObj2 &CpVel2,
						Array2d& uv1, Array2d& uv2, 
						const BoundingBoxType& bb,
						const double deltaDist,
						const double upperTime = DeltaT) {
		struct PatchPair{
			ParamBound1 pb1;
			ParamBound2 pb2;
			Array2d tIntv;
			PatchPair(const ParamBound1& c1, const ParamBound2& c2, 
					const Array2d& t = Array2d(0,DeltaT)): pb1(c1), pb2(c2), tIntv(t) {}
			bool operator<(PatchPair const &o) const { return tIntv[0] > o.tIntv[0]; }
			double calcWidth() const{
				const double w1 = pb1.width(), w2 = pb2.width();
				return std::max(std::max(w1, w2), tIntv[1]-tIntv[0]);
			}
		};

		std::priority_queue<PatchPair> heap;
		ParamBound1 initParam1;
		ParamBound2 initParam2;
		Array2d initTimeIntv(0,upperTime), colTime;
		if (primitiveCheck(CpPos1, CpVel1, CpPos2, CpVel2, initParam1, initParam2, colTime, bb, initTimeIntv))
			heap.emplace(initParam1, initParam2, colTime);
		while (!heap.empty()) {
			auto const cur = heap.top();
			heap.pop();

			// Meets the precision requirement
			if (cur.calcWidth() < deltaDist) {
				uv1 = cur.pb1.centerParam();
				uv2 = cur.pb2.centerParam();
				return cur.tIntv[0];
			}

			// Divide the current patch into four-to-four pieces
			for (int i = 0; i < 4; i++) {
				ParamBound1 divUvB1(cur.pb1.interpSubpatchParam(i));
				for (int j = 0; j < 4; j++) {
					ParamBound2 divUvB2(cur.pb2.interpSubpatchParam(j));
					if (primitiveCheck(CpPos1, CpVel1, CpPos2, CpVel2, divUvB1, divUvB2, colTime, bb, cur.tIntv)){
						heap.emplace(divUvB1, divUvB2, colTime);
					}
				}
			}
		}

		return -1;
	}
};
