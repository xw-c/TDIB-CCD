# pragma once
#include "mathOps.h"
// modification: hi-intvR calculated symmetrically
// modification: individual OBB axes
// modification: tlower --> time interval
template<typename ParamObj1, typename ParamObj2, typename ParamBound1, typename ParamBound2>
class SolverTD{
public:
	static bool primitiveCheck(const ParamObj1 &CpPos1, const ParamObj1 &CpVel1, 
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
			ch1.clear(); ch2.clear();

			// for(const auto& l:lines1)std::cout<<"lines1:  "<<l.k<<" "<<l.b<<"\n";
			robustCH(lines1, ch1, true, Array2d(0,upperTime));
			// for(const auto& l:ch1)std::cout<<"ch1:  "<<l.k<<" "<<l.b<<"\n";
			// for(const auto& l:lines2)std::cout<<"lines2:  "<<l.k<<" "<<l.b<<"\n";
			robustCH(lines2, ch2, false, Array2d(0,upperTime));
			// for(const auto& l:ch2)std::cout<<"ch2:  "<<l.k<<" "<<l.b<<"\n";
			const auto intvT = robustHullIntersect(ch1, ch2, Array2d(0,upperTime));
			// 	std::cout<<timeIntv.transpose()<<"\n"<<intvT.transpose()<<"\n";
			// 	std::cin.get();
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
		// for(const auto&l:feasibleIntvs){
		// 	if(l[0]==l[1]){
		// 		std::cout<<"timeIntv:"<<timeIntv.transpose()<<"\n";
		// 		for(const auto&l1:feasibleIntvs)std::cout<<"intv:"<<l1.transpose()<<"\n";
		// 		std::cin.get();
		// 	}
		// }
		
		if (feasibleIntvs.size()==0) {
			//这意味着整段时间都有碰撞
			colTime = timeIntv;
			return true; 
		}
		double minT = 0, maxT = upperTime;
		std::sort(feasibleIntvs.begin(), feasibleIntvs.end(), 
			[](const Array2d& intv1, const Array2d& intv2){
				return (intv1(0)<intv2(0));
			});
		if(feasibleIntvs[0](0)<=0){
			minT = feasibleIntvs[0](1);
			for(int i=1;i<feasibleIntvs.size();i++)
				if(feasibleIntvs[i](0)<minT) //不能加等，因为无碰撞给的是开区间，如果有),(的情况加等号会把这个情况漏掉
					minT=std::max(minT, feasibleIntvs[i](1));
				else break;
		}
		if(minT > maxT){ colTime = Array2d(-1,-1); return false; }
		
		std::sort(feasibleIntvs.begin(), feasibleIntvs.end(), 
			[](const Array2d& intv1, const Array2d& intv2){
				return (intv1(1)>intv2(1));
			});
		if(feasibleIntvs[0](1)>=upperTime){
			maxT = feasibleIntvs[0](0);
			for(int i=1;i<feasibleIntvs.size();i++)
				if(feasibleIntvs[i](1)>maxT) //不能加等，因为无碰撞给的是开区间，如果有),(的情况加等号会把这个情况漏掉
					maxT=std::min(maxT, feasibleIntvs[i](0));
				else break;
		}
		if(minT > maxT){ colTime = Array2d(-1,-1); return false; }
		// if(minT==maxT){
		// 	for(const auto&l:feasibleIntvs)std::cout<<"intv:"<<l.transpose()<<"\n";
		// }
		colTime = Array2d(minT + timeIntv[0], maxT + timeIntv[0]); 
		return true;
	}
	static void robustCH(std::vector<Line>& lines, std::vector<Line>& ch, 
				const bool getMaxCH, const Array2d& tIntv) {
		if(!getMaxCH)std::reverse(lines.begin(),lines.end());
		lines.erase(std::unique(lines.begin(), lines.end()), lines.end()); // 去重
		// std::cout<<lines.size()<<"\n";
		ch.clear();
		ch.push_back(lines[0]);
		int alpha = 1;
		while(alpha < lines.size()){
			// std::cout<<id<<"  "<<pts.size()<<"\n";
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
				if((getMaxCH&&chStart>0)||(!getMaxCH&&chStart<0))
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
		// TBD: intersect point \in [0,deltaT];
		if(ch.empty()){
			std::cout<<"empty CH!\n";
			exit(-1);
		}
	}
	static Array2d robustHullIntersect(const std::vector<Line>& ch1, const std::vector<Line>& ch2, 
							const Array2d& tIntv) {
		int id1=0, id2=0;
		double intvL=-1, intvR=-1;
		if(ch1[0].k*tIntv[0]+ch1[0].b<ch2[0].k*tIntv[0]+ch2[0].b)intvL=tIntv[0];
		else{
		//寻找intersection左端点
			while(id1<ch1.size()&&id2<ch2.size()){
				if(ch1[id1].k>=ch2[id2].k){
					break;
				}
				double hifp1, hifp2;
				if(id1<ch1.size()-1)
					hifp1=-((ch1[id1].k-ch1[id1+1].k)*(ch1[id1].b-ch2[id2].b)
							-(ch1[id1].k-ch2[id2].k)*(ch1[id1].b-ch1[id1+1].b));
				else 
					hifp1=tIntv[1]*(ch1[id1].k-ch2[id2].k)+(ch1[id1].b-ch2[id2].b);
				if(id2<ch2.size()-1)
					hifp2=(ch2[id2].k-ch2[id2+1].k)*(ch1[id1].b-ch2[id2].b)
							-(ch1[id1].k-ch2[id2].k)*(ch2[id2].b-ch2[id2+1].b);
				else
					hifp2=tIntv[1]*(ch1[id1].k-ch2[id2].k)+(ch1[id1].b-ch2[id2].b);
				// std::cout<<"finding intvL "<<id1<<" "<<id2<<", hifp1="<<hifp1<<"  hifp2="<<hifp2<<"\n";
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
			// std::cout<<"left: "<<id1<<"  "<<id2<<"  "<<intvL<<'\n';
			// if(intvL==-1)return Array2d(-1,-1);
			if(intvL==-1||intvL>=tIntv[1])return Array2d(-1,-1);
		}

		id1 = ch1.size()-1, id2 = ch2.size()-1;
		if((ch1[id1].k-ch2[id2].k)*tIntv[1]+(ch1[id1].b-ch2[id2].b)<0)intvR=tIntv[1];
		else{
		//寻找intersection左端点
			while(id1>=0&&id2>=0){
				if(ch1[id1].k<=ch2[id2].k){
					std::cerr<<"end at strange slopes?\n";
					exit(-1);
				}
				double hifp1, hifp2;
				if(id1>0)
					hifp1=(ch1[id1].k-ch1[id1-1].k)*(ch1[id1].b-ch2[id2].b)
							-(ch1[id1].k-ch2[id2].k)*(ch1[id1].b-ch1[id1-1].b);
				else 
					hifp1=tIntv[0]*(ch1[id1].k-ch2[id2].k)+(ch1[id1].b-ch2[id2].b);
				if(id2>0)
					hifp2=-((ch2[id2].k-ch2[id2-1].k)*(ch1[id1].b-ch2[id2].b)
							-(ch1[id1].k-ch2[id2].k)*(ch2[id2].b-ch2[id2-1].b));
				else
					hifp2=tIntv[0]*(ch1[id1].k-ch2[id2].k)+(ch1[id1].b-ch2[id2].b);
				// std::cout<<"finding intvL "<<id1<<" "<<id2<<", hifp1="<<hifp1<<"  hifp2="<<hifp2<<"\n";
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
			// std::cout<<"left: "<<id1<<"  "<<id2<<"  "<<intvL<<'\n';
			if(intvR==-1){
				std::cerr<<"intvL done but no intvR?\n";
				exit(-1);
			}
			if(intvR<=intvL)
				return Array2d(-1,-1);
		}

		
		if(intvL>intvR||intvL<tIntv[0]||intvR>tIntv[1]){
			std::cout<<"error intersection!\n";
			std::cout<<intvL<<" "<<intvR<<", in range"<<tIntv[0]<<" "<<tIntv[0]<<"\n";
			exit(-1);
		}
		else return Array2d(intvL,intvR);
	}			
public:
	static double solveCCD(const ParamObj1 &CpPos1, const ParamObj1 &CpVel1, 
						const ParamObj2 &CpPos2, const ParamObj2 &CpVel2,
						Array2d& uv1, Array2d& uv2, 
						const double upperTime = DeltaT,
						const double deltaDist = MinL1Dist) {
		struct PatchPair{
			ParamBound1 pb1;
			ParamBound2 pb2;
			Array2d tIntv;
			PatchPair(const ParamBound1& c1, const ParamBound2& c2, 
					const Array2d& t = Array2d(0,DeltaT)): pb1(c1), pb2(c2), tIntv(t) {}
			bool operator<(PatchPair const &o) const { return tIntv[0] > o.tIntv[0]; }
			double calcL1Dist(const ParamObj1 &CpPos1, const ParamObj1 &CpVel1, 
							const ParamObj2 &CpPos2, const ParamObj2 &CpVel2) const{
				auto ptPos1 = CpPos1.divideBezierPatch(pb1);
				auto ptVel1 = CpVel1.divideBezierPatch(pb1);
				auto ptPos2 = CpPos2.divideBezierPatch(pb2);
				auto ptVel2 = CpVel2.divideBezierPatch(pb2);
				for(int i=0;i<ParamObj1::cntCp;i++)
					ptPos1[i]+=ptVel1[i]*tIntv[0];
				for(int i=0;i<ParamObj2::cntCp;i++)
					ptPos2[i]+=ptVel2[i]*tIntv[0];
				double d1=calcAAExtent<ParamObj1>(ptPos1);
				double d2=calcAAExtent<ParamObj2>(ptPos2);
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
		cnt=1;
		while (!heap.empty()) {
			auto const cur = heap.top();
			heap.pop();
			cnt++;
			// if(SHOWANS) std::cout<<cnt<<"\n";

			// Decide whether the algorithm converges
			if (cur.calcL1Dist(CpPos1, CpVel1, CpPos2, CpVel2) < deltaDist) {
				uv1 = cur.pb1.centerParam();
				uv2 = cur.pb2.centerParam();
				const auto endTime = steady_clock::now();
				if(SHOWANS)
					std::cout << "min time: "<<  cur.tIntv[0] 
						<< "\nused seconds: " << duration(endTime - initialTime).count()
						<< std::endl;
				return cur.tIntv[0];
			}

			// Divide the current patch into four-to-four pieces
			for (int i = 0; i < 4; i++) {
				ParamBound1 divUvB1(cur.pb1.interpSubpatchParam(i));
				for (int j = 0; j < 4; j++) {
					ParamBound2 divUvB2(cur.pb2.interpSubpatchParam(j));
					if (primitiveCheck(CpPos1, CpVel1, CpPos2, CpVel2, divUvB1, divUvB2, colTime, cur.tIntv)){
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
