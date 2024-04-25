# pragma once
#include "mathOps.h"
#define DEBUG 0
template<typename ParamObj1, typename ParamObj2, typename ParamBound1, typename ParamBound2>
class SolverRobustTD{
public:
	static bool primitiveCheck(const ParamObj1 &CpPos1, const ParamObj1 &CpVel1, 
						const ParamObj2 &CpPos2, const ParamObj2 &CpVel2,
						const ParamBound1 &divUvB1, const ParamBound2 &divUvB2,
						Array2d& colTime, Array2d& colTimeError,
						const BoundingBoxType& bb,
						const Array2d& initTimeIntv = Array2d(0,DeltaT), 
						const Array2d& timeError = Array2d(0,0)) {
		auto ptPos1 = CpPos1.divideBezierPatch(divUvB1);
		auto ptVel1 = CpVel1.divideBezierPatch(divUvB1);
		auto ptPos2 = CpPos2.divideBezierPatch(divUvB2);
		auto ptVel2 = CpVel2.divideBezierPatch(divUvB2);
		Array2d timeIntv(initTimeIntv[0]-1e-6, initTimeIntv[1]+1e-6);

		std::vector<Vector3d> axes;
		double boundA;
		auto setRobustAxes =[&](const double &t){	
			if(bb==BoundingBoxType::AABB){
				axes = {Vector3d::Unit(0), Vector3d::Unit(1), Vector3d::Unit(2)};
				boundA = 1;
			}
			else if(bb==BoundingBoxType::OBB){
				Vector3d lu1 = ParamObj1::axisU(ptPos1) + t*ParamObj1::axisU(ptVel1);//u延展的方向
				Vector3d lv1tmp = ParamObj1::axisV(ptPos1) + t*ParamObj1::axisV(ptVel1);//v延展的方向
				while(lu1.cwiseAbs().maxCoeff()>=1)lu1*=0.1;
				while(lv1tmp.cwiseAbs().maxCoeff()>=1)lv1tmp*=0.1;
				double maxCoeff1 = (lu1.cwiseAbs().cwiseMax(lv1tmp.cwiseAbs())).maxCoeff();
				if(maxCoeff1==0){
					std::cerr<<"degenerate bounding box!\n";
					exit(-1);
				}
				Vector3d ln1 = lu1.cross(lv1tmp);
				Vector3d lv1 = ln1.cross(lu1);

				Vector3d lu2 = ParamObj2::axisU(ptPos2) + t*ParamObj2::axisU(ptVel2);//u延展的方向
				Vector3d lv2tmp = ParamObj2::axisV(ptPos2) + t*ParamObj2::axisU(ptVel2);//v延展的方向
				while(lu2.cwiseAbs().maxCoeff()>=1)lu2*=0.1;
				while(lv2tmp.cwiseAbs().maxCoeff()>=1)lv2tmp*=0.1;
				double maxCoeff2 = (lu2.cwiseAbs().cwiseMax(lv2tmp.cwiseAbs())).maxCoeff();
				if(maxCoeff2==0){
					std::cerr<<"degenerate bounding box!\n";
					exit(-1);
				}
				Vector3d ln2 = lu2.cross(lv2tmp);
				Vector3d lv2 = ln2.cross(lu2);

				axes = {lu1,lv1,ln1,lu2,lv2,ln2, 
					lu1.cross(lu2), lu1.cross(lv2), lu1.cross(ln2), 
					lv1.cross(lu2), lv1.cross(lv2), lv1.cross(ln2), 
					ln1.cross(lu2), ln1.cross(lv2), ln1.cross(ln2)};

				boundA = std::max(maxCoeff1, maxCoeff2);
			}
		};
		setRobustAxes(initTimeIntv[0]);
		// 可能要加kt的地方
		int orderCH1 = 7 * ParamObj1::order + 22;
		int orderCH2 = 7 * ParamObj2::order + 22;
		double boundMHCoeff1 = (6*std::pow(4,ParamObj1::order+2))*std::pow(boundA, 6);
		double boundMHCoeff2 = (6*std::pow(4,ParamObj2::order+2))*std::pow(boundA, 6);
		double boundP1 = 0, boundP2 = 0;
		for(int i=0;i<ParamObj1::cntCp;i++){
			boundP1 = std::max(boundP1, (CpPos1.ctrlp[i].cwiseAbs()+CpVel1.ctrlp[i].cwiseAbs()*timeIntv[1]).maxCoeff());
		}
		for(int i=0;i<ParamObj2::cntCp;i++){
			boundP2 = std::max(boundP2, (CpPos2.ctrlp[i].cwiseAbs()+CpVel2.ctrlp[i].cwiseAbs()*timeIntv[1]).maxCoeff());
		}
		double boundCHProj1 = boundMHCoeff1 * boundP1;
		double boundCHProj2 = boundMHCoeff2 * boundP2;
		double errorCHProj1 = boundCHProj1 * MachineEps * orderCH1;
		double errorCHProj2 = boundCHProj2 * MachineEps * orderCH2;

		std::vector<Array2dError> feasibleIntvs;
		feasibleIntvs.clear();

		auto AxisCheck=[&](std::vector<Line> lines1, std::vector<Line> lines2){
			std::vector<Line > ch1, ch2;
			ch1.clear(); ch2.clear();

			// if(DEBUG)std::cout<<"min coeffs:"<<lines1[10].k<<" "<<lines1[10].b<<"\n";
			// if(DEBUG)std::cout<<"min coeffs:"<<lines2[10].k<<" "<<lines2[10].b<<"\n";
			// for(const auto& l:lines1)std::cout<<"lines1:  "<<l.k<<" "<<l.b<<"\n";
			double errorCH1 = robustCH(lines1, ch1, true, timeIntv);
			// for(const auto& l:ch1)std::cout<<"ch1:  "<<l.k<<" "<<l.b<<"\n";
			// for(const auto& l:lines2)std::cout<<"lines2:  "<<l.k<<" "<<l.b<<"\n";
			double errorCH2 = robustCH(lines2, ch2, false, timeIntv);
			// for(const auto& l:ch2)std::cout<<"ch2:  "<<l.k<<" "<<l.b<<"\n";
			// if(DEBUG)std::cout<<"error parts of 1: "<<errorCHProj1<<" "<<errorCH1<<"\n";
			// if(DEBUG)std::cout<<"error parts of 2: "<<errorCHProj2<<" "<<errorCH2<<"\n";
			errorCH1 += errorCHProj1, errorCH2 += errorCHProj2;
			// std::cout<<"CH errors: "<<errorCH1<<"  "<<errorCH2<<"\n";

			// 更简单一些的errorCH版本
			double maxbias=lines1[0].b, minbias=lines1[0].b;
			for(auto l:lines1){maxbias=std::max(maxbias,l.b);minbias=std::min(minbias,l.b);}
			double boundBias1 = boundCHProj1 + lines1.size() * 
								std::max(4*(maxbias-minbias), 
										(lines1.back().k-lines1.front().k)*timeIntv[1]+(maxbias-minbias));
			maxbias=lines2[0].b, minbias=lines2[0].b;
			for(auto l:lines2){maxbias=std::max(maxbias,l.b);minbias=std::min(minbias,l.b);}
			double boundBias2 = boundCHProj2 + lines2.size() * 
								std::max(4*(maxbias-minbias), 
										(lines2.back().k-lines2.front().k)*timeIntv[1]+(maxbias-minbias));
			// std::cout<<"error parts of 2: "<<boundCHProj2<<" "<<boundBias2<<"\n";
			int orderBias1 = orderCH1 + 2, orderBias2 = orderCH2 + 2;
			// for(const auto& l:ch1)std::cout<<"ch1:  "<<l.k<<" "<<l.b<<"\n";
			// for(const auto& l:ch2)std::cout<<"ch2:  "<<l.k<<" "<<l.b<<"\n";
			const auto intvT = robustHullIntersect(ch1, ch2, boundBias1, boundBias2, orderBias1, orderBias2, timeIntv);
			// if(DEBUG)std::cout<<intvT.first.transpose()<<"\n"<<intvT.second.transpose()<<"\n";
			// 	std::cin.get();
			if(intvT.first[0]!=-1)feasibleIntvs.push_back(intvT);
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
				// std::cin.get();
		// 	}
		// }
		if (feasibleIntvs.size()==0) {
			//这意味着整段时间都有碰撞
			colTime = initTimeIntv;
			return true; 
		}
		double minT = initTimeIntv[0], maxT = initTimeIntv[1];
		double minTerror, maxTerror;
		std::sort(feasibleIntvs.begin(), feasibleIntvs.end(), 
			[](const Array2dError& intv1, const Array2dError& intv2){
				return (intv1.first(0)-intv2.first(0)<intv2.second(0)-intv1.second(0));
				// return (intv1(0)<intv2(0));
			});
		if(DEBUG)
			for(const auto&l1:feasibleIntvs)
				std::cout<<"intv:"<<l1.first.transpose()<<" with error "<<l1.second.transpose()<<"\n";
		if(feasibleIntvs[0].first(0)-initTimeIntv[0]<=-feasibleIntvs[0].second(0)){
		// if(feasibleIntvs[0](0)<=initTimeIntv[0]){
			minT = feasibleIntvs[0].first(1);
			minTerror = feasibleIntvs[0].second(1);
			for(int i=1;i<feasibleIntvs.size();i++)
				if(feasibleIntvs[i].first(0)-minT<-minTerror-feasibleIntvs[i].second(0)){ //不能加等，因为无碰撞给的是开区间，如果有),(的情况加等号会把这个情况漏掉
					if(minT-feasibleIntvs[i].first(1)<minTerror-feasibleIntvs[i].second(1))
						minT=feasibleIntvs[i].first(1), minTerror = feasibleIntvs[i].second(1);
				}
				else break;
		}
		if(DEBUG)std::cout<<minT<<", "<<minTerror<<"\n";
		if(minT > maxT){ colTime = Array2d(-1,-1); return false; }
		
		std::sort(feasibleIntvs.begin(), feasibleIntvs.end(), 
			[](const Array2dError& intv1, const Array2dError& intv2){
				return (intv1.first(1)-intv2.first(1)>intv1.second(1)-intv2.second(1));
			});
		if(feasibleIntvs[0].first(1)-initTimeIntv[1]>=feasibleIntvs[0].second(1)){
			maxT = feasibleIntvs[0].first(0),
			maxTerror = feasibleIntvs[0].second(0);
			for(int i=1;i<feasibleIntvs.size();i++)
				if(feasibleIntvs[i].first(1)-maxT>feasibleIntvs[i].second(1)+maxTerror) {//不能加等，因为无碰撞给的是开区间，如果有),(的情况加等号会把这个情况漏掉
					if(maxT-feasibleIntvs[i].first(0)<-maxTerror+feasibleIntvs[i].second(0))
						maxT=feasibleIntvs[i].first(0),maxTerror=feasibleIntvs[i].second(0);
				}
				else break;
		}
		if(DEBUG)std::cout<<maxT<<", "<<maxTerror<<"\n";
		if(minT > maxT){ colTime = Array2d(-1,-1); return false; }
		// if(DEBUG)std::cin.get();
		// if(minT==maxT){
		// 	for(const auto&l:feasibleIntvs)std::cout<<"intv:"<<l.transpose()<<"\n";
		// }
		colTime = Array2d(minT+minTerror, maxT-maxTerror); 
		return true;
	}
	static double robustCH(std::vector<Line>& lines, std::vector<Line>& ch, 
				const bool getMaxCH, const Array2d& tIntv) {
		if(!getMaxCH)std::reverse(lines.begin(),lines.end());
		lines.erase(std::unique(lines.begin(), lines.end()), lines.end()); // 去重
		// std::cout<<lines.size()<<"\n";
		std::vector<double> errorLines;

		ch.clear();
		errorLines.clear();
		ch.push_back(lines[0]);
		errorLines.push_back(0);
		int alpha = 1;
		while(alpha < lines.size()){
			// std::cout<<id<<"  "<<pts.size()<<"\n";
			int beta = ch.size()-1;
			double errorLine = 0;
			while(beta > 0){
				double chfp = (ch[beta].k-ch[beta-1].k)*(lines[alpha].b-ch[beta-1].b)
							-(lines[alpha].k-ch[beta-1].k)*(ch[beta].b-ch[beta-1].b);
				double boundCH = std::abs(lines[alpha].k-ch[beta-1].k)*(std::abs(lines[alpha].b-ch[beta-1].b)+std::abs(ch[beta].b-ch[beta-1].b));
				double errorCH = boundCH * MachineEps * 9;
				//尽可能pop出去
				if(chfp>=-errorCH){
					errorLine += errorLines.back();
					ch.pop_back();
					errorLines.pop_back();
					beta--;
					if(chfp<errorCH){
						errorLine += (std::abs(lines[alpha].b-ch[beta-1].b)+std::abs(ch[beta].b-ch[beta-1].b)) * MachineEps * 13;
					}
				}
				else break;
			}
			if(beta==0){
				double chStart = tIntv[0]*(lines[alpha].k-ch[0].k)+(lines[alpha].b-ch[0].b);
				double boundCHStart = tIntv[0]*std::abs(lines[alpha].k-ch[0].k)+std::abs(lines[alpha].b-ch[0].b);
				double errorCHStart = boundCHStart * MachineEps * 6;
				if(!getMaxCH) chStart = - chStart;
				if(chStart>=-errorCHStart){
					ch.pop_back();
					if(chStart<errorCHStart){
						double errorLineStart = errorCHStart;
						errorLine += errorLineStart;
					}
				}
			}
			if(ch.empty() || errorLine!=0){ch.push_back(lines[alpha]); errorLines.push_back(errorLine);}
			else{
				double chEnd = tIntv[1]*(lines[alpha].k-ch[beta].k)+(lines[alpha].b-ch[beta].b);
				double boundCHEnd = tIntv[1]*std::abs(lines[alpha].k-ch[beta].k)+std::abs(lines[alpha].b-ch[beta].b);
				double errorCHEnd = boundCHEnd * MachineEps * 6;
				if(!getMaxCH) chEnd = - chEnd;
				if(chEnd>errorCHEnd){
					ch.push_back(lines[alpha]);
					errorLines.push_back(errorLine);
				}
				else if(chEnd>-errorCHEnd){
					double errorLineEnd = errorCHEnd;
					errorLines.back() = std::max(errorLines.back(), errorLineEnd);
				}
			}
			alpha++;
		}

		if(ch.empty()){
			std::cout<<"empty CH!\n";
			exit(-1);
		}
		return *std::max_element(errorLines.begin(), errorLines.end());
	}
	static Array2dError robustHullIntersect(const std::vector<Line>& ch1, const std::vector<Line>& ch2, 
							const double& boundBias1, const double& boundBias2, 
							const int& orderBias1, const int& orderBias2, 
							const Array2d& tIntv) {
		int id1=0, id2=0;
		double intvL=-1, intvR=-1;

		int orderBias = std::max(orderBias1, orderBias2);
		double errorIntvL = 0, errorIntvR = 0;

		double boundStart = std::abs((ch1[0].k-ch2[0].k)*tIntv[0])+(boundBias1+boundBias2);
		double errorStart = boundStart*MachineEps*(orderBias+2);
		if((ch1[0].k-ch2[0].k)*tIntv[0]+(ch1[0].b-ch2[0].b)<-errorStart)intvL=tIntv[0];//尽可能不要
		// if((ch1[0].k-ch2[0].k)*tIntv[0]+(ch1[0].b-ch2[0].b)<0)intvL=tIntv[0];
		else{
			//寻找intersection左端点
			while(id1<ch1.size()&&id2<ch2.size()){
				if(ch1[id1].k-ch2[id2].k>=-SqrtMachineEps){//尽可能删，可以漏
					break;
				}
				double hifp1, hifp2, errorHIfp1, errorHIfp2;
				if(id1<ch1.size()-1)
					// hifp1=-((ch1[id1].k-ch1[id1+1].k)*(ch1[id1].b-ch2[id2].b)
					// 		-(ch1[id1].k-ch2[id2].k)*(ch1[id1].b-ch1[id1+1].b));
					hifp1=(ch1[id1+1].k-ch2[id2].k)*(ch1[id1].b-ch2[id2].b)
							-(ch1[id1].k-ch2[id2].k)*(ch1[id1+1].b-ch2[id2].b),
					errorHIfp1=(std::abs(ch1[id1+1].k-ch2[id2].k)+std::abs(ch1[id1].k-ch2[id2].k))*(boundBias1+boundBias2)
								*MachineEps*(orderBias+9);
				else 
					hifp1=tIntv[1]*(ch1[id1].k-ch2[id2].k)+(ch1[id1].b-ch2[id2].b),
					errorHIfp1=(std::abs(tIntv[1]*(ch1[id1].k-ch2[id2].k))+(boundBias1+boundBias2))
								*MachineEps*(orderBias+2);
				if(id2<ch2.size()-1)
					// hifp2=(ch2[id2].k-ch2[id2+1].k)*(ch1[id1].b-ch2[id2].b)
					// 		-(ch1[id1].k-ch2[id2].k)*(ch2[id2].b-ch2[id2+1].b);
					hifp2=(ch1[id1].k-ch2[id2+1].k)*(ch1[id1].b-ch2[id2].b)
							-(ch1[id1].k-ch2[id2].k)*(ch1[id1].b-ch2[id2+1].b),
					errorHIfp2=(std::abs(ch1[id1].k-ch2[id2+1].k)+std::abs(ch1[id1].k-ch2[id2].k))*(boundBias1+boundBias2)
								*MachineEps*(orderBias+9);
				else
					hifp2=tIntv[1]*(ch1[id1].k-ch2[id2].k)+(ch1[id1].b-ch2[id2].b),
					errorHIfp2=(std::abs(tIntv[1]*(ch1[id1].k-ch2[id2].k))+(boundBias1+boundBias2))
								*MachineEps*(orderBias+2);
				// std::cout<<"finding intvL "<<id1<<" "<<id2<<", hifp1="<<hifp1<<"  hifp2="<<hifp2<<"\n";
				if(hifp1<-errorHIfp1){//尽可能往后找，可以漏不能多
					if(hifp2<-errorHIfp2){
						intvL = -(ch1[id1].b-ch2[id2].b)/(ch1[id1].k-ch2[id2].k);
						// std::cout<<"error of hi-L: "<<errorIntvL<<" ";
						errorIntvL += (boundBias1+boundBias2)/(ch2[id2].k-ch1[id1].k)*MachineEps*(orderBias+2);
						// std::cout<<errorIntvL<<"\n";
						// std::cout<<"id: "<<id1<<" "<<id2<<"\n";
						// std::cout<<"hi: "<<hifp1<<" "<<hifp2<<"\n";
						// std::cout<<"error: "<<errorHIfp1<<" "<<errorHIfp2<<"\n";
						// std::cout<<"res: "<<intvL<<" "<<errorIntvL<<"\n";
						break;
					}
					else {
						id2++;
						if(id2<ch2.size()&&hifp2<errorHIfp2)//ch2进入考虑下一条
							errorIntvL += 2*(boundBias1+boundBias2)*MachineEps*(orderBias+11);
						else if(id2<ch2.size()&&errorIntvL!=0){
							std::cerr<<"nonzero errorIntvL at id2?\n";
							for(const auto& l:ch1)std::cout<<"ch1:  "<<l.k<<" "<<l.b<<"\n";
							for(const auto& l:ch2)std::cout<<"ch2:  "<<l.k<<" "<<l.b<<"\n";
							std::cout<<"id: "<<id1<<" "<<id2<<"\n";
							std::cout<<"hi: "<<hifp1<<" "<<hifp2<<"\n";
							std::cout<<"error: "<<errorHIfp1<<" "<<errorHIfp2<<"\n";
							std::cout<<"res: "<<intvL<<" "<<errorIntvL<<"\n";
							exit(-1);
						}
					}
				}
				else {
					id1++;
					// if(hifp2<0);
					// else id2++;
					if(id1<ch1.size()&&hifp1<errorHIfp1)//ch2进入考虑下一条
						errorIntvL += 2*(boundBias1+boundBias2)*MachineEps*(orderBias+11);
					else if(id1<ch1.size()&&errorIntvL!=0){
						std::cerr<<"nonzero errorIntvL at id1?\n";
						for(const auto& l:ch1)std::cout<<"ch1:  "<<l.k<<" "<<l.b<<"\n";
						for(const auto& l:ch2)std::cout<<"ch2:  "<<l.k<<" "<<l.b<<"\n";
						std::cout<<"id: "<<id1<<" "<<id2<<"\n";
						std::cout<<"hi: "<<hifp1<<" "<<hifp2<<"\n";
						std::cout<<"error: "<<errorHIfp1<<" "<<errorHIfp2<<"\n";
						std::cout<<"res: "<<intvL<<" "<<errorIntvL<<"\n";
						exit(-1);
					}
				}
			}
			// std::cout<<"left: "<<id1<<"  "<<id2<<"  "<<intvL<<'\n';
			// if(intvL==-1)return Array2d(-1,-1);
			if(intvL==-1||intvL-tIntv[1]>=-errorIntvL)return Array2dError(Array2d(-1,-1),Array2d(-1,-1));
		}

		id1 = ch1.size()-1, id2 = ch2.size()-1;
		double boundEnd = std::abs((ch1[id1].k-ch2[id2].k)*tIntv[1])+(boundBias1+boundBias2);
		double errorEnd = boundStart*MachineEps*(orderBias+2);
		if((ch1[id1].k-ch2[id2].k)*tIntv[1]+(ch1[id1].b-ch2[id2].b)<-errorEnd)intvR=tIntv[1];
		else{
			//寻找intersection右端点
			while(id1>=0&&id2>=0){
				if(ch1[id1].k-ch2[id2].k<=0){//SqrtMachineEps){
					std::cerr<<"end at strange slopes?\n";
					std::cerr<<intvL<<" "<<id1<<" "<<id2<<"\n";
					for(const auto& l:ch1)std::cout<<"ch1:  "<<l.k<<" "<<l.b<<"\n";
					for(const auto& l:ch2)std::cout<<"ch2:  "<<l.k<<" "<<l.b<<"\n";
					exit(-1);
				}
				double hifp1, hifp2, errorHIfp1, errorHIfp2;
				if(id1>0)
					// hifp1=(ch1[id1].k-ch1[id1-1].k)*(ch1[id1].b-ch2[id2].b)
					// 		-(ch1[id1].k-ch2[id2].k)*(ch1[id1].b-ch1[id1-1].b);
					hifp1=(ch1[id1].k-ch2[id2].k)*(ch1[id1-1].b-ch2[id2].b)
							-(ch1[id1-1].k-ch2[id2].k)*(ch1[id1].b-ch2[id2].b),
					errorHIfp1=(std::abs(ch1[id1-1].k-ch2[id2].k)+std::abs(ch1[id1].k-ch2[id2].k))*(boundBias1+boundBias2)
							*MachineEps*(orderBias+9);
				else 
					hifp1=tIntv[0]*(ch1[id1].k-ch2[id2].k)+(ch1[id1].b-ch2[id2].b),
					errorHIfp1=(std::abs(tIntv[0]*(ch1[id1].k-ch2[id2].k))+(boundBias1+boundBias2))
								*MachineEps*(orderBias+2);
				if(id2>0)
					// hifp2=-((ch2[id2].k-ch2[id2-1].k)*(ch1[id1].b-ch2[id2].b)
					// 		-(ch1[id1].k-ch2[id2].k)*(ch2[id2].b-ch2[id2-1].b));
					hifp2=(ch1[id1].k-ch2[id2].k)*(ch1[id1].b-ch2[id2-1].b)
							-(ch1[id1].k-ch2[id2-1].k)*(ch1[id1].b-ch2[id2].b),
					errorHIfp2=(std::abs(ch1[id1].k-ch2[id2-1].k)+std::abs(ch1[id1].k-ch2[id2].k))*(boundBias1+boundBias2)
								*MachineEps*(orderBias+9);
				else
					hifp2=tIntv[0]*(ch1[id1].k-ch2[id2].k)+(ch1[id1].b-ch2[id2].b),
					errorHIfp2=(std::abs(tIntv[0]*(ch1[id1].k-ch2[id2].k))+(boundBias1+boundBias2))
								*MachineEps*(orderBias+2);
				// std::cout<<"finding intvL "<<id1<<" "<<id2<<", hifp1="<<hifp1<<"  hifp2="<<hifp2<<"\n";
				if(hifp1<-errorHIfp1){
					if(hifp2<-errorHIfp2){
						intvR = -(ch1[id1].b-ch2[id2].b)/(ch1[id1].k-ch2[id2].k);
						errorIntvR += (boundBias1+boundBias2)/(ch1[id1].k-ch2[id2].k)*MachineEps*(orderBias+2);
						break;
					}
					else {
						id2--;
						if(id2>=0&&hifp2<errorHIfp2)//ch2进入考虑下一条
							errorIntvR += 2*(boundBias1+boundBias2)*MachineEps*(orderBias+11);
						else if(id2>=0&&errorIntvR!=0){
							std::cerr<<"nonzero errorIntvR at id2?\n";
							for(const auto& l:ch1)std::cout<<"ch1:  "<<l.k<<" "<<l.b<<"\n";
							for(const auto& l:ch2)std::cout<<"ch2:  "<<l.k<<" "<<l.b<<"\n";
							std::cout<<"id: "<<id1<<" "<<id2<<"\n";
							std::cout<<"hi: "<<hifp1<<" "<<hifp2<<"\n";
							std::cout<<"error: "<<errorHIfp1<<" "<<errorHIfp2<<"\n";
							std::cout<<"res: "<<intvL<<" "<<errorIntvL<<"\n";
							std::cout<<"res: "<<intvR<<" "<<errorIntvR<<"\n";
							exit(-1);
						}
					}
				}
				else {
					id1--;
					// if(hifp2<0);
					// else id2--;
					if(id1>=0&&hifp1<errorHIfp1)//ch2进入考虑下一条
						errorIntvR += 2*(boundBias1+boundBias2)*MachineEps*(orderBias+11);
					else if(id1>=0&&errorIntvR!=0){
						std::cerr<<"nonzero errorIntvR at id1?\n";
						for(const auto& l:ch1)std::cout<<"ch1:  "<<l.k<<" "<<l.b<<"\n";
						for(const auto& l:ch2)std::cout<<"ch2:  "<<l.k<<" "<<l.b<<"\n";
						std::cout<<"id: "<<id1<<" "<<id2<<"\n";
						std::cout<<"hi: "<<hifp1<<" "<<hifp2<<"\n";
						std::cout<<"error: "<<errorHIfp1<<" "<<errorHIfp2<<"\n";
						std::cout<<"res: "<<intvL<<" "<<errorIntvL<<"\n";
						std::cout<<"res: "<<intvR<<" "<<errorIntvR<<"\n";
						exit(-1);
					}
				}
			}
			// std::cout<<"left: "<<id1<<"  "<<id2<<"  "<<intvL<<'\n';
			if(intvR==-1){
				std::cerr<<"intvL done but no intvR?\n";
				exit(-1);
			}
			if(intvR-intvL<=errorIntvR+errorIntvL)
				return Array2dError(Array2d(-1,-1),Array2d(-1,-1));
		}

		
		if(intvL>intvR||intvL<tIntv[0]||intvR>tIntv[1]){
			std::cout<<"error intersection!\n";
			std::cout<<intvL<<" "<<intvR<<", in range"<<tIntv[0]<<" "<<tIntv[0]<<"\n";
			exit(-1);
		}
		else return Array2dError(Array2d(intvL,intvR),Array2d(errorIntvL, errorIntvR));
	}	
public:
	static double solveCCD(const ParamObj1 &CpPos1, const ParamObj1 &CpVel1, 
						const ParamObj2 &CpPos2, const ParamObj2 &CpVel2,
						Array2d& uv1, Array2d& uv2, 
						const BoundingBoxType& bb,
						const double upperTime = DeltaT,
						const double deltaDist = MinL1Dist) {
		struct PatchPair{
			ParamBound1 pb1;
			ParamBound2 pb2;
			Array2d tIntv, tError;
			PatchPair(const ParamBound1& c1, const ParamBound2& c2, 
					const Array2d& t = Array2d(0,DeltaT),
					const Array2d& te = Array2d(0,0)): pb1(c1), pb2(c2), tIntv(t), tError(te) {}
			bool operator<(PatchPair const &o) const { return tIntv[0] - o.tIntv[0] > tError[0] - o.tError[0]; }
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
		Array2d initTimeIntv(0,upperTime), initTimeError(0,0), colTime, colTimeError;
		if (primitiveCheck(CpPos1, CpVel1, CpPos2, CpVel2, initParam1, initParam2, colTime, colTimeError, bb, initTimeIntv, initTimeError))
			heap.emplace(initParam1, initParam2, colTime, colTimeError);
		if(DEBUG){
			std::cout<<colTime.transpose()<<"\n"<<colTimeError.transpose()<<"\n";
			std::cin.get();
		}
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
				return cur.tIntv[0]; // 减掉tError应该也看不出来？
			}

			// Divide the current patch into four-to-four pieces
			for (int i = 0; i < 4; i++) {
				ParamBound1 divUvB1(cur.pb1.interpSubpatchParam(i));
				for (int j = 0; j < 4; j++) {
					ParamBound2 divUvB2(cur.pb2.interpSubpatchParam(j));
					if (primitiveCheck(CpPos1, CpVel1, CpPos2, CpVel2, divUvB1, divUvB2, colTime, colTimeError, bb, initTimeIntv, initTimeError)){
						heap.emplace(divUvB1, divUvB2, colTime, colTimeError);
						if(DEBUG){
							std::cout<<colTime.transpose()<<"\n"<<colTimeError.transpose()<<"\n";
							std::cin.get();
						}
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
