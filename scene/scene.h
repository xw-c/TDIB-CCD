#pragma once
#include"paramMesh.h"
#include"solverBase.h"
#include"solverBaseManifold.h"
#include"solverTD.h"
#include"solverTDManifold.h"
#include"utilOps.h"
#include "triBezier.h"
#include "triRatBezier.h"
#include "recBezier.h"
#include "recRatBezier.h"
// template<typename ObjType, typename ParamType>
// void randomTest(const double denom){
// 	// std::ofstream fp("posDiff.txt");
// 	// std::ofstream ft("timeDiff.txt");
// 	ObjType obj1, obj2, vel1, vel2;
// 	std::srand(0);
// 	int hasCol = 0;
// 	double t;
// 	Array2d uv1,uv2;
// 	using steady_clock = std::chrono::steady_clock;
// 	using duration = std::chrono::duration<double>;
// 	const auto initialTime = steady_clock::now();
// 	for(int kase = 0;kase<Kase;kase++){
// 		generatePatchPair<ObjType>(obj1.ctrlp, vel1.ctrlp, obj2.ctrlp, vel2.ctrlp, denom);
// 		switch(solverType){
// 			case SolverType::TDIntv:
// 				t = SolverTD<ObjType,ObjType,ParamType,ParamType>::solveCCD(obj1,vel1,obj2,vel2,uv1,uv2,DeltaT);
// 				break;
// 			case SolverType::BaseIntv:
// 				t = SolverBase<ObjType,ObjType,ParamType,ParamType>::solveCCD(obj1,vel1,obj2,vel2,uv1,uv2,DeltaT);
// 				break;
// 			default:
// 				std::cerr<<"solver not implemented!\n";
// 				exit(-1);
// 		}
// 		if(t>=0)hasCol++;
// 		// if(kase==21){
// 		// 	saveDoFs<ObjType>(obj1.ctrlp, vel1.ctrlp, obj2.ctrlp, vel2.ctrlp);
// 		// 	break;
// 		// }
// 		// ft<<t<<"\n";
// 		if(SHOWANS) std::cout<<calcDist(obj1,vel1,obj2,vel2,uv1,uv2,t)<<"\n";
// 		// std::cout<<kase<<": "<<duration(steady_clock::now() - initialTime).count()<<"s\n";
// 	}
// 	const auto endTime = steady_clock::now();
// 	std::cout << hasCol<<" pairs have collided.\n";
// 	std::cout << "used seconds: " <<
// 		duration(endTime - initialTime).count()/Kase
// 		<< std::endl;
// 	// fp.close();
// 	// ft.close();
// }

template<typename ObjType, typename ParamType>
void randomTest(const std::string& solverType, const double& deltaDist, const int& kase, const double& velMag, const std::string& outputFile){
	ObjType obj1, obj2, vel1, vel2;
	std::srand(0);
	int hasCol = 0;
	double t;
	Array2d uv1,uv2;
	std::ofstream file("validation/"+outputFile+".txt");
	file << std::fixed << std::setprecision(10);

	using steady_clock = std::chrono::steady_clock;
	using duration = std::chrono::duration<double>;
	const auto initialTime = steady_clock::now();
	for(int k = 0; k < kase; k ++){
		generatePatchPair<ObjType>(obj1.ctrlp, vel1.ctrlp, obj2.ctrlp, vel2.ctrlp, velMag);
		if(solverType=="td")
			t = SolverTD<ObjType,ObjType,ParamType,ParamType>::solveCCD(obj1,vel1,obj2,vel2,uv1,uv2,DeltaT,deltaDist);
		else if(solverType=="base")
			t = SolverBase<ObjType,ObjType,ParamType,ParamType>::solveCCD(obj1,vel1,obj2,vel2,uv1,uv2,DeltaT,deltaDist);
		else{
			std::cerr<<"solver not implemented!\n";
			exit(-1);
		}
		if(SHOWANS) std::cout<<cnt<<"\n";
		if(t>=0)hasCol++;
		// if(k==3){
		// 	saveDoFs<ObjType,ObjType>(obj1.ctrlp, vel1.ctrlp, obj2.ctrlp, vel2.ctrlp);
		// 	exit(-1);
		// }
		// ft<<t<<"\n";
		file<<t<<"\n";
		// if(t>0)file<<t<<"  "<<uv1[0]<<"  "<<uv1[1]<<"  "<<uv2[0]<<"  "<<uv2[1]<<"\n";
		if(SHOWANS) std::cout<<"case "<<k<<" done "<<calcDist(obj1,vel1,obj2,vel2,uv1,uv2,t)<<"\n";
		
		// std::cout<<kase<<": "<<duration(steady_clock::now() - initialTime).count()<<"s\n";
	}
	const auto endTime = steady_clock::now();
	std::cout << hasCol<<" pairs have collided.\n";
	std::cout << "used seconds: " <<
		duration(endTime - initialTime).count()/kase
		<< std::endl;
	file.close();
	// fp.close();
	// ft.close();
}


template<typename ObjType, typename ParamType>
void planeTest(const std::string& solverType, const double& deltaDist, const int& kase, const double& velMag, const std::string& outputFile){
	ObjType obj1, obj2, vel1, vel2;
	int hasCol = 0;
	double t;
	Array2d uv1,uv2;
	std::multiset<CCDIntv<ParamType, ParamType> > solutSet;
	solutSet.clear();

	std::ofstream file("validation/"+outputFile+".txt");
	file << std::fixed << std::setprecision(10);
	
	obj1.ctrlp = {
		Vector3d(0,0,0), Vector3d(0,1,0), Vector3d(0,2,0), Vector3d(0,3,0),
		Vector3d(1,0,0), Vector3d(1,1,0), Vector3d(1,2,0), Vector3d(1,3,0),
		Vector3d(2,0,0), Vector3d(2,1,0), Vector3d(2,2,0), Vector3d(2,3,0),
		Vector3d(3,0,0), Vector3d(3,1,0), Vector3d(3,2,0), Vector3d(3,3,0)
	}, 
	vel1.ctrlp = {Vector3d(0,0,1), Vector3d(0,0,1), Vector3d(0,0,1), Vector3d(0,0,1), 
		Vector3d(0,0,1), Vector3d(0,0,1), Vector3d(0,0,1), Vector3d(0,0,1),
		Vector3d(0,0,1), Vector3d(0,0,1), Vector3d(0,0,1), Vector3d(0,0,1),
		Vector3d(0,0,1), Vector3d(0,0,1), Vector3d(0,0,1), Vector3d(0,0,1)};
	obj2.ctrlp = {
		Vector3d(0,0,1.1), Vector3d(0,1,1.1), Vector3d(0,2,1.1), Vector3d(0,3,1.1),
		Vector3d(1,0,1.1), Vector3d(1,1,1.1), Vector3d(1,2,1.1), Vector3d(1,3,1.1),
		Vector3d(2,0,1.1), Vector3d(2,1,1.1), Vector3d(2,2,1.1), Vector3d(2,3,1.1),
		Vector3d(3,0,1.1), Vector3d(3,1,1.1), Vector3d(3,2,1.1), Vector3d(3,3,1.1)
	}, 
	vel2.ctrlp = {Vector3d(0,0,-1), Vector3d(0,0,-1), Vector3d(0,0,-1), Vector3d(0,0,-1),
		Vector3d(0,0,-1), Vector3d(0,0,-1), Vector3d(0,0,-1), Vector3d(0,0,-1),
		Vector3d(0,0,-1), Vector3d(0,0,-1), Vector3d(0,0,-1), Vector3d(0,0,-1),
		Vector3d(0,0,-1), Vector3d(0,0,-1), Vector3d(0,0,-1), Vector3d(0,0,-1)};
	using steady_clock = std::chrono::steady_clock;
	using duration = std::chrono::duration<double>;
	const auto initialTime = steady_clock::now();
	if(solverType=="td"){
		SolverTDManifold<ObjType,ObjType,ParamType,ParamType> solver;
		t = solver.solveCCD(obj1,vel1,obj2,vel2,solutSet,DeltaT,deltaDist);
	}else if(solverType=="base"){
		SolverBaseManifold<ObjType,ObjType,ParamType,ParamType> solver;
		t = solver.solveCCD(obj1,vel1,obj2,vel2,solutSet,DeltaT,deltaDist);}
	else{
		std::cerr<<"solver not implemented!\n";
		exit(-1);
	}
	if(SHOWANS) std::cout<<" done "<<calcDist(obj1,vel1,obj2,vel2,uv1,uv2,t)<<"\n";
	Vector3d const p1 = obj1.evaluatePatchPoint(uv1);
	Vector3d const v1 = vel1.evaluatePatchPoint(uv1);
	Vector3d const p2 = obj2.evaluatePatchPoint(uv2);
	Vector3d const v2 = vel2.evaluatePatchPoint(uv2);
	Vector3d const pt1=(v1*t+p1), pt2=(v2*t+p2);
	std::cout<<"uv at:"<<uv1.transpose()<<"     "<<uv2.transpose()<<"\npos at: "<<p1.transpose()<<"     "<<p2.transpose()<<"\n";
	std::cout<<"delta: "<<(pt2-pt1).norm()<<"\n";
	const auto endTime = steady_clock::now();
	std::cout << "used seconds: " <<
		duration(endTime - initialTime).count()
		<< std::endl;
	
	for(const auto& r:solutSet){
		auto v=(r.aabb1[0]+r.aabb1[1])/2.;
		file<<r.t<<" "<<v[0]<<" "<<v[1]<<"\n";
	}
	file.close();
	// fp.close();
	// ft.close();
}


template<typename ObjType, typename ParamType>
void validate(const std::string& solverType, const double& deltaDist, const int& kase, const double& velMag, const std::string& outputFile){
	ObjType obj1, obj2, vel1, vel2;
	int hasCol = 0;
	double t;
	Array2d uv1,uv2;
	std::multiset<CCDIntv<ParamType, ParamType> > solutSet;
	solutSet.clear();
	
	obj1.ctrlp = {
		Vector3d(0,0,0), Vector3d(0,1,1), Vector3d(0,2,2), Vector3d(0,3,3),
		Vector3d(1,0,1), Vector3d(1,1,2), Vector3d(1,2,3), Vector3d(1,3,4),
		Vector3d(2,0,2), Vector3d(2,1,3), Vector3d(2,2,4), Vector3d(2,3,5),
		Vector3d(3,0,3), Vector3d(3,1,4), Vector3d(3,2,5), Vector3d(3,3,6)
	}, 
	vel1.ctrlp = {Vector3d(0,0,1), Vector3d(0,0,1), Vector3d(0,0,1), Vector3d(0,0,1), 
		Vector3d(0,0,1), Vector3d(0,0,1), Vector3d(0,0,1), Vector3d(0,0,1),
		Vector3d(0,0,1), Vector3d(0,0,1), Vector3d(0,0,1), Vector3d(0,0,1),
		Vector3d(0,0,1), Vector3d(0,0,1), Vector3d(0,0,1), Vector3d(0,0,1)};
	obj2.ctrlp = {
		Vector3d(0,0,7.1), Vector3d(0,1,7.1), Vector3d(0,2,7.1), Vector3d(0,3,7.1),
		Vector3d(1,0,7.1), Vector3d(1,1,7.1), Vector3d(1,2,7.1), Vector3d(1,3,7.1),
		Vector3d(2,0,7.1), Vector3d(2,1,7.1), Vector3d(2,2,7.1), Vector3d(2,3,7.1),
		Vector3d(3,0,7.1), Vector3d(3,1,7.1), Vector3d(3,2,7.1), Vector3d(3,3,7.1)
	}, 
	vel2.ctrlp = {Vector3d(0,0,-1), Vector3d(0,0,-1), Vector3d(0,0,-1), Vector3d(0,0,-1),
		Vector3d(0,0,-1), Vector3d(0,0,-1), Vector3d(0,0,-1), Vector3d(0,0,-1),
		Vector3d(0,0,-1), Vector3d(0,0,-1), Vector3d(0,0,-1), Vector3d(0,0,-1),
		Vector3d(0,0,-1), Vector3d(0,0,-1), Vector3d(0,0,-1), Vector3d(0,0,-1)};

	readinDoFs<ObjType,ObjType>(obj1.ctrlp, vel1.ctrlp, obj2.ctrlp, vel2.ctrlp);
	using steady_clock = std::chrono::steady_clock;
	using duration = std::chrono::duration<double>;
	const auto initialTime = steady_clock::now();
	if(solverType=="td")
		t = SolverTD<ObjType,ObjType,ParamType,ParamType>::solveCCD(obj1,vel1,obj2,vel2,uv1,uv2,DeltaT,deltaDist);
	else if(solverType=="base")
		t = SolverBase<ObjType,ObjType,ParamType,ParamType>::solveCCD(obj1,vel1,obj2,vel2,uv1,uv2,DeltaT,deltaDist);
	else{
		std::cerr<<"solver not implemented!\n";
		exit(-1);
	}
	if(SHOWANS) std::cout<<" done "<<calcDist(obj1,vel1,obj2,vel2,uv1,uv2,t)<<"\n";
	// Vector3d const p1 = obj1.evaluatePatchPoint(uv1);
	// Vector3d const v1 = vel1.evaluatePatchPoint(uv1);
	// Vector3d const p2 = obj2.evaluatePatchPoint(uv2);
	// Vector3d const v2 = vel2.evaluatePatchPoint(uv2);
	// Vector3d const pt1=(v1*t+p1), pt2=(v2*t+p2);
	// std::cout<<"uv at:"<<uv1.transpose()<<"     "<<uv2.transpose()<<"\npos at: "<<p1.transpose()<<"     "<<p2.transpose()<<"\n";
	// std::cout<<"delta: "<<(pt2-pt1).norm()<<"\n";
	const auto endTime = steady_clock::now();
	std::cout << "used seconds: " <<
		duration(endTime - initialTime).count()
		<< std::endl;
	// fp.close();
	// ft.close();
}


// template<typename ObjType, typename ParamType>
// void planeTest(const std::string& solverType, const double& deltaDist, const int& kase, const double& velMag, const std::string& outputFile){
// 	ObjType obj1, obj2, vel1, vel2;
// 	int hasCol = 0;
// 	double t;
// 	Array2d uv1,uv2;
// 	obj1.ctrlp = {
// 		Vector3d(0,0,0), Vector3d(0,1,0), Vector3d(0,2,0), Vector3d(0,3,0),
// 		Vector3d(1,0,0), Vector3d(1,1,0), Vector3d(1,2,0), Vector3d(1,3,0),
// 		Vector3d(2,0,0), Vector3d(2,1,0), Vector3d(2,2,0), Vector3d(2,3,0),
// 		Vector3d(3,0,0), Vector3d(3,1,0), Vector3d(3,2,0), Vector3d(3,3,0)
// 	}, 
// 	vel1.ctrlp = {Vector3d(0,0,1), Vector3d(0,0,1), Vector3d(0,0,1), Vector3d(0,0,1), 
// 		Vector3d(0,0,1), Vector3d(0,0,1), Vector3d(0,0,1), Vector3d(0,0,1),
// 		Vector3d(0,0,1), Vector3d(0,0,1), Vector3d(0,0,1), Vector3d(0,0,1),
// 		Vector3d(0,0,1), Vector3d(0,0,1), Vector3d(0,0,1), Vector3d(0,0,1)};
// 	obj2.ctrlp = {
// 		Vector3d(0,0,1.1), Vector3d(0,1,1.1), Vector3d(0,2,1.1), Vector3d(0,3,1.1),
// 		Vector3d(1,0,1.1), Vector3d(1,1,1.1), Vector3d(1,2,1.1), Vector3d(1,3,1.1),
// 		Vector3d(2,0,1.1), Vector3d(2,1,1.1), Vector3d(2,2,1.1), Vector3d(2,3,1.1),
// 		Vector3d(3,0,1.1), Vector3d(3,1,1.1), Vector3d(3,2,1.1), Vector3d(3,3,1.1)
// 	}, 
// 	vel2.ctrlp = {Vector3d(0,0,-1), Vector3d(0,0,-1), Vector3d(0,0,-1), Vector3d(0,0,-1),
// 		Vector3d(0,0,-1), Vector3d(0,0,-1), Vector3d(0,0,-1), Vector3d(0,0,-1),
// 		Vector3d(0,0,-1), Vector3d(0,0,-1), Vector3d(0,0,-1), Vector3d(0,0,-1),
// 		Vector3d(0,0,-1), Vector3d(0,0,-1), Vector3d(0,0,-1), Vector3d(0,0,-1)};
// 	using steady_clock = std::chrono::steady_clock;
// 	using duration = std::chrono::duration<double>;
// 	const auto initialTime = steady_clock::now();
// 	if(solverType=="td")
// 		t = SolverTD<ObjType,ObjType,ParamType,ParamType>::solveCCD(obj1,vel1,obj2,vel2,uv1,uv2,DeltaT,deltaDist);
// 	else if(solverType=="base")
// 		t = SolverBase<ObjType,ObjType,ParamType,ParamType>::solveCCD(obj1,vel1,obj2,vel2,uv1,uv2,DeltaT,deltaDist);
// 	else{
// 		std::cerr<<"solver not implemented!\n";
// 		exit(-1);
// 	}
// 	if(SHOWANS) std::cout<<" done "<<calcDist(obj1,vel1,obj2,vel2,uv1,uv2,t)<<"\n";
// 	Vector3d const p1 = obj1.evaluatePatchPoint(uv1);
// 	Vector3d const v1 = vel1.evaluatePatchPoint(uv1);
// 	Vector3d const p2 = obj2.evaluatePatchPoint(uv2);
// 	Vector3d const v2 = vel2.evaluatePatchPoint(uv2);
// 	Vector3d const pt1=(v1*t+p1), pt2=(v2*t+p2);
// 	std::cout<<"uv at:"<<uv1.transpose()<<"     "<<uv2.transpose()<<"\npos at: "<<p1.transpose()<<"     "<<p2.transpose()<<"\n";
// 	std::cout<<"delta: "<<(pt2-pt1).norm()<<"\n";
// 	const auto endTime = steady_clock::now();
// 	std::cout << "used seconds: " <<
// 		duration(endTime - initialTime).count()
// 		<< std::endl;
// 	// fp.close();
// 	// ft.close();
// }

template<typename Mesh1, typename Mesh2, typename Func>
double ccd(const Mesh1& mesh1, const Mesh1& vel1,
			const Mesh2& mesh2, const Mesh2& vel2,
			Func&& paramCCD, 
			const double upperTime = DeltaT){
	double minTime = upperTime;
	Array2d uv1, uv2;
	if(mesh1.cntPatches!=vel1.cntPatches||mesh2.cntPatches!=vel2.cntPatches){
		std::cerr<<"dofs do not match!\n";
		exit(-1);
	}
	//问一下gpt：1、两个for auto patch同时遍历，2、为什么Func后面有两个&，3、为什么传参不能不写DeltaT
	for(int i = 0; i < mesh1.cntPatches; i++)
		for(int j = 0; j < mesh2.cntPatches; j++){
			double t = paramCCD(mesh1.patches[i], vel1.patches[i], mesh2.patches[j], vel2.patches[j], 
								uv1, uv2, minTime);
			std::cout<<i<<" "<<j<<": "<<t<<"\n";
			if(t>=0){
				if(t==0)return 0;
				minTime = std::min(minTime, t);
			}
		}
	return minTime;
}


template<typename ParamObj1, typename ParamObj2, typename ParamBound1, typename ParamBound2>
double ccd_ratBezier_manifold(const std::string& solverType, const std::string& outputFile,
			const RatParamMesh<ParamObj1>& mesh1, const RatParamMesh<ParamObj1>& vel1,
			const RatParamMesh<ParamObj2>& mesh2, const RatParamMesh<ParamObj2>& vel2,
			const double upperTime, const double& deltaDist){
	double minTimeUB = upperTime;
	std::multiset<CCDIntv<ParamBound1, ParamBound2> > solutSet;
	solutSet.clear();
	for(int i = 0; i < mesh1.cntPatches; i++)
		for(int j = 0; j < mesh2.cntPatches; j++){
			// if(i==10&&j==1)DEBUG=1;
			double t;
			if(solverType=="td"){
				SolverTDManifold<ParamObj1,ParamObj2,ParamBound1,ParamBound2> solver;
				t = solver.solveCCD(mesh1.patches[i], vel1.patches[i], mesh2.patches[j], vel2.patches[j], 
								solutSet, minTimeUB, deltaDist);
			}
			else if(solverType=="base"){
				SolverBaseManifold<ParamObj1,ParamObj2,ParamBound1,ParamBound2> solver;
				t = solver.solveCCD(mesh1.patches[i], vel1.patches[i], mesh2.patches[j], vel2.patches[j], 
								solutSet, minTimeUB, deltaDist);
			}
			else{
				std::cerr<<"solver not implemented!\n";
				exit(-1);
			}
			std::cout<<i<<" "<<j<<": "<<t<<solutSet.size()<<"\n";
			if(t>=0){
				if(t==0)return 0;
				minTimeUB = std::min(minTimeUB, t);
			}
		}
	std::ofstream output("manifold/"+outputFile+".txt");
	// output<<solutSet.size()<<" solutions\n";
	// for(const auto& r:solutSet){
	// 	output<<"time = "<<r.t<<"\npatch 1: "<<((r.aabb1[0]+r.aabb1[1])/2.).transpose()
	// 			<<"\npatch 2:"<<((r.aabb2[0]+r.aabb2[1])/2.).transpose()<<"\n";
	// }
	// output<<solutSet.size()<<"\n";
	for(const auto& r:solutSet){
		auto v=(r.aabb1[0]+r.aabb1[1])/2.;
		output<<r.tIntv.transpose()<<" "<<v[0]<<" "<<v[1]<<"\n";
		// output<<r.tIntv.transpose()<<": "<<v.transpose()<<", "<<((r.aabb2[0]+r.aabb2[1])/2.).transpose()<<"\n";
	}
	// auto primitiveMaxDist = [&](const CCDRoot& r1, const CCDRoot&r2){
	// 	Vector3d aaExtent1 = (r1.aabb1[1]-r2.aabb1[0]).cwiseMax(r2.aabb1[1]-r1.aabb1[0]),
	// 	aaExtent2 = (r1.aabb2[1]-r2.aabb2[0]).cwiseMax(r2.aabb2[1]-r1.aabb2[0]);
	// 	return std::max(aaExtent1.maxCoeff(), aaExtent2.maxCoeff());
	// };
	// int i=0;
	// for(const auto& r1:solutSet){
	// 	int j=0;
	// 	for(const auto& r2:solutSet)
	// 		output<<i<<", "<<j++<<"  "<<primitiveMaxDist(r1, r2)<<"\n";
	// 	i++;
	// }

	// for(const auto& r1:solutSet){
	// 	for(const auto& r2:solutSet){
	// 		Vector3d aaExtent1 = (r1.aabb1[1]-r2.aabb1[0]).cwiseMax(r2.aabb1[1]-r1.aabb1[0]),
	// 		aaExtent2 = (r1.aabb2[1]-r2.aabb2[0]).cwiseMax(r2.aabb2[1]-r1.aabb2[0]);
	// 		std::cout<<std::max(aaExtent1.maxCoeff(), aaExtent2.maxCoeff())<<"\n";
	// 	}
	// }
	output.close();
	return minTimeUB;
	// double t = solver.solveCCD(mesh1.patches[2], vel1.patches[2], mesh2.patches[0], vel2.patches[0], 
	// 				solutSet, minTime);
	// std::cout<<"\ngt: r="<<2-2./sqrt(5)<<"  z="<<4./sqrt(5)<<"\n\n";
	// std::cout<<solutSet.size()<<" solutions\n";
	// for(const auto& r:solutSet){
	// 	std::cout<<"time = "<<r.t<<"\npatch 1: "<<((r.aabb1[0]+r.aabb1[1])/2.).transpose()
	// 			<<"\npatch 2:"<<((r.aabb2[0]+r.aabb2[1])/2.).transpose()<<"\n";
	// }
	// return t;
}

// typedef RatParamMesh<RecQuadRatBezier> RatParamMesh<RecQuadRatBezier>;
// typedef RatParamMesh<TriQuadRatBezier> RatParamMesh<TriQuadRatBezier>;

void generateTorus(RatParamMesh<RecQuadRatBezier>& torus, RatParamMesh<RecQuadRatBezier>& vel){
	double r=2, a=1;
	//顶点位置没有乘权重
	// std::array<double, 9> weight = {1,1,2,1,1,2,2,2,4};
	std::array<double, 9> weight = {1,1./sqrt(2),1,1./sqrt(2),0.5,1./sqrt(2),1,1./sqrt(2),1};
	std::array<Vector3d, 9> pos1 = { Vector3d(r+a,0,0), Vector3d(r+a,0,a), Vector3d(r,0,a), 
				Vector3d(r+a,r+a,0), Vector3d(r+a,r+a,a), Vector3d(r,r,a), 
				Vector3d(0,r+a,0), Vector3d(0,r+a,a), Vector3d(0,r,a)	},

	pos2 = { Vector3d(r-a,0,0), Vector3d(r-a,0,a), Vector3d(r,0,a), 
				Vector3d(r-a,r-a,0), Vector3d(r-a,r-a,a), Vector3d(r,r,a), 
				Vector3d(0,r-a,0), Vector3d(0,r-a,a), Vector3d(0,r,a)	},

	pos3 = { Vector3d(r-a,0,0), Vector3d(r-a,0,-a), Vector3d(r,0,-a), 
				Vector3d(r-a,r-a,0), Vector3d(r-a,r-a,-a), Vector3d(r,r,-a), 
				Vector3d(0,r-a,0), Vector3d(0,r-a,-a), Vector3d(0,r,-a)	},

	pos4 = { Vector3d(r+a,0,0), Vector3d(r+a,0,-a), Vector3d(r,0,-a), 
				Vector3d(r+a,r+a,0), Vector3d(r+a,r+a,-a), Vector3d(r,r,-a), 
				Vector3d(0,r+a,0), Vector3d(0,r+a,-a), Vector3d(0,r,-a)	};

	torus.patches.clear();
	torus.patches.emplace_back(pos1,weight);
	torus.patches.emplace_back(pos2,weight);
	torus.patches.emplace_back(pos3,weight);
	torus.patches.emplace_back(pos4,weight);
	RatParamMesh<RecQuadRatBezier> torus2=torus;
	for(int i=0;i<3;i++){
		torus2.rotateObj(PI*0.5, Vector3d::Unit(2));
		torus.patches.insert(torus.patches.end(),torus2.patches.begin(),torus2.patches.end());
	}
	torus.cntPatches = 16;
	vel = RatParamMesh<RecQuadRatBezier>(16, weight);
}
void generateCone(RatParamMesh<TriQuadRatBezier>& cone, RatParamMesh<TriQuadRatBezier>& vel){
	double r=2, h=4;
	//顶点位置没有乘权重
	std::array<double, 6> weight = {1,1./sqrt(2),1, 1,1,1};
	std::array<Vector3d, 6> pos1 = { Vector3d(r,0,0), Vector3d(r,r,0), Vector3d(0,r,0), 
				Vector3d(r/2.,0,h/2.), Vector3d(0,r/2.,h/2.), Vector3d(0,0,h)},

	pos2 = { Vector3d(0,r,0), Vector3d(-r,r,0), Vector3d(-r,0,0), 
				Vector3d(0,r/2.,h/2.), Vector3d(-r/2.,0,h/2.), Vector3d(0,0,h)},

	pos3 = { Vector3d(-r,0,0), Vector3d(-r,-r,0), Vector3d(0,-r,0), 
				Vector3d(-r/2.,0,h/2.), Vector3d(0,-r/2.,h/2.), Vector3d(0,0,h)},

	pos4 = { Vector3d(0,-r,0), Vector3d(r,-r,0), Vector3d(r,0,0), 
				Vector3d(0,-r/2.,h/2.), Vector3d(r/2.,0,h/2.), Vector3d(0,0,h)};

	cone.patches.clear();
	cone.patches.emplace_back(pos1,weight);
	cone.patches.emplace_back(pos2,weight);
	cone.patches.emplace_back(pos3,weight);
	cone.patches.emplace_back(pos4,weight);

	cone.cntPatches = 4;

	vel = RatParamMesh<TriQuadRatBezier>(4, weight);
}
void torusTest(const std::string& solverType, const double& deltaDist, const std::string& outputFile){
	RatParamMesh<RecQuadRatBezier> torusPos, torusVel;
	generateTorus(torusPos,torusVel);

	RatParamMesh<RecQuadRatBezier> torusPos2=torusPos, torusVel2= torusVel;
	torusPos2.moveObj(Vector3d(0,0,6));
	torusVel2.moveObj(Vector3d(0,0,-5));// gt: 0.8s

	using steady_clock = std::chrono::steady_clock;
	using duration = std::chrono::duration<double>;
	const auto initialTime = steady_clock::now();
	// double t=ccd(torusPos, torusVel, conePos, coneVel, 
	// 			SolverBase<RecQuadRatBezier,TriQuadRatBezier,RecParamBound, TriParamBound>::solveCCD);
	double t = ccd_ratBezier_manifold<RecQuadRatBezier,RecQuadRatBezier,RecParamBound, RecParamBound>(solverType, outputFile, torusPos, torusVel, torusPos2, torusVel2, DeltaT, deltaDist);
	const auto endTime = steady_clock::now();
	std::cout<<"colTime: "<<t<<"\n";
	std::cout << "used seconds: " <<
		duration(endTime - initialTime).count()
		<< std::endl;
}
// 这个圆锥有点问题，这应该不是一个圆锥而是曲四棱锥
// void manifoldTest(){
// 	RatParamMesh<TriQuadRatBezier> conePos, coneVel;
// 	generateCone(conePos,coneVel);
// 	conePos.writeObj("cone.obj", 0.1,0.1);
// 	RatParamMesh<RecQuadRatBezier> torusPos, torusVel;
// 	generateTorus(torusPos,torusVel);
// 	torusPos.writeObj("torus.obj");

// 	// // gt: 0.5s
// 	// torusPos.moveObj(Vector3d(0,0,2.5+sqrt(5)));
// 	// torusVel.moveObj(Vector3d(0,0,-5));

// 	// gt: 0.2s
// 	torusPos.moveObj(Vector3d(2,0,6));
// 	torusVel.moveObj(Vector3d(0,0,-5));

// 	using steady_clock = std::chrono::steady_clock;
// 	using duration = std::chrono::duration<double>;
// 	const auto initialTime = steady_clock::now();
// 	// double t=ccd(torusPos, torusVel, conePos, coneVel, 
// 	// 			SolverBase<RecQuadRatBezier,TriQuadRatBezier,RecParamBound, TriParamBound>::solveCCD);
// 	double t=ccd_ratBezier_manifold<RecQuadRatBezier,TriQuadRatBezier,RecParamBound, TriParamBound>(torusPos, torusVel, conePos, coneVel);
// 	const auto endTime = steady_clock::now();
// 	std::cout<<"colTime: "<<t<<"\n";
// 	std::cout << "used seconds: " <<
// 		duration(endTime - initialTime).count()
// 		<< std::endl;
// }

// take down the case settings
/*
void testTeapot(){
	RecBezierMesh obj1("teapot32.txt"), obj2 = obj1;
	RecBezierMesh vel1(obj1.cntPatches), vel2 = vel1;
	// obj.writeObj("teapot32.obj");
	Vector3d displace(10,0,0);
	obj2.moveObj(displace);
	obj1.rotateObj(PI, Vector3d(0,1,0)+Vector3d::Random());
	vel1.moveObj(displace);
	
	using steady_clock = std::chrono::steady_clock;
	using duration = std::chrono::duration<double>;
	const auto initialTime = steady_clock::now();
	double t=ccd(obj1, vel1, obj2, vel2, recBezierCCD);
		const auto endTime = steady_clock::now();
	std::cout<<"colTime: "<<t<<"\n";
	std::cout << "used seconds: " <<
		duration(endTime - initialTime).count()
		<< std::endl;

	obj1.writeObj("teapot1.obj");
	obj2.writeObj("teapot2.obj");
	obj1.moveObj(displace * t);
	obj1.writeObj("teapotCol1.obj");
	obj2.writeObj("teapotCol2.obj");
	
	//
	// RecBezierMesh testobj1(1), testobj2(1);
	// testobj1.patches[0]=obj1.patches[28];
	// testobj2.patches[0]=obj2.patches[28];
	// testobj2.patches[1]=obj2.patches[30];
	// testobj2.patches[2]=obj2.patches[31];
	// testobj1.moveObj(Vector3d(0,0,5));
	// testobj2.moveObj(Vector3d(0,0,5));
	// testobj1.writeObj("test1.obj");
	// testobj2.writeObj("test2.obj");

}
void testTorus(){
	RatParamMesh<RecQuadRatBezier> obj1, vel1;
	generateTorusComponent(obj1, vel1);
	obj1.writeObj("torus.obj");

	RatParamMesh<RecQuadRatBezier> obj2 = obj1, vel2 = vel1;
	// obj.writeObj("teapot32.obj");
	Vector3d displace(6,6,0);
	obj2.moveObj(displace);
	obj1.rotateObj(PI, Vector3d(0,1,0)+Vector3d::Random());
	vel1.moveObj(displace);
	
	using steady_clock = std::chrono::steady_clock;
	using duration = std::chrono::duration<double>;
	const auto initialTime = steady_clock::now();
	double t=ccd(obj1, vel1, obj2, vel2, recRatBezierCCD);
		const auto endTime = steady_clock::now();
	std::cout<<"colTime: "<<t<<"\n";
	std::cout << "used seconds: " <<
		duration(endTime - initialTime).count()
		<< std::endl;

	obj1.moveObj(displace*t);
	obj1.writeObj("torusCol1.obj");
	obj2.writeObj("torusCol2.obj");
}

void testBunny(){
	TriLinearMesh obj1("bunny292.obj"), vel1(obj1.cntPatches);
	obj1.writeObj("bunny.obj");

	TriLinearMesh obj2 = obj1, vel2 = vel1;
	// obj.writeObj("teapot32.obj");
	Vector3d displace(6,6,0);
	obj2.moveObj(displace*0.5);
	obj1.rotateObj(PI, Vector3d(0,1,0)+Vector3d::Random());
	vel1.moveObj(displace);
	std::cout<<"preprocess done!\n";
	using steady_clock = std::chrono::steady_clock;
	using duration = std::chrono::duration<double>;
	const auto initialTime = steady_clock::now();
	double t=ccd(obj1, vel1, obj2, vel2, triLinearCCD);
	const auto endTime = steady_clock::now();
	std::cout<<"colTime: "<<t<<"\n";
	std::cout << "used seconds: " <<
		duration(endTime - initialTime).count()
		<< std::endl;

	obj1.moveObj(displace*t);
	obj1.writeObj("bunnyCol1.obj");
	obj2.writeObj("bunnyCol2.obj");
}

void sortupBunny(){
	std::ifstream in("low-bunny.obj");
	int num;
	std::string tit; 
	std::vector<Vector3d> verts;
	std::vector<Vector3i> faces;
	in>>tit;
	while(tit=="v"){
		Vector3d vert;
		in >> vert[0] >> vert[1] >> vert[2] >> tit;
		verts.push_back(vert/20.);
	}
	while(tit=="vn"){
		Vector3d vert;
		in >> vert[0] >> vert[1] >> vert[2] >> tit;
	}
	while(tit=="f"){
		Vector3i face;
		int placeholder;
		in >> face[0]>>tit>>placeholder;
		in >> face[1]>>tit>>placeholder;
		in >> face[2]>>tit>>placeholder;
		faces.push_back(face);
		in>>tit;
	}
	in.close();
	std::ofstream out("bunny292.obj");
	out<<verts.size()<<"  "<<faces.size()<<"\n";
	for(auto&vert:verts)
		out<<"v "<<vert[0]<<" "<<vert[1]<<" "<<vert[2]<<"\n";
	for(auto&face:faces)
		out<<"f "<<face[0]<<" "<<face[1]<<" "<<face[2]<<"\n";
	out.close();
}

void parabolaBunnyTorus(){
	// read in torus, bunny, both zero weights
	TriLinearMesh bunnyPos("bunny292.obj"), bunnyVel(bunnyPos.cntPatches);
	RatParamMesh<RecQuadRatBezier> torusPos, torusVel;
	generateTorusComponent(torusPos, torusVel);
	// torusPos.writeObj("bunny-torus/0.obj");

	using steady_clock = std::chrono::steady_clock;
	using duration = std::chrono::duration<double>;
	std::srand(0);
	Vector3d displace = Vector3d::Unit(0);//Random();
	Vector3d swirlAxisBunny = Vector3d::Random().normalized();
	Vector3d swirlAxisTorus = Vector3d::Random().normalized();
	displace[1] = 0;
	displace.normalize();
	displace *= 8;
	Vector3d initVel = displace*2;
	double swirlSpeed = PI/2.;
	bunnyPos.setOrigin(displace);
	torusPos.setOrigin(-displace);
	torusPos.rotateObj(-PI/2.,Vector3d::Unit(0),torusPos.getOrigin());
	bunnyPos.writeObj("bunny-torus/0-bunny.obj");
	torusPos.writeObj("bunny-torus/0-torus.obj");

	constexpr double totalTime = 1., deltaT = 0.02;
	constexpr int totalFrame = static_cast<int>(totalTime/deltaT);
	bool hasCol[totalFrame];
	double timeCost[totalFrame];
	double firstCol = totalTime;
	
	TriLinearMesh newBunnyPos = bunnyPos;
	RatParamMesh<RecQuadRatBezier> newTorusPos = torusPos;
	std::unique_ptr<SolverType> solver = std::make_unique<SolverTD>();
	for(int fr = 0; fr < totalFrame; fr++){
		Vector3d accel = (fr+1)*deltaT*Vector3d(0,-9.8,0);

		newBunnyPos.moveObj(deltaT*(-initVel+accel));
		newBunnyPos.rotateObj(swirlSpeed*deltaT,swirlAxisBunny,newBunnyPos.getOrigin());
		bunnyVel.setVel(bunnyPos, newBunnyPos, deltaT);

		newTorusPos.moveObj(deltaT*(initVel+accel));
		newTorusPos.rotateObj(-swirlSpeed*deltaT,swirlAxisTorus,newTorusPos.getOrigin());
		torusVel.setVel(torusPos, newTorusPos, deltaT);

		const auto initialTime = steady_clock::now();
		double t = ccd(bunnyPos, bunnyVel, torusPos, torusVel, 
				solver<TriLinearBezier,RecQuadRatBezier,TriParamBound,RecParamBound>->solveCCD, deltaT);
		const auto endTime = steady_clock::now();
		if(t < deltaT){
			hasCol[fr] = true;
			if(firstCol >= totalTime){
				firstCol = t + fr * deltaT;
				bunnyPos.moveObj(bunnyVel, t);
				torusPos.moveObj(torusVel, t);
				bunnyPos.writeObj("bunny-torus/col-bunny.obj");
				torusPos.writeObj("bunny-torus/col-torus.obj");
			}
		}
		else hasCol[fr] = false;
		timeCost[fr] = duration(endTime - initialTime).count();
		std::cout<<"frame "<<fr<<" costs "<<timeCost[fr]<<"s.\n";
		
		bunnyPos = newBunnyPos;
		torusPos = newTorusPos;

		// bunnyPos.writeObj("./bunny-torus/"+std::to_string(fr)+"-bunny.obj");
		// torusPos.writeObj("./bunny-torus/"+std::to_string(fr)+"-torus.obj");
	}
	std::cout<<"First collision time: "<< firstCol<<"s.\n";
	std::ofstream f("./bunny-torus/time_cost.txt");
	for(int i=0;i<totalFrame;i++)
		f<<hasCol[i]<<" "<<timeCost[i]<<"\n";
	f.close();
}

void parabolaPotCup(){
	// read in torus, bunny, both zero weights
	RecBezierMesh teapotPos("teapot32.txt"), teacupPos("teacup26.txt");
	RecBezierMesh teapotVel(teapotPos.cntPatches), teacupVel(teacupPos.cntPatches);
	// torusPos.writeObj("pot-cup/0.obj");

	using steady_clock = std::chrono::steady_clock;
	using duration = std::chrono::duration<double>;
	std::srand(0);
	Vector3d displace = Vector3d::Unit(0);//Random();
	Vector3d swirlAxisPot = Vector3d::Random().normalized();
	Vector3d swirlAxisCup = Vector3d::Random().normalized();
	displace[1] = 0;
	displace.normalize();
	displace *= 8;
	Vector3d initVel = displace*2;
	double swirlSpeed = PI/2.;
	teapotPos.setOrigin(-displace);
	teacupPos.setOrigin(displace);
	// teapotPos.rotateObj(-PI/2,(Vector3d::Random()).normalized(),teapotPos.getOrigin());
	// teacupPos.rotateObj(PI/2.,Vector3d::Unit(0),teacupPos.getOrigin());
	teapotPos.rotateObj(-PI/2,(Vector3d::Random()).normalized(),teapotPos.getOrigin());
	teacupPos.rotateObj(PI/2.,(Vector3d::Random()).normalized(),teacupPos.getOrigin());
	teapotPos.writeObj("pot-cup/0-pot.obj");
	teacupPos.writeObj("pot-cup/0-cup.obj");

	constexpr double totalTime = 1., deltaT = 0.02;
	constexpr int totalFrame = static_cast<int>(totalTime/deltaT);
	bool hasCol[totalFrame];
	double timeCost[totalFrame];
	double firstCol = totalTime;
	
	RecBezierMesh newTeapotPos = teapotPos;
	RecBezierMesh newTeacupPos = teacupPos;
	for(int fr = 0; fr < totalFrame; fr++){
		Vector3d accel = (fr+1)*deltaT*Vector3d(0,-9.8,0);

		newTeapotPos.moveObj(deltaT*(initVel+accel));
		newTeapotPos.rotateObj(-swirlSpeed*deltaT,swirlAxisPot,newTeapotPos.getOrigin());
		teapotVel.setVel(teapotPos, newTeapotPos, deltaT);

		newTeacupPos.moveObj(deltaT*(-initVel+accel));
		newTeacupPos.rotateObj(swirlSpeed*deltaT,swirlAxisCup,newTeacupPos.getOrigin());
		teacupVel.setVel(teacupPos, newTeacupPos, deltaT);

		const auto initialTime = steady_clock::now();
		double t = ccd(teapotPos, teapotVel, teacupPos, teacupVel, 
					solveCCD<RecCubicBezier,RecCubicBezier,RecParamBound,RecParamBound>, deltaT);
		std::cout<<t<<"\n";
		const auto endTime = steady_clock::now();
		if(t < deltaT){
			hasCol[fr] = true;
			if(firstCol >= totalTime){
				firstCol = t + fr * deltaT;
				teapotPos.moveObj(teapotVel, t);
				teacupPos.moveObj(teacupVel, t);
				// teapotPos.writeObj("pot-cup/col-pot.obj");
				// teacupPos.writeObj("pot-cup/col-cup.obj");
			}
		}
		else hasCol[fr] = false;
		timeCost[fr] = duration(endTime - initialTime).count();
		std::cout<<"frame "<<fr<<" costs "<<timeCost[fr]<<"s.\n";
		
		teapotPos = newTeapotPos;
		teacupPos = newTeacupPos;

		// teapotPos.writeObj("./pot-cup/"+std::to_string(fr)+"-pot.obj");
		// teacupPos.writeObj("./pot-cup/"+std::to_string(fr)+"-cup.obj");
	}
	std::cout<<"First collision time: "<< firstCol<<"s.\n";
	std::ofstream f("./pot-cup/time_cost.txt");
	for(int i=0;i<totalFrame;i++)
		f<<hasCol[i]<<" "<<timeCost[i]<<"\n";
	f.close();
}


template<typename ObjType, typename ParamType>
static void singleTest(){
	ObjType obj1, obj2, vel1, vel2;
	using steady_clock = std::chrono::steady_clock;
	using duration = std::chrono::duration<double>;

	readinDoFs<ObjType>(obj1.ctrlp, vel1.ctrlp, obj2.ctrlp, vel2.ctrlp);

	Array2d uv1,uv2;
	const auto initialTime = steady_clock::now();
	
	std::unique_ptr<SolverType> solver;
	switch (solverType){
		case SolverType::TDIntv:
			sovler = std::make_unique<SolverTD>();
			break;
		case SolverType::BaseIntv:
			sovler = std::make_unique<SolverBase>();
			break;
		default:
			std::cerr<<"";
			exit(-1);
	}

	double t = solver->solveCCD<ObjType, ObjType, ParamType, ParamType>(obj1,vel1,obj2,vel2,uv1,uv2, DeltaT);

	const auto endTime = steady_clock::now();
	std::cout<<cnt<<"\n";
	std::cout << "used seconds: " <<
		duration(endTime - initialTime).count()
		<< std::endl;
}

void boundaryTest(){
	enum class Kase {FF, EF, VF, EE, EV, VV};
	constexpr Kase kase = Kase::EE;
	RecBezierMesh obj(2);
	switch (kase){
		case Kase::EE:
			obj.patches[0].ctrlp = {
				Vector3d(0, 0, 0), Vector3d(1, 0, 0), Vector3d(2, 0, 0), Vector3d(3, 0, 0),
				Vector3d(0, 1, 0), Vector3d(1, 1, 0), Vector3d(2, 1, 0), Vector3d(3, 1, 0),
				Vector3d(0, 2, 0), Vector3d(1, 2, 0), Vector3d(2, 2, 0), Vector3d(3, 2, 0),
				Vector3d(0, 3, 0), Vector3d(1, 3, 0), Vector3d(2, 3, 0), Vector3d(3, 3, 0)
			}, 
			obj.patches[1].ctrlp = {
				Vector3d(0, 0, 0), Vector3d(1, 0, 0), Vector3d(2, 0, 0), Vector3d(3, 0, 0),
				Vector3d(0, 0, 1), Vector3d(1, 0, 1), Vector3d(2, 0, 1), Vector3d(3, 0, 1),
				Vector3d(0, 0, 2), Vector3d(1, 0, 2), Vector3d(2, 0, 2), Vector3d(3, 0, 2),
				Vector3d(0, 0, 3), Vector3d(1, 0, 3), Vector3d(2, 0, 3), Vector3d(3, 0, 3)
			};
			for(auto &p:obj.patches[0].ctrlp)p[0]-=2;
			for(auto &p:obj.patches[1].ctrlp)p[0]+=2;
			obj.writeObj("boundary/EE.obj");
			break;
		case Kase::EV:
			obj.patches[0].ctrlp = {
				Vector3d(0, 0, 0), Vector3d(1, 0, 0), Vector3d(2, 0, 0), Vector3d(3, 0, 0),
				Vector3d(0, 1, 0), Vector3d(1, 1, 0), Vector3d(2, 1, 0), Vector3d(3, 1, 0),
				Vector3d(0, 2, 0), Vector3d(1, 2, 0), Vector3d(2, 2, 0), Vector3d(3, 2, 0),
				Vector3d(0, 3, 0), Vector3d(1, 3, 0), Vector3d(2, 3, 0), Vector3d(3, 3, 0)
			}, 
			obj.patches[1].ctrlp = {
				Vector3d(0, 0, 0), Vector3d(1, 0, 0), Vector3d(2, 0, 0), Vector3d(3, 0, 0),
				Vector3d(0, 0, 1), Vector3d(1, 0, 1), Vector3d(2, 0, 1), Vector3d(3, 0, 1),
				Vector3d(0, 0, 2), Vector3d(1, 0, 2), Vector3d(2, 0, 2), Vector3d(3, 0, 2),
				Vector3d(0, 0, 3), Vector3d(1, 0, 3), Vector3d(2, 0, 3), Vector3d(3, 0, 3)
			};
			for(auto &p:obj.patches[0].ctrlp)p[0]-=2;
			for(auto &p:obj.patches[1].ctrlp)p[0]+=2;
			obj.writeObj("boundary/EE.txt");
			break;
		default:
			std::cerr<<"no such kase.\n";
	}
}
*/