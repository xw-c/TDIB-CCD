#pragma once
#include"paramMesh.h"
#include"solverBase.h"
#include"solverBaseManifold.h"
#include"solverTD.h"
#include"solverInitTD.h"
#include"solverRobustTD.h"
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
void randomTest(const std::string& solverType, const BoundingBoxType & bb,
				const double& deltaDist, const int& kase, const double& velMag, 
				const std::string& outputFile){
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
			t = SolverTD<ObjType,ObjType,ParamType,ParamType>::solveCCD(obj1,vel1,obj2,vel2,uv1,uv2,bb,DeltaT,deltaDist);
		else if(solverType=="base")
			t = SolverBase<ObjType,ObjType,ParamType,ParamType>::solveCCD(obj1,vel1,obj2,vel2,uv1,uv2,bb,DeltaT,deltaDist);
		else{
			std::cerr<<"solver not implemented!\n";
			exit(-1);
		}
		if(t>=0)hasCol++;
		// if(k==19){
		// 	saveDoFs<ObjType,ObjType>(obj1.ctrlp, vel1.ctrlp, obj2.ctrlp, vel2.ctrlp);
			// file<<t<<"\n";
		// }
		// ft<<t<<"\n";
		// if(t>0)file<<t<<"  "<<uv1[0]<<"  "<<uv1[1]<<"  "<<uv2[0]<<"  "<<uv2[1]<<"\n";
		std::cout<<"case "<<k<<" done "<<calcDist(obj1,vel1,obj2,vel2,uv1,uv2,t)<<"\n";
		
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

void FNCase(const std::string& solverType, const BoundingBoxType & bb,
				const double& deltaDist, const int& kase, const double& velMag, 
				const std::string& outputFile){
	RecLinearBezier obj1, obj2, vel1, vel2;
	double t;
	Array2d uv1,uv2;

	// 一个由CH计算导致的fn
	obj1.ctrlp = {
		Vector3d(-1e-20,0,0), Vector3d(-1e-20,1,0),
		Vector3d(-1,0,1), Vector3d(-1,1,1)
	}, 
	vel1.ctrlp = {
		Vector3d(0,0,0), Vector3d(0,0,0), 
		Vector3d(1,0,0), Vector3d(1,0,0)
	},
	obj2.ctrlp = {
		Vector3d(1,0,0), Vector3d(1,1,0),
		Vector3d(-1e-30,0,1), Vector3d(-1e-30,1,1)
	}, 
	vel2.ctrlp = {
		Vector3d(0,0,0), Vector3d(0,0,0), 
		Vector3d(0,0,0), Vector3d(0,0,0)
	};

	// readinDoFs<ObjType,ObjType>(obj1.ctrlp, vel1.ctrlp, obj2.ctrlp, vel2.ctrlp);
	using steady_clock = std::chrono::steady_clock;
	using duration = std::chrono::duration<double>;
	const auto initialTime = steady_clock::now();
	if(solverType=="td")
		t = SolverTD<RecLinearBezier,RecLinearBezier,RecParamBound,RecParamBound>::solveCCD(obj1,vel1,obj2,vel2,uv1,uv2,bb,DeltaT,deltaDist);
	else if(solverType=="robust")
			t = SolverRobustTD<RecLinearBezier,RecLinearBezier,RecParamBound,RecParamBound>::solveCCD(obj1,vel1,obj2,vel2,uv1,uv2,bb,DeltaT,deltaDist);
	else if(solverType=="base")
		t = SolverBase<RecLinearBezier,RecLinearBezier,RecParamBound,RecParamBound>::solveCCD(obj1,vel1,obj2,vel2,uv1,uv2,bb,DeltaT,deltaDist);
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


template<typename ObjType, typename ParamType>
void validate(const std::string& solverType, const BoundingBoxType & bb,
				const double& deltaDist, const int& kase, const double& velMag, 
				const std::string& outputFile){
	ObjType obj1, obj2, vel1, vel2;
	int hasCol = 0;
	double t;
	Array2d uv1,uv2;
	std::multiset<CCDIntv<ParamType, ParamType> > solutSet;
	solutSet.clear();
	
	// {obj1.ctrlp = {
	// 	Vector3d(0,0,0), Vector3d(0,1,1), Vector3d(0,2,2), Vector3d(0,3,3),
	// 	Vector3d(1,0,1), Vector3d(1,1,2), Vector3d(1,2,3), Vector3d(1,3,4),
	// 	Vector3d(2,0,2), Vector3d(2,1,3), Vector3d(2,2,4), Vector3d(2,3,5),
	// 	Vector3d(3,0,3), Vector3d(3,1,4), Vector3d(3,2,5), Vector3d(3,3,6)
	// }, 
	// vel1.ctrlp = {Vector3d(0,0,1), Vector3d(0,0,1), Vector3d(0,0,1), Vector3d(0,0,1), 
	// 	Vector3d(0,0,1), Vector3d(0,0,1), Vector3d(0,0,1), Vector3d(0,0,1),
	// 	Vector3d(0,0,1), Vector3d(0,0,1), Vector3d(0,0,1), Vector3d(0,0,1),
	// 	Vector3d(0,0,1), Vector3d(0,0,1), Vector3d(0,0,1), Vector3d(0,0,1)};
	// obj2.ctrlp = {
	// 	Vector3d(0,0,7.1), Vector3d(0,1,7.1), Vector3d(0,2,7.1), Vector3d(0,3,7.1),
	// 	Vector3d(1,0,7.1), Vector3d(1,1,7.1), Vector3d(1,2,7.1), Vector3d(1,3,7.1),
	// 	Vector3d(2,0,7.1), Vector3d(2,1,7.1), Vector3d(2,2,7.1), Vector3d(2,3,7.1),
	// 	Vector3d(3,0,7.1), Vector3d(3,1,7.1), Vector3d(3,2,7.1), Vector3d(3,3,7.1)
	// }, 
	// vel2.ctrlp = {Vector3d(0,0,-1), Vector3d(0,0,-1), Vector3d(0,0,-1), Vector3d(0,0,-1),
	// 	Vector3d(0,0,-1), Vector3d(0,0,-1), Vector3d(0,0,-1), Vector3d(0,0,-1),
	// 	Vector3d(0,0,-1), Vector3d(0,0,-1), Vector3d(0,0,-1), Vector3d(0,0,-1),
	// 	Vector3d(0,0,-1), Vector3d(0,0,-1), Vector3d(0,0,-1), Vector3d(0,0,-1)};}

	readinDoFs<ObjType,ObjType>(obj1.ctrlp, vel1.ctrlp, obj2.ctrlp, vel2.ctrlp);
	using steady_clock = std::chrono::steady_clock;
	using duration = std::chrono::duration<double>;
	const auto initialTime = steady_clock::now();
	if(solverType=="td")
		t = SolverTD<ObjType,ObjType,ParamType,ParamType>::solveCCD(obj1,vel1,obj2,vel2,uv1,uv2,bb,DeltaT,deltaDist);
	else if(solverType=="robust")
			t = SolverRobustTD<ObjType,ObjType,ParamType,ParamType>::solveCCD(obj1,vel1,obj2,vel2,uv1,uv2,bb,DeltaT,deltaDist);
	else if(solverType=="base")
		t = SolverBase<ObjType,ObjType,ParamType,ParamType>::solveCCD(obj1,vel1,obj2,vel2,uv1,uv2,bb,DeltaT,deltaDist);
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

// template<typename Mesh1, typename Mesh2, typename Func>
// double ccd(const Mesh1& mesh1, const Mesh1& vel1,
// 			const Mesh2& mesh2, const Mesh2& vel2,
// 			Func&& paramCCD, 
// 			const double upperTime = DeltaT){
// 	double minTime = upperTime;
// 	Array2d uv1, uv2;
// 	if(mesh1.cntPatches!=vel1.cntPatches||mesh2.cntPatches!=vel2.cntPatches){
// 		std::cerr<<"dofs do not match!\n";
// 		exit(-1);
// 	}
// 	//问一下gpt：1、两个for auto patch同时遍历，2、为什么Func后面有两个&，3、为什么传参不能不写DeltaT
// 	for(int i = 0; i < mesh1.cntPatches; i++)
// 		for(int j = 0; j < mesh2.cntPatches; j++){
// 			double t = paramCCD(mesh1.patches[i], vel1.patches[i], mesh2.patches[j], vel2.patches[j], 
// 								uv1, uv2, minTime);
// 			std::cout<<i<<" "<<j<<": "<<t<<"\n";
// 			if(t>=0){
// 				if(t==0)return 0;
// 				minTime = std::min(minTime, t);
// 			}
// 		}
// 	return minTime;
// }

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


// template<typename ObjType, typename ParamType>
// static void singleTest(){
// 	ObjType obj1, obj2, vel1, vel2;
// 	using steady_clock = std::chrono::steady_clock;
// 	using duration = std::chrono::duration<double>;

// 	readinDoFs<ObjType>(obj1.ctrlp, vel1.ctrlp, obj2.ctrlp, vel2.ctrlp);

// 	Array2d uv1,uv2;
// 	const auto initialTime = steady_clock::now();
	
// 	std::unique_ptr<SolverType> solver;
// 	switch (solverType){
// 		case SolverType::TDIntv:
// 			sovler = std::make_unique<SolverTD>();
// 			break;
// 		case SolverType::BaseIntv:
// 			sovler = std::make_unique<SolverBase>();
// 			break;
// 		default:
// 			std::cerr<<"";
// 			exit(-1);
// 	}

// 	double t = solver->solveCCD<ObjType, ObjType, ParamType, ParamType>(obj1,vel1,obj2,vel2,uv1,uv2, DeltaT);

// 	const auto endTime = steady_clock::now();
// 	std::cout<<cnt<<"\n";
// 	std::cout << "used seconds: " <<
// 		duration(endTime - initialTime).count()
// 		<< std::endl;
// }

// void boundaryTest(){
// 	enum class Kase {FF, EF, VF, EE, EV, VV};
// 	constexpr Kase kase = Kase::EE;
// 	RecBezierMesh obj(2);
// 	switch (kase){
// 		case Kase::EE:
// 			obj.patches[0].ctrlp = {
// 				Vector3d(0, 0, 0), Vector3d(1, 0, 0), Vector3d(2, 0, 0), Vector3d(3, 0, 0),
// 				Vector3d(0, 1, 0), Vector3d(1, 1, 0), Vector3d(2, 1, 0), Vector3d(3, 1, 0),
// 				Vector3d(0, 2, 0), Vector3d(1, 2, 0), Vector3d(2, 2, 0), Vector3d(3, 2, 0),
// 				Vector3d(0, 3, 0), Vector3d(1, 3, 0), Vector3d(2, 3, 0), Vector3d(3, 3, 0)
// 			}, 
// 			obj.patches[1].ctrlp = {
// 				Vector3d(0, 0, 0), Vector3d(1, 0, 0), Vector3d(2, 0, 0), Vector3d(3, 0, 0),
// 				Vector3d(0, 0, 1), Vector3d(1, 0, 1), Vector3d(2, 0, 1), Vector3d(3, 0, 1),
// 				Vector3d(0, 0, 2), Vector3d(1, 0, 2), Vector3d(2, 0, 2), Vector3d(3, 0, 2),
// 				Vector3d(0, 0, 3), Vector3d(1, 0, 3), Vector3d(2, 0, 3), Vector3d(3, 0, 3)
// 			};
// 			for(auto &p:obj.patches[0].ctrlp)p[0]-=2;
// 			for(auto &p:obj.patches[1].ctrlp)p[0]+=2;
// 			obj.writeObj("boundary/EE.obj");
// 			break;
// 		case Kase::EV:
// 			obj.patches[0].ctrlp = {
// 				Vector3d(0, 0, 0), Vector3d(1, 0, 0), Vector3d(2, 0, 0), Vector3d(3, 0, 0),
// 				Vector3d(0, 1, 0), Vector3d(1, 1, 0), Vector3d(2, 1, 0), Vector3d(3, 1, 0),
// 				Vector3d(0, 2, 0), Vector3d(1, 2, 0), Vector3d(2, 2, 0), Vector3d(3, 2, 0),
// 				Vector3d(0, 3, 0), Vector3d(1, 3, 0), Vector3d(2, 3, 0), Vector3d(3, 3, 0)
// 			}, 
// 			obj.patches[1].ctrlp = {
// 				Vector3d(0, 0, 0), Vector3d(1, 0, 0), Vector3d(2, 0, 0), Vector3d(3, 0, 0),
// 				Vector3d(0, 0, 1), Vector3d(1, 0, 1), Vector3d(2, 0, 1), Vector3d(3, 0, 1),
// 				Vector3d(0, 0, 2), Vector3d(1, 0, 2), Vector3d(2, 0, 2), Vector3d(3, 0, 2),
// 				Vector3d(0, 0, 3), Vector3d(1, 0, 3), Vector3d(2, 0, 3), Vector3d(3, 0, 3)
// 			};
// 			for(auto &p:obj.patches[0].ctrlp)p[0]-=2;
// 			for(auto &p:obj.patches[1].ctrlp)p[0]+=2;
// 			obj.writeObj("boundary/EE.txt");
// 			break;
// 		default:
// 			std::cerr<<"no such kase.\n";
// 	}
// }