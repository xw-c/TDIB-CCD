#pragma once
#include "paramMesh.h"
#include "solverTrad.h"
#include "solverTD.h"
#include "utils.h"
#include "triBezier.h"
#include "triRatBezier.h"
#include "recBezier.h"
#include "recRatBezier.h"


template<typename ObjType>
static void generatePatchPair(std::array<Vector3d, ObjType::cntCp> &CpPos1, std::array<Vector3d, ObjType::cntCp> &CpVel1,
									std::array<Vector3d, ObjType::cntCp> &CpPos2, std::array<Vector3d, ObjType::cntCp> &CpVel2, const double& velMagnitude = 1){
	Vector3d dir=Vector3d::Random().normalized();
	for (int i = 0; i < ObjType::cntCp; i++) {
		CpPos1[i] = Vector3d::Random() - dir;
		CpVel1[i] = Vector3d::Random() + dir * velMagnitude;
		CpPos2[i] = Vector3d::Random() + dir;
		CpVel2[i] = Vector3d::Random() - dir * velMagnitude;
	}
}

template<typename ObjType, typename ParamType>
void randomTest(const SolverType& solver, const BoundingBoxType & bb,
				const double& deltaDist, const int& kase){
	ObjType pos1, pos2, vel1, vel2;
	std::srand(0);
	int hasCol = 0;
	double t;
	Array2d uv1,uv2;

	using steady_clock = std::chrono::steady_clock;
	using duration = std::chrono::duration<double>;
	const auto initialTime = steady_clock::now();
	for(int k = 0; k < kase; k ++){
		generatePatchPair<ObjType>(pos1.ctrlp, vel1.ctrlp, pos2.ctrlp, vel2.ctrlp);
		if(solver==SolverType::TDIntv)
			t = SolverTD<ObjType,ObjType,ParamType,ParamType>::solveCCD(pos1,vel1,pos2,vel2,uv1,uv2,bb,deltaDist);
		else if(solver==SolverType::TradIntv)
			t = SolverTrad<ObjType,ObjType,ParamType,ParamType>::solveCCD(pos1,vel1,pos2,vel2,uv1,uv2,bb,deltaDist);
		else{
			std::cerr<<"solver not implemented!\n";
			exit(-1);
		}
		if(t>=0)hasCol++;
		if(!(k/100))std::cout<<"case "<<k<<" done.\n";
	}
	const auto endTime = steady_clock::now();
	std::cout << hasCol<<" pairs have collided.\n";
	std::cout << "average seconds: " <<
		duration(endTime - initialTime).count()/kase
		<< std::endl;
}

template<typename ObjType, typename ParamType>
void singleTest(const SolverType& solver, const BoundingBoxType & bb,
				const double& deltaDist){
	ObjType pos1, pos2, vel1, vel2;
	readinDoFs<ObjType,ObjType>(obj1.ctrlp, vel1.ctrlp, obj2.ctrlp, vel2.ctrlp);
	int hasCol = 0;
	double t;
	Array2d uv1,uv2;

	using steady_clock = std::chrono::steady_clock;
	using duration = std::chrono::duration<double>;
	const auto initialTime = steady_clock::now();
	if(solver==SolverType::TDIntv)
		t = SolverTD<ObjType,ObjType,ParamType,ParamType>::solveCCD(pos1,vel1,pos2,vel2,uv1,uv2,bb,deltaDist);
	else if(solver==SolverType::TradIntv)
		t = SolverTrad<ObjType,ObjType,ParamType,ParamType>::solveCCD(pos1,vel1,pos2,vel2,uv1,uv2,bb,deltaDist);
	else{
		std::cerr<<"solver not implemented!\n";
		exit(-1);
	}
	const auto endTime = steady_clock::now();
	std::cout << "used seconds: " <<
		duration(endTime - initialTime).count()
		<< std::endl;
	Vector3d const p1 = pos1.evaluatePatchPoint(uv1);
	Vector3d const v1 = vel1.evaluatePatchPoint(uv1);
	Vector3d const p2 = pos2.evaluatePatchPoint(uv2);
	Vector3d const v2 = vel2.evaluatePatchPoint(uv2);
	Vector3d const pt1=(v1*t+p1), pt2=(v2*t+p2);
	std::cout<<"distance residual: "<<(pt2-pt1).norm()<<"\n";
}
