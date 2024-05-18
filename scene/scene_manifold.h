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

template<typename ObjType, typename ParamType>
void planeTest(const std::string& solverType, const BoundingBoxType & bb, const double& deltaDist, const std::string& outputFile){
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
		t = solver.solveCCD(obj1,vel1,obj2,vel2,solutSet,bb,DeltaT,deltaDist);
	}else if(solverType=="base"){
		SolverBaseManifold<ObjType,ObjType,ParamType,ParamType> solver;
		t = solver.solveCCD(obj1,vel1,obj2,vel2,solutSet,bb,DeltaT,deltaDist);}
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
		file<<r.tIntv[0]<<" "<<v[0]<<" "<<v[1]<<"\n";
	}
	file.close();
	// fp.close();
	// ft.close();
}

template<typename ParamObj1, typename ParamObj2, typename ParamBound1, typename ParamBound2>
double ccd_ratBezier_manifold(const std::string& solverType, const std::string& outputFile,
			const RatParamMesh<ParamObj1>& mesh1, const RatParamMesh<ParamObj1>& vel1,
			const RatParamMesh<ParamObj2>& mesh2, const RatParamMesh<ParamObj2>& vel2,
			const BoundingBoxType & bb, 
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
								solutSet,bb, minTimeUB, deltaDist);
			}
			else if(solverType=="base"){
				SolverBaseManifold<ParamObj1,ParamObj2,ParamBound1,ParamBound2> solver;
				t = solver.solveCCD(mesh1.patches[i], vel1.patches[i], mesh2.patches[j], vel2.patches[j], 
								solutSet,bb, minTimeUB, deltaDist);
			}
			else{
				std::cerr<<"solver not implemented!\n";
				exit(-1);
			}
			std::cout<<i<<" "<<j<<": "<<t<<" "<<solutSet.size()<<"\n";
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

void torusTest(const std::string& solverType, const BoundingBoxType & bb, const double& deltaDist, const std::string& outputFile){
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
	double t = ccd_ratBezier_manifold<RecQuadRatBezier,RecQuadRatBezier,RecParamBound, RecParamBound>(solverType, outputFile, torusPos, torusVel, torusPos2, torusVel2,bb, DeltaT, deltaDist);
	const auto endTime = steady_clock::now();
	std::cout<<"colTime: "<<t<<"\n";
	std::cout << "used seconds: " <<
		duration(endTime - initialTime).count()
		<< std::endl;
}