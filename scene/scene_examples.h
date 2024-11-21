#pragma once
#include "paramMesh.h"
#include "solverTrad.h"
#include "solverTD.h"
#include "utils.h"
#include "triBezier.h"
#include "triRatBezier.h"
#include "recBezier.h"
#include "recRatBezier.h"

void generateTorusComponent(RatParamMesh<RecQuadRatBezier>& torus, RatParamMesh<RecQuadRatBezier>& vel){
	double r=2, a=1;
	std::array<double, 9> weight = {1,1,2,1,1,2,2,2,4};
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
void readInBunny(ParamMesh<TriLinearBezier>& bunnyPos, ParamMesh<TriLinearBezier>& bunnyVel){
	std::ifstream in("bunny292.obj");
	int cntVerts = 148, cntFaces = 292;
	bunnyVel = ParamMesh<TriLinearBezier>(292);
	bunnyPos.cntPatches = cntFaces;

	std::vector<Vector3d> verts;
	std::vector<Vector3i> faces;
	verts.resize(cntVerts);
	faces.resize(cntFaces);
	char c;
	for(auto& vert: verts)
		in >> c >> vert[0] >> vert[1] >> vert[2];
	for(auto& face: faces){
		in >> c >> face[0] >> face[1] >> face[2];
		std::array<Vector3d, 3> pts = {verts[face[0]-1], verts[face[1]-1], verts[face[2]-1]};
		bunnyPos.patches.emplace_back(pts);
	}
	in.close();
}
template<typename Mesh1, typename Mesh2, typename Func>
double ccd(const Mesh1& mesh1, const Mesh1& vel1,
			const Mesh2& mesh2, const Mesh2& vel2,
			Func&& paramCCD,  
			const BoundingBoxType & bb = BoundingBoxType::OBB,
			const double upperTime = DeltaT,
			const double& deltaDist = 1e-6){
	double minTime = upperTime;
	Array2d uv1, uv2;
	if(mesh1.cntPatches!=vel1.cntPatches||mesh2.cntPatches!=vel2.cntPatches){
		std::cerr<<"Error: numbers of dof for pos and vel do not match!\n";
		exit(-1);
	}
	for(int i = 0; i < mesh1.cntPatches; i++)
		for(int j = 0; j < mesh2.cntPatches; j++){
			double t = paramCCD(mesh1.patches[i], vel1.patches[i], mesh2.patches[j], vel2.patches[j], 
								uv1, uv2, bb, minTime, deltaDist);
			if(t>=0){
				if(t==0)return 0;
				minTime = std::min(minTime, t);
			}
		}
	return minTime;
}

void parabolaBunnyTorus(){
	ParamMesh<TriLinearBezier> bunnyPos, bunnyVel;
	RatParamMesh<RecQuadRatBezier> torusPos, torusVel;
	readInBunny(bunnyPos, bunnyVel);
	generateTorusComponent(torusPos, torusVel);

	using steady_clock = std::chrono::steady_clock;
	using duration = std::chrono::duration<double>;

	std::srand(0);
	Vector3d displace = Vector3d::Unit(0);
	Vector3d swirlAxisBunny = Vector3d::Random().normalized();
	Vector3d swirlAxisTorus = Vector3d::Random().normalized();
	displace *= 8;
	Vector3d initVel = displace*2;
	double swirlSpeed = PI/2.;
	bunnyPos.setOrigin(displace);
	torusPos.setOrigin(-displace);
	torusPos.rotateObj(-PI/2.,Vector3d::Unit(0),torusPos.getOrigin());
	bunnyPos.writeObj("start-bunny.obj");
	torusPos.writeObj("start-torus.obj");

	constexpr double totalTime = 1., deltaT = 0.02;
	constexpr int totalFrame = static_cast<int>(totalTime/deltaT);
	bool hasCol[totalFrame];
	double timeCost[totalFrame];
	double firstCol = totalTime;
	
	TriLinearMesh newBunnyPos = bunnyPos;
	RatParamMesh<RecQuadRatBezier> newTorusPos = torusPos;
	for(int fr = 0; fr < totalFrame; fr++){
		Vector3d accel = (fr+1)*deltaT*Vector3d(0,-9.8,0);

		newBunnyPos.moveObj(deltaT*(-initVel+accel));
		newBunnyPos.rotateObj(swirlSpeed*deltaT,swirlAxisBunny,newBunnyPos.getOrigin());
		bunnyVel.setVel(bunnyPos, newBunnyPos, deltaT);

		newTorusPos.moveObj(deltaT*(initVel+accel));
		newTorusPos.rotateObj(-swirlSpeed*deltaT,swirlAxisTorus,newTorusPos.getOrigin());
		torusVel.setVel(torusPos, newTorusPos, deltaT);
		const auto initialTime = steady_clock::now();
		double t = ccd(bunnyPos, bunnyVel, torusPos, torusVel, SolverTD<TriLinearBezier,RecQuadRatBezier,TriParamBound,RecParamBound>::solveCCD, BoundingBoxType::OBB,deltaT);
		const auto endTime = steady_clock::now();
		if(t < deltaT){
			hasCol[fr] = true;
			if(firstCol >= totalTime){
				firstCol = t + fr * deltaT;
				bunnyPos.moveObj(bunnyVel, t);
				torusPos.moveObj(torusVel, t);
				bunnyPos.writeObj("col-bunny.obj");
				torusPos.writeObj("col-torus.obj");
			}
		}
		else hasCol[fr] = false;
		timeCost[fr] = duration(endTime - initialTime).count();
		std::cout<<"frame "<<fr<<" costs "<<timeCost[fr]<<"s.\n";
		bunnyPos = newBunnyPos;
		torusPos = newTorusPos;
	}
	bunnyPos.writeObj("end-bunny.obj");
	torusPos.writeObj("end-torus.obj");
	std::cout<<"First collision time: "<< firstCol<<"s.\n";
	// std::ofstream f("./bunny-torus/time_cost.txt");
	// for(int i=0;i<totalFrame;i++)
	// 	f<<hasCol[i]<<" "<<timeCost[i]<<"\n";
	// f.close();
}
