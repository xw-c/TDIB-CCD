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
/*
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
*/
class TriLinearMesh{
public:
	int cntPatches;
	// int cntVerts, cntFaces;
	std::vector<TriLinearBezier> patches;
	// std::vector<Vector3d> verts;
	// std::vector<Vector3i> faces;

	// 用于初始化速度
	TriLinearMesh(const int _cntFaces): cntPatches(_cntFaces)
	// , cntVerts(_cntVerts), cntFaces(_cntFaces)
	{
		patches.resize(cntPatches);
		for(auto& patch: patches){
			for(auto& pt:patch.ctrlp)
				pt.setZero();
		}
	};
	// 读入模型
	TriLinearMesh(const std::string& filename){
		std::ifstream in(filename);
		int cntV;
		in >> cntV >> cntPatches;

		std::vector<Vector3d> verts;
		std::vector<Vector3i> faces;
		verts.resize(cntV);
		faces.resize(cntPatches);
		char c;
		for(auto& vert: verts)
			in >> c >> vert[0] >> vert[1] >> vert[2];
		for(auto& face: faces){
			in >> c >> face[0] >> face[1] >> face[2];
			std::array<Vector3d, 3> pts = {verts[face[0]-1], verts[face[1]-1], verts[face[2]-1]};
			patches.emplace_back(pts);
		}
		in.close();
	};
	void convert2Mesh(std::vector<Vector3d>&verts, std::vector<Vector3i>&faces) const {
		int cntVerts=0;
		for(auto& patch: patches){
			for(auto& pt: patch.ctrlp)
				verts.push_back(pt);
			faces.emplace_back(cntVerts+1, cntVerts+2, cntVerts+3);
			cntVerts+=3;
		}
	}
	void writeObj(const std::string& filename) const {
		std::vector<Vector3d> verts;
		std::vector<Vector3i> faces;
		convert2Mesh(verts, faces);
		std::ofstream out(filename);
		for(auto&vert:verts)
			out<<"v "<<vert[0]<<" "<<vert[1]<<" "<<vert[2]<<"\n";
		for(auto&face:faces)
			out<<"f "<<face[0]<<" "<<face[1]<<" "<<face[2]<<"\n";
		out.close();
	}
	void moveObj(const Vector3d& dis){
		for(auto& patch: patches)
			for(auto& pt:patch.ctrlp)
				pt+=dis;
	}
	void moveObj(const TriLinearMesh& vel, const double& t){
		for(int i=0;i<cntPatches;i++)
			for(int j=0;j<TriLinearBezier::cntCp;j++)
				patches[i].ctrlp[j]+=vel.patches[i].ctrlp[j]*t;
	}
	Vector3d getOrigin(){
		Vector3d center = Vector3d::Zero();
		for(const auto& patch:patches){
			for(const auto& p:patch.ctrlp){
				center+=p;
			}
		}
		center /= cntPatches * TriLinearBezier::cntCp;
		return center;
	}
	void setOrigin(const Vector3d& dis){
		Vector3d diff = getOrigin() - dis;
		for(auto& patch:patches){
			for(auto& p:patch.ctrlp){
				p -= diff;
			}
		}
	}
	void rotateObj(const double& angle, const Vector3d& axis, const Vector3d& origin = Vector3d::Zero()){
		//似乎带仿射
		Eigen::AngleAxisd rot(angle, axis);
		for(auto& patch: patches)
			for(auto& pt:patch.ctrlp)
				pt = (rot.matrix()*(pt-origin)+origin).eval();
	}
	void setVel(const TriLinearMesh& p1, const TriLinearMesh& p2, const double& dt){
		for(int i=0;i<cntPatches;i++)
			for(int j=0;j<TriLinearBezier::cntCp;j++)
				patches[i].ctrlp[j] = (p2.patches[i].ctrlp[j] - p1.patches[i].ctrlp[j]) / dt;
	}
};

void generateTorusComponent(RatParamMesh<RecQuadRatBezier>& torus, RatParamMesh<RecQuadRatBezier>& vel){
	double r=2, a=1;
	//顶点位置没有乘权重
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
		std::cerr<<"dofs do not match!\n";
		exit(-1);
	}
	//问一下gpt：1、两个for auto patch同时遍历，2、为什么Func后面有两个&，3、为什么传参不能不写DeltaT
	for(int i = 0; i < mesh1.cntPatches; i++)
		for(int j = 0; j < mesh2.cntPatches; j++){
			double t = paramCCD(mesh1.patches[i], vel1.patches[i], mesh2.patches[j], vel2.patches[j], 
								uv1, uv2, bb, minTime, deltaDist);
			std::cout<<i<<" "<<j<<": "<<t<<"\n";
			if(t>=0){
				if(t==0)return 0;
				minTime = std::min(minTime, t);
			}
		}
	return minTime;
}

void parabolaBunnyTorus(){
	// read in torus, bunny, both zero weights
	TriLinearMesh bunnyPos("bunny292.obj"), bunnyVel(bunnyPos.cntPatches);
	RatParamMesh<RecQuadRatBezier> torusPos, torusVel;
	generateTorusComponent(torusPos, torusVel);
	torusPos.writeObj("bunny-torus/0.obj");

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
	// std::unique_ptr<SolverType> solver = std::make_unique<SolverTD<TriLinearBezier,RecQuadRatBezier,TriParamBound,RecParamBound> >();
	for(int fr = 0; fr < totalFrame; fr++){
		Vector3d accel = (fr+1)*deltaT*Vector3d(0,-9.8,0);

		newBunnyPos.moveObj(deltaT*(-initVel+accel));
		newBunnyPos.rotateObj(swirlSpeed*deltaT,swirlAxisBunny,newBunnyPos.getOrigin());
		bunnyVel.setVel(bunnyPos, newBunnyPos, deltaT);

		newTorusPos.moveObj(deltaT*(initVel+accel));
		newTorusPos.rotateObj(-swirlSpeed*deltaT,swirlAxisTorus,newTorusPos.getOrigin());
		torusVel.setVel(torusPos, newTorusPos, deltaT);

		const auto initialTime = steady_clock::now();
		double t = ccd(bunnyPos, bunnyVel, torusPos, torusVel, SolverTD<TriLinearBezier,RecQuadRatBezier,TriParamBound,RecParamBound>::solveCCD);
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

void compareInLocking(const std::string& solverType, const BoundingBoxType & bb, 
					const double& deltaDist, const std::string& outputFile){
	ParamMesh<RecCubicBezier> cloth(100), clothVel(100);
	Array2d uv1, uv2;
	using steady_clock = std::chrono::steady_clock;
	using duration = std::chrono::duration<double>;
	std::ifstream readin;
	// 不同时间步长、不同solver、不同文件的cloth
	const int fileTotal = 1, stepTotal=1;
	double costs[fileTotal][stepTotal];
	int frames[fileTotal]={28};//,120,140,160,180,207,213
	double timeSteps[stepTotal]={0.002};//0.001,0.002,0.005,0.01,0.02
	for(int frcnt=0;frcnt<fileTotal;frcnt++){
		int fr = frames[frcnt];
		std::string filename = "locking-"+std::to_string(fr)+".dat";
		readin.open(filename, std::ios::binary);
		for(int id=0;id<100;id++){
			for (int i = 0; i < 16; i++)
				for(int k=0;k<3;k++)
					readin.read(reinterpret_cast<char *>(&cloth.patches[id].ctrlp[i][k]), sizeof(double));
			for (int i = 0; i < 16; i++)
				for(int k=0;k<3;k++)
					readin.read(reinterpret_cast<char *>(&clothVel.patches[id].ctrlp[i][k]), sizeof(double));
		}
		readin.close();
		for(int stcnt=0;stcnt<stepTotal;stcnt++){
			auto h=timeSteps[stcnt];
			std::cout<<"frame"<< fr <<"time step "<<h<<": ";
			double t, minT = 0.02;
			const auto initialTime = steady_clock::now();
			for(int id1=0;id1<100;id1++){
				for(int id2=0;id2<100;id2++){
					int i=id1/10,j=id1%10;
					int ii=id2/10,jj=id2%10;
					if (ii == i || ii == i - 1 || ii == i + 1)
						if(jj == j || jj == j - 1 || jj == j + 1)
							continue;
					if (ii >= i && jj >= j) 
						continue;
					if(solverType=="td")
						t = SolverTD<RecCubicBezier,RecCubicBezier,RecParamBound,RecParamBound>::solveCCD(cloth.patches[id1],clothVel.patches[id1],cloth.patches[id2],clothVel.patches[id2],uv1,uv2,bb,minT,deltaDist);
					else if(solverType=="base")
						t = SolverBase<RecCubicBezier,RecCubicBezier,RecParamBound,RecParamBound>::solveCCD(cloth.patches[id1],clothVel.patches[id1],cloth.patches[id2],clothVel.patches[id2],uv1,uv2,bb,minT,deltaDist);
					else{
						std::cerr<<"solver not implemented!\n";
						exit(-1);
					}
					if(t!=-1)minT=std::min(minT,t);
				}
			}
			const auto endTime = steady_clock::now();
			costs[frcnt][stcnt]=duration(endTime - initialTime).count();
			std::cout<<"toi " << minT <<", costs "<<duration(endTime - initialTime).count()<<"s.\n";
		}
	}
	std::cout<<std::fixed<<std::setprecision(2);
	for(int i=0;i<fileTotal;i++){
		std::cout<<"\\hline\nFr. "<<frames[i];
		for(int j=0;j<stepTotal;j++)
			std::cout<<"&"<<costs[i][j];
		std::cout<<"\\\\\n";
	}
}


void compareInTeapot(const std::string& solverType, const BoundingBoxType & bb, 
					const double& deltaDist, const std::string& outputFile){
	ParamMesh<RecCubicBezier> cloth(400), clothVel(400);
	ParamMesh<RecCubicBezier> teapot(32), teapotVel(32);

	Array2d uv1, uv2;
	using steady_clock = std::chrono::steady_clock;
	using duration = std::chrono::duration<double>;
	std::ifstream readin;
	readin.open("teapot-collider.dat", std::ios::binary);
	for(int id=0;id<32;id++)
		for (int i = 0; i < 16; i++)
			for(int k=0;k<3;k++)
				readin.read(reinterpret_cast<char *>(&teapot.patches[id].ctrlp[i][k]), sizeof(double));
	for(int id=0;id<32;id++)
		for (int i = 0; i < 16; i++)
			for(int k=0;k<3;k++)
				readin.read(reinterpret_cast<char *>(&teapotVel.patches[id].ctrlp[i][k]), sizeof(double));
	readin.close();
	// 不同时间步长、不同solver、不同文件的cloth
	const int fileTotal = 6, stepTotal=5;
	double costs[fileTotal][stepTotal];
	int frames[fileTotal]={94,113,126,278,356,389};
	double timeSteps[stepTotal]={0.001,0.002,0.005,0.01,0.02};
	for(int frcnt=0;frcnt<fileTotal;frcnt++){
		int fr = frames[frcnt];
		std::string filename = "teapot-cloth-"+std::to_string(fr)+".dat";
		readin.open(filename, std::ios::binary);
		for(int id=0;id<400;id++){
			for (int i = 0; i < 16; i++)
				for(int k=0;k<3;k++)
					readin.read(reinterpret_cast<char *>(&cloth.patches[id].ctrlp[i][k]), sizeof(double));
			for (int i = 0; i < 16; i++)
				for(int k=0;k<3;k++)
					readin.read(reinterpret_cast<char *>(&clothVel.patches[id].ctrlp[i][k]), sizeof(double));
		}
		readin.close();
		// cloth.writeObj("teapot-cloth.obj",0.05,0.05);
		// teapot.writeObj("teapot-collider.obj",0.05,0.05);
		for(int stcnt=0;stcnt<stepTotal;stcnt++){
			auto h=timeSteps[stcnt];
			std::cout<<"frame"<< fr <<"time step "<<h<<": ";
			double t, minT=h;
			const auto initialTime = steady_clock::now();
			for(int id1=0;id1<400;id1++){
				for(int id2=0;id2<32;id2++){
					// if(SHOWANS)std::cout<<i<<", "<<j<<": "<<ii<<", "<<jj<<":\n";
					if(solverType=="td")
						t = SolverTD<RecCubicBezier,RecCubicBezier,RecParamBound,RecParamBound>::solveCCD(cloth.patches[id1],clothVel.patches[id1],teapot.patches[id2],teapotVel.patches[id2],uv1,uv2,bb,minT,deltaDist);
					else if(solverType=="base")
						t = SolverBase<RecCubicBezier,RecCubicBezier,RecParamBound,RecParamBound>::solveCCD(cloth.patches[id1],clothVel.patches[id1],teapot.patches[id2],teapotVel.patches[id2],uv1,uv2,bb,minT,deltaDist);
					else{
						std::cerr<<"solver not implemented!\n";
						exit(-1);
					}
					if(t!=-1)minT=std::min(minT,t);
				}
			}
			const auto endTime = steady_clock::now();
			costs[frcnt][stcnt]=duration(endTime - initialTime).count();
			std::cout<<"toi " << minT <<", costs "<<duration(endTime - initialTime).count()<<"s.\n";
		}
	}
	std::cout<<std::fixed<<std::setprecision(2);
	for(int i=0;i<fileTotal;i++){
		std::cout<<"\\hline\nFr. "<<frames[i];
		for(int j=0;j<stepTotal;j++)
			std::cout<<"&"<<costs[i][j];
		std::cout<<"\\\\\n";
	}
}

/*
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

*/