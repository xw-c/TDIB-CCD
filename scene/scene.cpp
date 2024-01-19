#include"triMesh.h"
#include"recBezierMesh.h"
#include"recRatBezierMesh.h"
#include"paramCCD.h"

template<typename Mesh1, typename Mesh2, typename Func>
double ccd(const Mesh1& mesh1, const Mesh1& vel1,
			const Mesh2& mesh2, const Mesh2& vel2,
			Func&& paramCCD, 
			const double upperTime = DeltaT){
	double minTime = upperTime;
	Array2d uv1, uv2;
	if(mesh1.cntPatches!=vel1.cntPatches||mesh2.cntPatches!=vel2.cntPatches){
		std::cerr<<"dofs donot match!\n";
		exit(-1);
	}
	//问一下gpt：1、两个for auto patch同时遍历，2、为什么Func后面有两个&，3、为什么传参不能不写DeltaT
	for(int i = 0; i < mesh1.cntPatches; i++)
		for(int j = 0; j < mesh2.cntPatches; j++){
			double t = paramCCD(mesh1.patches[i], vel1.patches[i], mesh2.patches[j], vel2.patches[j], 
								uv1, uv2, minTime);
			// std::cout<<i<<" "<<j<<": "<<t<<"\n";
			if(t>=0){
				if(t==0)return 0;
				minTime = std::min(minTime, t);
			}
		}
	return minTime;
}

// take down the case settings
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
	RecRatBezierMesh obj1, vel1;
	generateTorusComponent(obj1, vel1);
	obj1.writeObj("torus.obj");

	RecRatBezierMesh obj2 = obj1, vel2 = vel1;
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

// void sortupBunny(){
// 	std::ifstream in("low-bunny.obj");
// 	int num;
// 	std::string tit; 
// 	std::vector<Vector3d> verts;
// 	std::vector<Vector3i> faces;
// 	in>>tit;
// 	while(tit=="v"){
// 		Vector3d vert;
// 		in >> vert[0] >> vert[1] >> vert[2] >> tit;
// 		verts.push_back(vert/20.);
// 	}
// 	while(tit=="vn"){
// 		Vector3d vert;
// 		in >> vert[0] >> vert[1] >> vert[2] >> tit;
// 	}
// 	while(tit=="f"){
// 		Vector3i face;
// 		int placeholder;
// 		in >> face[0]>>tit>>placeholder;
// 		in >> face[1]>>tit>>placeholder;
// 		in >> face[2]>>tit>>placeholder;
// 		faces.push_back(face);
// 		in>>tit;
// 	}
// 	in.close();
// 	std::ofstream out("bunny292.obj");
// 	out<<verts.size()<<"  "<<faces.size()<<"\n";
// 	for(auto&vert:verts)
// 		out<<"v "<<vert[0]<<" "<<vert[1]<<" "<<vert[2]<<"\n";
// 	for(auto&face:faces)
// 		out<<"f "<<face[0]<<" "<<face[1]<<" "<<face[2]<<"\n";
// 	out.close();
// }

void parabolaBunnyTorus(){
	// read in torus, bunny, both zero weights
	TriLinearMesh bunnyPos("bunny292.obj"), bunnyVel(bunnyPos.cntPatches);
	RecRatBezierMesh torusPos, torusVel;
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
	RecRatBezierMesh newTorusPos = torusPos;
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
				solveCCD<TriLinearBezier,RecQuadRatBezier,TriParamBound,RecParamBound>, deltaT);
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
		double t = ccd(teapotPos, teapotVel, teacupPos, teacupVel, recBezierCCD, deltaT);
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

int main(){
	parabolaPotCup();
	// testBunny();
	// sortupBunny();
}