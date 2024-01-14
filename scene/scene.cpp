#include"triMesh.h"
#include"writeBezierObj.h"
#include"writeNURBSObj.h"
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
			// if(i==j)continue;
			if(i%100==0&&j==0)std::cout<<i<<" "<<j<<" : "<<minTime<<" \n";
			double t = paramCCD(mesh1.patches[i], vel1.patches[i], mesh2.patches[j], vel2.patches[j], 
								uv1, uv2, BoundingBoxType::OBB, minTime);
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

int main(){
	testBunny();
	// sortupBunny();
}