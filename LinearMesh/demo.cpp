#include "../core/argsParser.h"
#include "linearGeom.h"
#include "linearSolver.h"
template<typename ObjType1, typename ObjType2>
static void saveDoFs(const std::array<Vector3d, ObjType1::cntCp>& CpPos1, 
					const std::array<Vector3d, ObjType1::cntCp>& CpVel1,
					const std::array<Vector3d, ObjType2::cntCp>& CpPos2, 
					const std::array<Vector3d, ObjType2::cntCp>& CpVel2,
					const std::string& filename="DoFs.txt"){
	std::ofstream f(filename);
	f<< std::fixed << std::setprecision(10);
	for(auto item : CpPos1)
		f<<item.transpose()<<"\n";
	for(auto item : CpVel1)
		f<<item.transpose()<<"\n";
	for(auto item : CpPos2)
		f<<item.transpose()<<"\n";
	for(auto item : CpVel2)
		f<<item.transpose()<<"\n";
	f.close();
}
static void generateEEPair(std::array<Vector3d, Edge::cntCp> &CpPos1, std::array<Vector3d, Edge::cntCp> &CpVel1,
						std::array<Vector3d, Edge::cntCp> &CpPos2, std::array<Vector3d, Edge::cntCp> &CpVel2){
	Vector3d dir=Vector3d::Random();
	for (int i = 0; i < Edge::cntCp; i++) {
		CpPos1[i] = Vector3d::Random() - dir;
		CpVel1[i] = Vector3d::Random() + dir;
		CpPos2[i] = Vector3d::Random() + dir;
		CpVel2[i] = Vector3d::Random() - dir;
	}
}
static void generateVFPair(Vector3d &CpPos1, Vector3d &CpVel1,
						std::array<Vector3d, Face::cntCp> &CpPos2, std::array<Vector3d, Face::cntCp> &CpVel2){
	Vector3d dir=Vector3d::Random();
	CpPos1 = Vector3d::Random() - dir;
	CpVel1 = Vector3d::Random() + dir;
	for (int i = 0; i < Face::cntCp; i++) {
		CpPos2[i] = Vector3d::Random() + dir;
		CpVel2[i] = Vector3d::Random() - dir;
	}
}
static void randomEETest(const BoundingBoxType & bb, const double& deltaDist, const int& kase, const std::string& outputFile){
	Edge obj1, obj2, vel1, vel2;
	std::srand(0);
	int hasCol = 0;
	double deltaT = 1.,t, u1, u2;
	// std::ofstream file(outputFile+".txt");
	// file << std::fixed << std::setprecision(10);

	// obj1.ctrlp ={Vector3d(1,0,0), Vector3d(-1,0,0)},
	// obj2.ctrlp ={Vector3d(0,1,2), Vector3d(0,-1,2)},
	// vel1.ctrlp ={Vector3d(0,0,1), Vector3d(0,0,1)},
	// vel2.ctrlp ={Vector3d(0,0,-1), Vector3d(0,0,-1)};

	using steady_clock = std::chrono::steady_clock;
	using duration = std::chrono::duration<double>;
	const auto initialTime = steady_clock::now();
	for(int k = 0; k < kase; k ++){
		generateEEPair(obj1.ctrlp, vel1.ctrlp, obj2.ctrlp, vel2.ctrlp);
		t = LinearSolverTD::solveEETest(obj1,vel1,obj2,vel2,u1,u2,bb,deltaT,deltaDist);
		saveDoFs<Edge, Edge>(obj1.ctrlp, vel1.ctrlp, obj2.ctrlp, vel2.ctrlp);
		std::cout<<u1<<" "<<u2<<" time "<<t<<"\n";
		if(t>=0){
			hasCol++;
			auto colPos1 = (obj1.lerp(u1)+t*vel1.lerp(u1));
			auto colPos2 = (obj2.lerp(u2)+t*vel2.lerp(u2));
			std::cout<<"case "<<k<<" done, seperation: "<<(colPos1-colPos2).norm()<<"\n";
		}
	}
	const auto endTime = steady_clock::now();
	std::cout << hasCol<<" pairs have collided.\n";
	std::cout << "used seconds: " <<
		duration(endTime - initialTime).count()/kase
		<< std::endl;
	// file.close();
}

static void randomVFTest(const BoundingBoxType & bb, const double& deltaDist, const int& kase, const std::string& outputFile){
	Vector3d obj1, vel1;
	Face obj2, vel2;
	std::srand(0);
	int hasCol = 0;
	double deltaT = 1., t;
	Array2d uv;
	// std::ofstream file(outputFile+".txt");
	// file << std::fixed << std::setprecision(10);

	// obj1 =Vector3d(0,0,0),
	// obj2.ctrlp ={Vector3d(0,1,1), Vector3d(0,-1,1), Vector3d(1,0,1)},
	// vel1 =Vector3d(0,0,1),
	// vel2.ctrlp ={Vector3d(0,0,-1), Vector3d(0,0,-1), Vector3d(0,0,-1)};

	using steady_clock = std::chrono::steady_clock;
	using duration = std::chrono::duration<double>;
	const auto initialTime = steady_clock::now();
	for(int k = 0; k < kase; k ++){
		generateVFPair(obj1, vel1, obj2.ctrlp, vel2.ctrlp);
		t = LinearSolverTD::solveVFTest(obj1,vel1,obj2,vel2,uv,bb,deltaT,deltaDist);
		if(t>=0){
			hasCol++;
			auto colPos1 = (obj1+t*vel1);
			auto colPos2 = (obj2.evaluatePatchPoint(uv)+t*vel2.evaluatePatchPoint(uv));
			std::cout<<"case "<<k<<" done, seperation: "<<(colPos1-colPos2).norm()<<"\n";
		}
	}
	const auto endTime = steady_clock::now();
	std::cout << hasCol<<" pairs have collided.\n";
	std::cout << "used seconds: " <<
		duration(endTime - initialTime).count()/kase
		<< std::endl;
	// file.close();
}

inline std::unique_ptr<ArgsParser> BuildArgsParser()
{
	auto parser = std::make_unique<ArgsParser>();
	parser->addArgument<std::string>("experiment", 'e', "type of experiment (randee, randvf)", "randee");
	parser->addArgument<std::string>("bb", 'b', "type of bounding box (aabb, obb)", "obb");
	parser->addArgument<double>("delta", 'd', "distance for convergence criterion", 1e-6);
	parser->addArgument<int>("kase", 'k', "number of generated cases", 100);
	parser->addArgument<std::string>("output", 'o', "the output filename", "output");
	return parser;
}

int main(int argc, char *argv[]){
	auto parser = BuildArgsParser();
	parser->parse(argc, argv);

    const auto expType = std::any_cast<std::string>(parser->getValueByName("experiment"));	
    const auto bbType = std::any_cast<std::string>(parser->getValueByName("bb"));
    const auto deltaDist = std::any_cast<double>(parser->getValueByName("delta"));
    const auto kase = std::any_cast<int>(parser->getValueByName("kase"));
    const auto outputFile = std::any_cast<std::string>(parser->getValueByName("output"));

	BoundingBoxType bb;
	if(bbType=="obb")
		bb = BoundingBoxType::OBB;
	else if(bbType=="aabb")
		bb = BoundingBoxType::AABB;
	else{
		std::cerr<<"what bounding box?\n";
		exit(-1);
	}

	if(expType=="randee")
		randomEETest(bb, deltaDist, kase, outputFile);
	else if(expType=="randvf")
		randomVFTest(bb, deltaDist, kase, outputFile);
	else{
		std::cerr<<"what experiment?\n";
		exit(-1);
	}
}