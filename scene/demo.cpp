#include "scene_random.h"
#include "scene_examples.h"
#include "argsParser.h"
inline std::unique_ptr<ArgsParser> BuildArgsParser()
{
	auto parser = std::make_unique<ArgsParser>();
	parser->addArgument<std::string>("solver", 's', "type of ccd solver (trad, td)", "td");
	parser->addArgument<std::string>("experiment", 'e', "type of experiment (rand, single, bunny)", "rand");
	parser->addArgument<std::string>("bb", 'b', "type of bounding box (aabb, obb)", "obb");

	parser->addArgument<double>("delta", 'd', "distance for convergence criterion", 1e-6);
	parser->addArgument<int>("kase", 'k', "number of generated cases", 100);
	return parser;
}

int main(int argc, char *argv[]){
	auto parser = BuildArgsParser();
	parser->parse(argc, argv);

    const auto expType = std::any_cast<std::string>(parser->getValueByName("experiment"));	
    const auto solverType = std::any_cast<std::string>(parser->getValueByName("solver"));	
    const auto bbType = std::any_cast<std::string>(parser->getValueByName("bb"));
    const auto deltaDist = std::any_cast<double>(parser->getValueByName("delta"));
    const auto kase = std::any_cast<int>(parser->getValueByName("kase"));

	BoundingBoxType bb;
	if(bbType=="obb")
		bb = BoundingBoxType::OBB;
	else if(bbType=="aabb")
		bb = BoundingBoxType::AABB;
	else{
		std::cerr<<"Bounding box not implemented.\n";
		exit(-1);
	}
	SolverType solver;
	if(solverType=="td")
		solver = SolverType::TDIntv;
	else if(solverType=="trad")
		solver = SolverType::TradIntv;
	else{
		std::cerr<<"Bounding box not implemented.\n";
		exit(-1);
	}

	if(expType=="single")
		singleTest<RecCubicBezier, RecParamBound>(solver, bb, deltaDist);
	else if(expType=="rand")
		randomTest<RecCubicBezier, RecParamBound>(solver, bb, deltaDist, kase);
	else if(expType=="bunny")
		parabolaBunnyTorus();
	else{
		std::cerr<<"Experiment not implemented.\n";
		exit(-1);
	}
}