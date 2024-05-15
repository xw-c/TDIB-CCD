#include "scene.h"
#include "argsParser.h"
inline std::unique_ptr<ArgsParser> BuildArgsParser()
{
	auto parser = std::make_unique<ArgsParser>();
	parser->addArgument<std::string>("solver", 's', "type of ccd solver (base, td)", "td");
	parser->addArgument<std::string>("experiment", 'e', "type of experiment (rand, valid)", "rand");
	parser->addArgument<std::string>("bb", 'b', "type of bounding box (aabb, obb)", "obb");
	// parser->addArgument<std::string>("primitive", 'p', "type of primitive (tri/rec + 1/2/3 + /rat)", "rec3");

	parser->addArgument<double>("delta", 'd', "distance for convergence criterion", 1e-6);
	parser->addArgument<int>("kase", 'k', "number of generated cases", 50);
	parser->addArgument<double>("velocity", 'v', "magnitude of velocity", 1);

	// parser->addArgument<bool>("process", 'p', "show details of solving process", SHOWANS);
	parser->addArgument<std::string>("output", 'o', "the output filename", "output");
	return parser;
}



int main(int argc, char *argv[]){
	auto parser = BuildArgsParser();
	parser->parse(argc, argv);

    const auto solverType = std::any_cast<std::string>(parser->getValueByName("solver"));	
    const auto expType = std::any_cast<std::string>(parser->getValueByName("experiment"));	
    const auto bbType = std::any_cast<std::string>(parser->getValueByName("bb"));
    // const auto primType = std::any_cast<std::string>(parser->getValueByName("primitive"));	
    const auto deltaDist = std::any_cast<double>(parser->getValueByName("delta"));
    const auto kase = std::any_cast<int>(parser->getValueByName("kase"));
    const auto velMag = std::any_cast<double>(parser->getValueByName("velocity"));
    const auto outputFile = std::any_cast<std::string>(parser->getValueByName("output"));

	BoundingBoxType bb;
	if(bbType=="obb")
		bb = BBDefault = BoundingBoxType::OBB;
	else if(bbType=="aabb")
		bb = BBDefault = BoundingBoxType::AABB;
	else{
		std::cerr<<"what bounding box?\n";
		exit(-1);
	}

	if(expType=="valid")
		validate<RecCubicBezier, RecParamBound>(solverType, bb, deltaDist, kase, velMag, outputFile);
	else if(expType=="rand")
		randomTest<RecCubicBezier, RecParamBound>(solverType, bb, deltaDist, kase, velMag, outputFile);
			else if(expType=="fn")
		FNCase(solverType, bb, deltaDist, kase, velMag, outputFile);
	else if(expType=="all"){
		randomTest<TriLinearBezier, TriParamBound>(solverType, bb, deltaDist, kase, velMag, outputFile);
		randomTest<TriQuadBezier, TriParamBound>(solverType, bb, deltaDist, kase, velMag, outputFile);
		randomTest<TriCubicBezier, TriParamBound>(solverType, bb, deltaDist, kase, velMag, outputFile);
		randomTest<RecLinearBezier, RecParamBound>(solverType, bb, deltaDist, kase, velMag, outputFile);
		randomTest<RecQuadBezier, RecParamBound>(solverType, bb, deltaDist, kase, velMag, outputFile);
		randomTest<RecCubicBezier, RecParamBound>(solverType, bb, deltaDist, kase, velMag, outputFile);
		// randomTest<RecQuadRatBezier, RecParamBound>(solverType, bb, deltaDist, kase, velMag, outputFile);
		// randomTest<RecCubicRatBezier, RecParamBound>(solverType, bb, deltaDist, kase, velMag, outputFile);
	}
	else{
		std::cerr<<"what experiment?\n";
		exit(-1);
	}

	// parabolaPotCup();
	// parabolaBunnyTorus();
	// torusTest(solverType, deltaDist, outputFile);
	// randomTest<TriLinearBezier, TriParamBound>(3);
	// randomTest<TriQuadBezier, TriParamBound>(2);
	// randomTest<TriCubicBezier, TriParamBound>(1.625);
	// randomTest<RecLinearBezier, RecParamBound>(3);
	// randomTest<RecQuadBezier, RecParamBound>(2);
	// randomTest<RecCubicBezier, RecParamBound>(1.625);
	// randomTest<RecQuadRatBezier, RecParamBound>(2);
	// randomTest<RecCubicRatBezier, RecParamBound>(1.625);
}