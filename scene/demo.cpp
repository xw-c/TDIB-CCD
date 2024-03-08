#include "scene.h"
#include "argsParser.h"
inline std::unique_ptr<ArgsParser> BuildArgsParser()
{
	auto parser = std::make_unique<ArgsParser>();
	parser->addArgument<std::string>("solver", 's', "type of ccd solver (base, td)", "base");
	// parser->addArgument<BoundingBoxType>("bb", 'b', "type of bounding box", BBDefault);
	// parser->addArgument<bool>("process", 'p', "show details of solving process", SHOWANS);
	parser->addArgument<double>("delta", 'd', "distance for convergence criterion", MinL1Dist);
	parser->addArgument<int>("kase", 'k', "number of generated cases", KaseDefault);
	parser->addArgument<double>("velocity", 'v', "magnitude of velocity", PullVelocity);
	parser->addArgument<std::string>("output", 'o', "the output filename", "output");
	// parser->addArgument<bool>("collision", 'p', "use self-collision of the patches or not", false);
	return parser;
}



int main(int argc, char *argv[]){
	auto parser = BuildArgsParser();
	parser->parse(argc, argv);

    const auto solverType = std::any_cast<std::string>(parser->getValueByName("solver"));	
    const auto deltaDist = std::any_cast<double>(parser->getValueByName("delta"));
    const auto kase = std::any_cast<int>(parser->getValueByName("kase"));
    const auto velMag = std::any_cast<double>(parser->getValueByName("velocity"));
    const auto outputFile = std::any_cast<std::string>(parser->getValueByName("output"));

	// planeTest<RecCubicBezier, RecParamBound>(solverType, deltaDist, kase, velMag, outputFile);
	randomTest<RecCubicBezier, RecParamBound>(solverType, deltaDist, kase, velMag, outputFile);

	// parabolaPotCup();
	// parabolaBunnyTorus();
	// torusTest();
	// randomTest<TriLinearBezier, TriParamBound>(3);
	// randomTest<TriQuadBezier, TriParamBound>(2);
	// randomTest<TriCubicBezier, TriParamBound>(1.625);
	// randomTest<RecLinearBezier, RecParamBound>(3);
	// randomTest<RecQuadBezier, RecParamBound>(2);
	// randomTest<RecCubicBezier, RecParamBound>(1.625);
	// randomTest<RecQuadRatBezier, RecParamBound>(2);
	// randomTest<RecCubicRatBezier, RecParamBound>(1.625);
}