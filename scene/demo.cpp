#include "scene.h"
int main(){
	// parabolaPotCup();
	// parabolaBunnyTorus();
	// randomTest<TriLinearBezier, TriParamBound>(3);
	// randomTest<TriQuadBezier, TriParamBound>(2);
	// randomTest<TriCubicBezier, TriParamBound>(1.625);
	// randomTest<RecLinearBezier, RecParamBound>(3);
	// randomTest<RecQuadBezier, RecParamBound>(2);
	randomTest<RecCubicBezier, RecParamBound>(1.625);
	// randomTest<RecQuadRatBezier, RecParamBound>(2);
	// randomTest<RecCubicRatBezier, RecParamBound>(1.625);
}