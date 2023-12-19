#include "triBezier.h"
#include "paramCCD.h"
#include "config.h"

void saveDoFs(){
	std::ofstream f("DoFs.txt");
	for(auto item : CpPos1.ctrlp)
		f<<item.transpose()<<"\n";
	for(auto item : CpPos2.ctrlp)
		f<<item.transpose()<<"\n";
	for(auto item : CpVel1.ctrlp)
		f<<item.transpose()<<"\n";
	for(auto item : CpVel2.ctrlp)
		f<<item.transpose()<<"\n";
	f.close();
}

void randomTest(){
	using steady_clock = std::chrono::steady_clock;
	using duration = std::chrono::duration<double>;
	int cntAABB=0, cntSample=0;
	const int Kase = 100;
	double ans[2][Kase];

	const auto initOBB = steady_clock::now();
	std::srand(0);
	for(int kase=0;kase<Kase;kase++){
		generateTriBeizer();
		ans[0][kase]=ccd(BB::OBB);
	}
	const auto endOBB = steady_clock::now();
	std::cout<<"OBB used seconds: "<<duration(endOBB - initOBB).count()/Kase<<"\n";

	const auto initAABB = steady_clock::now();
	// std::srand(0);
	for(int kase=0;kase<Kase;kase++){
		generateTriBeizer();
		ans[0][kase]=ccd(BB::AABB);
	}
	const auto endAABB = steady_clock::now();
	std::cout<<"AABB used seconds: "<<duration(endAABB - initAABB).count()/Kase<<"\n";
	std::cout<<"OBB used seconds: "<<duration(endOBB - initOBB).count()/Kase<<"\n";
}
void certificate(){
	std::srand(std::time(nullptr));
    generateTriBeizer();

    {
        double t = ccd(BB::OBB);

        Vector3d const p1 = CpPos1.blossomBicubicBezier(uv1);
        Vector3d const v1 = CpVel1.blossomBicubicBezier(uv1);
        Vector3d const p2 = CpPos2.blossomBicubicBezier(uv2);
        Vector3d const v2 = CpVel2.blossomBicubicBezier(uv2);
        Vector3d const pt1=(v1*t+p1), pt2=(v2*t+p2);
        std::cout<<"pos at: "<<pt1.transpose()<<"     "<<pt2.transpose()<<"\n";
        std::cout<<"delta: "<<(pt2-pt1).norm()<<"\n";
    }

    // {
    //     double t = ccd(BB::AABB);

    //     Vector3d const p1 = CpPos1.blossomBicubicBezier(uv1);
    //     Vector3d const v1 = CpVel1.blossomBicubicBezier(uv1);
    //     Vector3d const p2 = CpPos2.blossomBicubicBezier(uv2);
    //     Vector3d const v2 = CpVel2.blossomBicubicBezier(uv2);
    //     Vector3d const pt1=(v1*t+p1), pt2=(v2*t+p2);
    //     std::cout<<"pos at: "<<pt1.transpose()<<"     "<<pt2.transpose()<<"\n";
    //     std::cout<<"delta: "<<(pt2-pt1).norm()<<"\n";
    // }
}
int main(){
	ParamCCD<PatchType::TriBezier> ccdSolver(BoundingBoxType::OBB);
	ccdSolver.CCD();
}