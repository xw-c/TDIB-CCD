#include "paramCCD.h"
#include "sampleCCD.h"
#include "config.h"

void saveDoFs(){

	// std::ofstream f("DoFs.txt");
	// for(auto item : cpp1.ctrlp)
	// 	f<<item.transpose()<<"\n";
	// for(auto item : cpv1.ctrlp)
	// 	f<<item.transpose()<<"\n";
	// for(auto item : cpp2.ctrlp)
	// 	f<<item.transpose()<<"\n";
	// for(auto item : cpv2.ctrlp)
	// 	f<<item.transpose()<<"\n";
	// f.close();
}

// void randomTest(){
// 	using steady_clock = std::chrono::steady_clock;
// 	using duration = std::chrono::duration<double>;
// 	int cntAABB=0, cntSample=0;
// 	const int Kase = 100;
// 	double ans[2][Kase];

// 	const auto initOBB = steady_clock::now();
// 	std::srand(0);
// 	for(int kase=0;kase<Kase;kase++){
// 		generateTriBeizer();
// 		ans[0][kase]=ccd(BB::OBB);
// 	}
// 	const auto endOBB = steady_clock::now();
// 	std::cout<<"OBB used seconds: "<<duration(endOBB - initOBB).count()/Kase<<"\n";

// 	const auto initAABB = steady_clock::now();
// 	// std::srand(0);
// 	for(int kase=0;kase<Kase;kase++){
// 		generateTriBeizer();
// 		ans[0][kase]=ccd(BB::AABB);
// 	}
// 	const auto endAABB = steady_clock::now();
// 	std::cout<<"AABB used seconds: "<<duration(endAABB - initAABB).count()/Kase<<"\n";
// 	std::cout<<"OBB used seconds: "<<duration(endOBB - initOBB).count()/Kase<<"\n";
// }
void certificate(){
	std::srand(std::time(nullptr));
	RecBezierObj cpp1,cpp2,cpv1,cpv2;
	generatePatches(cpp1,cpp2,cpv1,cpv2);
	Array2d uv1,uv2;
    {
        double t = recBezierCCD(cpp1,cpv1,cpp2,cpv2, uv1,uv2,BoundingBoxType::OBB);

        Vector3d const p1 = cpp1.evaluatePatchPoint(uv1);
        Vector3d const v1 = cpv1.evaluatePatchPoint(uv1);
        Vector3d const p2 = cpp2.evaluatePatchPoint(uv2);
        Vector3d const v2 = cpv2.evaluatePatchPoint(uv2);
        Vector3d const pt1=(v1*t+p1), pt2=(v2*t+p2);
        std::cout<<"pos at: "<<pt1.transpose()<<"     "<<pt2.transpose()<<"\n";
        std::cout<<"delta: "<<(pt2-pt1).norm()<<"\n";
    }

    {
        double t = recBezierCCD(cpp1,cpv1,cpp2,cpv2, uv1,uv2,BoundingBoxType::AABB);

        Vector3d const p1 = cpp1.evaluatePatchPoint(uv1);
        Vector3d const v1 = cpv1.evaluatePatchPoint(uv1);
        Vector3d const p2 = cpp2.evaluatePatchPoint(uv2);
        Vector3d const v2 = cpv2.evaluatePatchPoint(uv2);
        Vector3d const pt1=(v1*t+p1), pt2=(v2*t+p2);
        std::cout<<"pos at: "<<pt1.transpose()<<"     "<<pt2.transpose()<<"\n";
        std::cout<<"delta: "<<(pt2-pt1).norm()<<"\n";
    }
	{
        double t = recSampleCCD(cpp1,cpv1,cpp2,cpv2, uv1,uv2,BoundingBoxType::DOP14);

        Vector3d const p1 = cpp1.evaluatePatchPoint(uv1);
        Vector3d const v1 = cpv1.evaluatePatchPoint(uv1);
        Vector3d const p2 = cpp2.evaluatePatchPoint(uv2);
        Vector3d const v2 = cpv2.evaluatePatchPoint(uv2);
        Vector3d const pt1=(v1*t+p1), pt2=(v2*t+p2);
        std::cout<<"pos at: "<<pt1.transpose()<<"     "<<pt2.transpose()<<"\n";
        std::cout<<"delta: "<<(pt2-pt1).norm()<<"\n";
    }
}
int main(){
	certificate();
	// RecBezierObj p1,p2,v1,v2;
	// generatePatches(p1,p2,v1,v2);
	// std::cout<<recBezierCCD(p1,p2,v1,v2, uv1,uv2,BoundingBoxType::OBB);
}