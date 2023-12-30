#include "paramCCD.h"
#include "sampleCCD.h"
#include "config.h"

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
	// certificate();
	TriBezierObj p1,p2,v1,v2;
	// generatePatches(p1,p2,v1,v2);
	// std::cout<<recBezierCCD(p1,p2,v1,v2, uv1,uv2,BoundingBoxType::OBB);
	RecBezierObj p;
	p.ctrlp={
		Vector3d(0, 0, 0), Vector3d(1, 0, 0), Vector3d(2, 0, 0), Vector3d(3, 0, 0.1),
		Vector3d(0, 1, 0), Vector3d(1, 1, 0), Vector3d(2, 1, 0), Vector3d(3, 1, 0.1),
		Vector3d(0, 2, 0), Vector3d(1, 2, 0), Vector3d(2, 2, 0), Vector3d(3, 2, 0.1),
		Vector3d(0, 3, 0), Vector3d(1, 3, 0), Vector3d(2, 3, 0), Vector3d(3, 3, 0.1),
	};
	std::cout<<p.evaluatePatchPoint(Array2d(0,1));
}