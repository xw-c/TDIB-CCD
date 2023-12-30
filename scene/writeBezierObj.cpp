#include "recBezier.h"
#include <vector>
#include <fstream>
using Eigen::Vector3i;

class BezierMesh{
	int cntPatches;
	std::vector<RecBezierObj> patches;

	std::vector<Vector3d> verts;
	std::vector<Vector3i> faces;

public:
	BezierMesh(const std::string& filename){
		std::ifstream in(filename);
		in>>cntPatches;
		patches.resize(cntPatches);
		int uOrder, vOrder;//暂时没管
		for(auto& patch: patches){
			in>>uOrder>>vOrder;
			for(auto& pt:patch.ctrlp)
				in>>pt[0]>>pt[1]>>pt[2];
		}
		in.close();
		convert2Mesh();
	};
	void convert2Mesh(){
		double du=0.2, dv=0.2;
		int cntVerts=0;
		for(auto& patch: patches){
			for(double u=0;u<1;u+=du)
				for(double v=0;v<1;v+=dv){
					Vector3d p1=patch.evaluatePatchPoint(Array2d(u,v)),
					p2=patch.evaluatePatchPoint(Array2d(u+du,v)),
					p3=patch.evaluatePatchPoint(Array2d(u,v+dv)),
					p4=patch.evaluatePatchPoint(Array2d(u+du,v+dv));
					verts.push_back(p1);
					verts.push_back(p2);
					verts.push_back(p3);
					verts.push_back(p4);
					faces.emplace_back(cntVerts+1, cntVerts+2, cntVerts+3);
					faces.emplace_back(cntVerts+4, cntVerts+3, cntVerts+2);
					cntVerts+=4;
				}
		}
	}
	void writeObj(const std::string& filename){
		std::ofstream out(filename);
		for(auto&vert:verts)
			out<<"v "<<vert[0]<<" "<<vert[1]<<" "<<vert[2]<<"\n";
		for(auto&face:faces)
			out<<"f "<<face[0]<<" "<<face[1]<<" "<<face[2]<<"\n";
		out.close();
	};
};
int main(){
	BezierMesh obj("teapot32.txt");
	obj.writeObj("teapot32.obj");
}