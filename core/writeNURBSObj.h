#pragma once
#include "recRationalBezier.h"
#include "config.h"
#include <vector>
#include <fstream>
using Eigen::Vector3i;

class RecRatBezierMesh{
public:
	int cntPatches;
	std::vector<RecQuadRatBezier> patches;

	// 用于初始化速度
	RecRatBezierMesh(){}
	RecRatBezierMesh(const int cnt, const std::array<double, 9>& weight){
		cntPatches = cnt;
		for(int i = 0; i < cnt; i++){
			patches.emplace_back(weight);
		}
	};
	// 读入模型
	// RecRatBezierMesh(const std::string& filename){
	// 	std::ifstream in(filename);
	// 	in>>cntPatches;
	// 	patches.resize(cntPatches);
	// 	int uOrder, vOrder;//暂时没管
	// 	for(auto& patch: patches){
	// 		in>>uOrder>>vOrder;
	// 		for(auto& pt:patch.ctrlp)
	// 			in>>pt[0]>>pt[1]>>pt[2];
	// 	}
	// 	in.close();
	// };

	// int cntPatches() { return cntPatches; }
	void convert2Mesh(std::vector<Vector3d>&verts, std::vector<Vector3i>&faces) const {
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
	void writeObj(const std::string& filename) const {
		std::vector<Vector3d> verts;
		std::vector<Vector3i> faces;
		convert2Mesh(verts, faces);
		std::ofstream out(filename);
		for(auto&vert:verts)
			out<<"v "<<vert[0]<<" "<<vert[1]<<" "<<vert[2]<<"\n";
		for(auto&face:faces)
			out<<"f "<<face[0]<<" "<<face[1]<<" "<<face[2]<<"\n";
		out.close();
	}
	void moveObj(const Vector3d& dis){
		for(auto& patch: patches)
			for(auto& pt:patch.ctrlp){
				Vector3d pos(pt[0]/pt[3], pt[1]/pt[3], pt[2]/pt[3]);
				pt.segment(0,3)=((pos+dis)*pt[3]).eval();
			}
	}
	void rotateObj(const double& angle, const Vector3d& axis, const Vector3d& origin = Vector3d::Zero()){
		//似乎带仿射
		Eigen::AngleAxisd rot(angle, axis);
		for(auto& patch: patches)
			for(auto& pt:patch.ctrlp){
				// pt.segment<3>(0) = (rot.matrix()*(pt-origin)+origin).eval();
				Vector3d pos(pt[0]/pt[3], pt[1]/pt[3], pt[2]/pt[3]);
				pt.segment(0,3)=((rot.matrix()*(pos-origin)+origin)*pt[3]).eval();
			}
	}
};


void generateTorusComponent(RecRatBezierMesh& torus, RecRatBezierMesh& vel){
	double r=2, a=1;
	//顶点位置没有乘权重
	std::array<double, 9> weight = {1,1,2,1,1,2,2,2,4};
	std::array<Vector3d, 9> pos1 = { Vector3d(r+a,0,0), Vector3d(r+a,0,a), Vector3d(r,0,a), 
				Vector3d(r+a,r+a,0), Vector3d(r+a,r+a,a), Vector3d(r,r,a), 
				Vector3d(0,r+a,0), Vector3d(0,r+a,a), Vector3d(0,r,a)	},

	pos2 = { Vector3d(r-a,0,0), Vector3d(r-a,0,a), Vector3d(r,0,a), 
				Vector3d(r-a,r-a,0), Vector3d(r-a,r-a,a), Vector3d(r,r,a), 
				Vector3d(0,r-a,0), Vector3d(0,r-a,a), Vector3d(0,r,a)	},

	pos3 = { Vector3d(r-a,0,0), Vector3d(r-a,0,-a), Vector3d(r,0,-a), 
				Vector3d(r-a,r-a,0), Vector3d(r-a,r-a,-a), Vector3d(r,r,-a), 
				Vector3d(0,r-a,0), Vector3d(0,r-a,-a), Vector3d(0,r,-a)	},

	pos4 = { Vector3d(r+a,0,0), Vector3d(r+a,0,-a), Vector3d(r,0,-a), 
				Vector3d(r+a,r+a,0), Vector3d(r+a,r+a,-a), Vector3d(r,r,-a), 
				Vector3d(0,r+a,0), Vector3d(0,r+a,-a), Vector3d(0,r,-a)	};

	torus.patches.clear();
	torus.patches.emplace_back(pos1,weight);
	torus.patches.emplace_back(pos2,weight);
	torus.patches.emplace_back(pos3,weight);
	torus.patches.emplace_back(pos4,weight);
	RecRatBezierMesh torus2=torus;
	for(int i=0;i<3;i++){
		torus2.rotateObj(PI*0.5, Vector3d::Unit(2));
		torus.patches.insert(torus.patches.end(),torus2.patches.begin(),torus2.patches.end());
	}
	torus.cntPatches = 16;
	vel = RecRatBezierMesh(16, weight);
}