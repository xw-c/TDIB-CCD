#pragma once
#include "recBezier.h"
#include <vector>
#include <fstream>
using Eigen::Vector3i;

class RecBezierMesh{
public:
	int cntPatches;
	std::vector<RecCubicBezier> patches;

	// 用于初始化速度
	RecBezierMesh(const int cnt){
		cntPatches = cnt;
		patches.resize(cntPatches);
		for(auto& patch: patches){
			for(auto& pt:patch.ctrlp)
				pt.setZero();
		}
	};
	// 读入模型
	RecBezierMesh(const std::string& filename){
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
	};

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
			for(auto& pt:patch.ctrlp)
				pt+=dis;
	}
	void moveObj(const RecBezierMesh& vel, const double& t){
		for(int i=0;i<cntPatches;i++)
			for(int j=0;j<RecCubicBezier::cntCp;j++)
				patches[i].ctrlp[j]+=vel.patches[i].ctrlp[j]*t;
	}
	Vector3d getOrigin(){
		Vector3d center = Vector3d::Zero();
		for(const auto& patch:patches){
			for(const auto& p:patch.ctrlp){
				center+=p;
			}
		}
		center /= cntPatches * RecCubicBezier::cntCp;
		return center;
	}
	void setOrigin(const Vector3d& dis){
		Vector3d diff = getOrigin() - dis;
		for(auto& patch:patches){
			for(auto& p:patch.ctrlp){
				p -= diff;
			}
		}
	}
	void rotateObj(const double& angle, const Vector3d& axis, const Vector3d& origin = Vector3d::Zero()){
		//似乎带仿射
		Eigen::AngleAxisd rot(angle, axis);
		for(auto& patch: patches)
			for(auto& pt:patch.ctrlp)
				pt = (rot.matrix()*(pt-origin)+origin).eval();
	}
	void setVel(const RecBezierMesh& p1, const RecBezierMesh& p2, const double& dt){
		for(int i=0;i<cntPatches;i++)
			for(int j=0;j<RecCubicBezier::cntCp;j++)
				patches[i].ctrlp[j] = (p2.patches[i].ctrlp[j] - p1.patches[i].ctrlp[j]) / dt;
	}
};