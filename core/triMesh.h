#pragma once
#include "triBezier.h"
#include <vector>
#include <fstream>
using Eigen::Vector3i;

class TriLinearMesh{
	public:
	int cntPatches;
	// int cntVerts, cntFaces;
	std::vector<TriLinearBezier> patches;
	// std::vector<Vector3d> verts;
	// std::vector<Vector3i> faces;

	// 用于初始化速度
	TriLinearMesh(const int _cntFaces): cntPatches(_cntFaces)
	// , cntVerts(_cntVerts), cntFaces(_cntFaces)
	{
		patches.resize(cntPatches);
		for(auto& patch: patches){
			for(auto& pt:patch.ctrlp)
				pt.setZero();
		}
	};
	// 读入模型
	TriLinearMesh(const std::string& filename){
		std::ifstream in(filename);
		int cntV;
		in >> cntV >> cntPatches;

		std::vector<Vector3d> verts;
		std::vector<Vector3i> faces;
		verts.resize(cntV);
		faces.resize(cntPatches);
		char c;
		for(auto& vert: verts)
			in >> c >> vert[0] >> vert[1] >> vert[2];
		for(auto& face: faces){
			in >> c >> face[0] >> face[1] >> face[2];
			std::array<Vector3d, 3> pts = {verts[face[0]-1], verts[face[1]-1], verts[face[2]-1]};
			patches.emplace_back(pts);
		}
		in.close();
	};
	void convert2Mesh(std::vector<Vector3d>&verts, std::vector<Vector3i>&faces) const {
		int cntVerts=0;
		for(auto& patch: patches){
			for(auto& pt: patch.ctrlp)
				verts.push_back(pt);
			faces.emplace_back(cntVerts+1, cntVerts+2, cntVerts+3);
			cntVerts+=3;
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
	void rotateObj(const double& angle, const Vector3d& axis, const Vector3d& origin = Vector3d::Zero()){
		//似乎带仿射
		Eigen::AngleAxisd rot(angle, axis);
		for(auto& patch: patches)
			for(auto& pt:patch.ctrlp)
				pt = (rot.matrix()*(pt-origin)+origin).eval();
	}
};