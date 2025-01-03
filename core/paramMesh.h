#pragma once
#include <vector>
#include <fstream>
#include <Eigen/Dense>
using Eigen::Array2d;
using Eigen::Vector3d;
using Eigen::Vector4d;
using Eigen::Vector3i;

template<typename PatchType>
class MeshBase{
public:
	int cntPatches;
	std::vector<PatchType> patches;
	
	Vector3d evaluatePatchNormal(const PatchType& patch, Array2d const &uv) const {
		double diff = 0.001;
		Vector3d p = patch.evaluatePatchPoint(uv),
		pdu = patch.evaluatePatchPoint(Array2d(uv[0]+diff, uv[1])),
		pdv = patch.evaluatePatchPoint(Array2d(uv[0], uv[1]+diff));
		return ((pdu-p).cross(pdv-p)).normalized();
	}

	// This function can sometimes affected by floating-point errors, resulting in extra or leaving patches
	void convert2Mesh(std::vector<Vector3d>&verts, std::vector<Vector3d>&norms, std::vector<Vector3i>&faces, 
					const double du=0.2, const double dv=0.2) const {
		int cntVerts=0;
		for(auto& patch: patches){
			for(double u=0;u<1-1e-10;u+=du)
				for(double v=0;v<patch.feasibleUpperV(u)-1e-10;v+=dv){
					Vector3d p1=patch.evaluatePatchPoint(Array2d(u,v)),
					p2=patch.evaluatePatchPoint(Array2d(u+du,v)),
					p3=patch.evaluatePatchPoint(Array2d(u,v+dv));
					Vector3d n1=evaluatePatchNormal(patch, Array2d(u,v)),
					n2=evaluatePatchNormal(patch, Array2d(u+du,v)),
					n3=evaluatePatchNormal(patch, Array2d(u,v+dv));
					verts.push_back(p1);
					verts.push_back(p2);
					verts.push_back(p3);
					norms.push_back(n1);
					norms.push_back(n2);
					norms.push_back(n3);
					faces.emplace_back(cntVerts+1, cntVerts+2, cntVerts+3);
					if(v<patch.feasibleUpperV(u+du)-1e-10){
						Vector3d p4=patch.evaluatePatchPoint(Array2d(u+du,v+dv));
						Vector3d n4=evaluatePatchNormal(patch, Array2d(u+du,v+dv));
						verts.push_back(p4);
						norms.push_back(n4);
						faces.emplace_back(cntVerts+4, cntVerts+3, cntVerts+2);
						cntVerts+=4;
					}
					else cntVerts+=3;
				}
		}
	}

	void writeObj(const std::string& filename, const double du=0.2, const double dv=0.2) const {
		std::vector<Vector3d> verts;
		std::vector<Vector3d> norms;
		std::vector<Vector3i> faces;
		convert2Mesh(verts, norms, faces, du, dv);
		std::ofstream out(filename);
		// out<<std::fixed<<std::setprecision(10);
		for(auto&vert:verts)
			out<<"v "<<vert[0]<<" "<<vert[1]<<" "<<vert[2]<<"\n";
		for(auto&norm:norms)
			out<<"vn "<<norm[0]<<" "<<norm[1]<<" "<<norm[2]<<"\n";
		for(auto&face:faces)
			out<<"f "<<face[0]<<"//"<<face[0]<<" "
				<<face[1]<<"//"<<face[1]<<" "
				<<face[2]<<"//"<<face[2]<<"\n";
		out.close();
	}
};


template<typename PatchType>
class ParamMesh: public MeshBase<PatchType>{
public:
    using MeshBase<PatchType>::cntPatches;
    using MeshBase<PatchType>::patches;

	ParamMesh(){}

	// Initialize zero velocities
	ParamMesh(const int cnt){
		cntPatches = cnt;
		patches.resize(cntPatches);
		for(auto& patch: patches){
			for(auto& pt:patch.ctrlp)
				pt.setZero();
		}
	};

	void moveObj(const Vector3d& dis){
		for(auto& patch: patches)
			for(auto& pt:patch.ctrlp)
				pt+=dis;
	}
	void moveObj(const ParamMesh& vel, const double& t){
		for(int i=0;i<cntPatches;i++)
			for(int j=0;j<PatchType::cntCp;j++)
				patches[i].ctrlp[j]+=vel.patches[i].ctrlp[j]*t;
	}
	Vector3d getOrigin(){
		Vector3d center = Vector3d::Zero();
		for(const auto& patch:patches){
			for(const auto& p:patch.ctrlp){
				center+=p;
			}
		}
		center /= cntPatches * PatchType::cntCp;
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
		// "axis" should be unit vector, otherwise obj would be stretched
		Eigen::AngleAxisd rot(angle, axis);
		for(auto& patch: patches)
			for(auto& pt:patch.ctrlp)
				pt = (rot.matrix()*(pt-origin)+origin).eval();
	}
	void setVel(const ParamMesh& p1, const ParamMesh& p2, const double& dt){
		for(int i=0;i<cntPatches;i++)
			for(int j=0;j<PatchType::cntCp;j++)
				patches[i].ctrlp[j] = (p2.patches[i].ctrlp[j] - p1.patches[i].ctrlp[j]) / dt;
	}
};

template<typename PatchType>
class RatParamMesh: public MeshBase<PatchType>{
public:
    using MeshBase<PatchType>::cntPatches;
    using MeshBase<PatchType>::patches;

	RatParamMesh(){}

	// Initialize zero velocities
	RatParamMesh(const int cnt, const std::array<double, PatchType::cntCp>& weight){
		cntPatches = cnt;
		for(int i = 0; i < cnt; i++){
			patches.emplace_back(weight);
		}
	}

	Vector3d getOrigin(){
		Vector3d center = Vector3d::Zero();
		for(const auto& patch:patches){
			for(const auto& p:patch.ctrlp){
				center+=p.segment(0,3)/p[3];
			}
		}
		center /= cntPatches * PatchType::cntCp;
		return center;
	}
	void setOrigin(const Vector3d& dis){
		Vector3d diff = dis - getOrigin();
		moveObj(diff);
	}
	void moveObj(const Vector3d& dis){
		for(auto& patch: patches)
			for(auto& pt:patch.ctrlp){
				Vector3d pos=pt.segment(0,3)/pt[3];
				pt.segment(0,3)=((pos+dis)*pt[3]).eval();
			}
	}
	void moveObj(const RatParamMesh<PatchType>& vel, const double& t){
		for(int i=0;i<cntPatches;i++)
			for(int j=0;j<PatchType::cntCp;j++){
				// 前提是vel的Weight与pos完全一致
				patches[i].ctrlp[j].segment(0,3)+=vel.patches[i].ctrlp[j].segment(0,3)*t;
			}
	}
	void rotateObj(const double& angle, const Vector3d& axis, const Vector3d& origin = Vector3d::Zero()){
		// "axis" should be unit vector, otherwise obj would be stretched
		Eigen::AngleAxisd rot(angle, axis);
		for(auto& patch: patches)
			for(auto& pt:patch.ctrlp){
				Vector3d pos(pt[0]/pt[3], pt[1]/pt[3], pt[2]/pt[3]);
				pt.segment(0,3)=((rot.matrix()*(pos-origin)+origin)*pt[3]).eval();
			}
	}
	void setVel(const RatParamMesh<PatchType>& p1, const RatParamMesh<PatchType>& p2, const double& dt){
		for(int i=0;i<cntPatches;i++)
			for(int j=0;j<PatchType::cntCp;j++)
				patches[i].ctrlp[j].segment(0,3) = 
				(p2.patches[i].ctrlp[j].segment(0,3) - p1.patches[i].ctrlp[j].segment(0,3)) / dt;
	}
};