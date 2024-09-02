#pragma once
#include "config.h"
template<typename ObjType1, typename ObjType2>
static void readinDoFs(std::array<Vector3d, ObjType1::cntCp>& CpPos1, 
					std::array<Vector3d, ObjType1::cntCp>& CpVel1,
					std::array<Vector3d, ObjType2::cntCp>& CpPos2, 
					std::array<Vector3d, ObjType2::cntCp>& CpVel2,
					const std::string& filename="DoFs.dat"){
	std::ifstream readin(filename, std::ios::binary);
	for (int i = 0; i < ObjType1::cntCp; i++)
		for(int k=0;k<3;k++)
			readin.read(reinterpret_cast<char *>(&CpPos1[i][k]), sizeof(double));
	for (int i = 0; i < ObjType1::cntCp; i++)
		for(int k=0;k<3;k++)
			readin.read(reinterpret_cast<char *>(&CpVel1[i][k]), sizeof(double));
	for (int i = 0; i < ObjType2::cntCp; i++)
		for(int k=0;k<3;k++)
			readin.read(reinterpret_cast<char *>(&CpPos2[i][k]), sizeof(double));
	for (int i = 0; i < ObjType2::cntCp; i++)
		for(int k=0;k<3;k++)
			readin.read(reinterpret_cast<char *>(&CpVel2[i][k]), sizeof(double));
	readin.close();
}
template<typename ObjType1, typename ObjType2>
static void saveDoFs(const std::array<Vector3d, ObjType1::cntCp>& CpPos1, 
					const std::array<Vector3d, ObjType1::cntCp>& CpVel1,
					const std::array<Vector3d, ObjType2::cntCp>& CpPos2, 
					const std::array<Vector3d, ObjType2::cntCp>& CpVel2,
					const std::string& filename="DoFs.dat"){
	std::ofstream outfile("filename.dat", std::ios::out|std::ios::trunc|std::ios::binary);
	for(int pid = 0; pid < ObjType1::cntCp; ++pid)
		for(int dim=0; dim < 3; ++dim) 
			outfile.write(reinterpret_cast<const char*>(&CpPos1[pid][dim]),sizeof(double));
	for(int pid = 0; pid < ObjType1::cntCp; ++pid)
		for(int dim=0; dim < 3; ++dim) 
			outfile.write(reinterpret_cast<const char*>(&CpVel1[pid][dim]),sizeof(double));
	for(int pid = 0; pid < ObjType2::cntCp; ++pid)
		for(int dim=0; dim < 3; ++dim) 
			outfile.write(reinterpret_cast<const char*>(&CpPos2[pid][dim]),sizeof(double));
	for(int pid = 0; pid < ObjType2::cntCp; ++pid)
		for(int dim=0; dim < 3; ++dim) 
			outfile.write(reinterpret_cast<const char*>(&CpVel2[pid][dim]),sizeof(double));
	outfile.close();
}

template<typename ObjType>
static void generatePatchPair(std::array<Vector3d, ObjType::cntCp> &CpPos1, std::array<Vector3d, ObjType::cntCp> &CpVel1,
									std::array<Vector3d, ObjType::cntCp> &CpPos2, std::array<Vector3d, ObjType::cntCp> &CpVel2, const double& velMagnitude = 1){
	Vector3d dir=Vector3d::Random().normalized();
	for (int i = 0; i < ObjType::cntCp; i++) {
		CpPos1[i] = Vector3d::Random() - dir;
		CpVel1[i] = Vector3d::Random() + dir * velMagnitude;
		CpPos2[i] = Vector3d::Random() + dir;
		CpVel2[i] = Vector3d::Random() - dir * velMagnitude;
	}
}
template<typename ObjType>
static void generatePatchPair(std::array<Vector4d, ObjType::cntCp> &CpPos1, std::array<Vector4d, ObjType::cntCp> &CpVel1,
									std::array<Vector4d, ObjType::cntCp> &CpPos2, std::array<Vector4d, ObjType::cntCp> &CpVel2, const double& velMagnitude = 1){
	Vector3d dir=Vector3d::Random().normalized();
	for (int i = 0; i < ObjType::cntCp; i++) {
		CpPos1[i][3]=((double)rand()+1)/(RAND_MAX);
		CpVel1[i][3]=CpPos1[i][3];
		CpPos2[i][3]=((double)rand()+1)/(RAND_MAX);
		CpVel2[i][3]=CpPos2[i][3];

		CpPos1[i].segment(0,3) = CpPos1[i][3] * (Vector3d::Random() - dir);
		CpVel1[i].segment(0,3) = CpVel1[i][3] * (Vector3d::Random() + dir * velMagnitude);
		CpPos2[i].segment(0,3) = CpPos2[i][3] * (Vector3d::Random() + dir);
		CpVel2[i].segment(0,3) = CpVel2[i][3] * (Vector3d::Random() - dir * velMagnitude);
	}

}