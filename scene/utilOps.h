#pragma once
#include "config.h"
template<typename ObjType1, typename ObjType2>
static void readinDoFs(std::array<Vector3d, ObjType1::cntCp>& CpPos1, 
					std::array<Vector3d, ObjType1::cntCp>& CpVel1,
					std::array<Vector3d, ObjType2::cntCp>& CpPos2, 
					std::array<Vector3d, ObjType2::cntCp>& CpVel2,
					const std::string& filename="DoFs.txt"){
	std::ifstream readin(filename);
	for (int i = 0; i < ObjType1::cntCp; i++)
		for(int k=0;k<3;k++)
			readin>>CpPos1[i](k);
	for (int i = 0; i < ObjType1::cntCp; i++)
		for(int k=0;k<3;k++)
			readin>>CpVel1[i](k);
			
	for (int i = 0; i < ObjType2::cntCp; i++)
		for(int k=0;k<3;k++)
			readin>>CpPos2[i](k);
	for (int i = 0; i < ObjType2::cntCp; i++)
		for(int k=0;k<3;k++)
			readin>>CpVel2[i](k);
	readin.close();
}
template<typename ObjType1, typename ObjType2>
static void saveDoFs(const std::array<Vector3d, ObjType1::cntCp>& CpPos1, 
					const std::array<Vector3d, ObjType1::cntCp>& CpVel1,
					const std::array<Vector3d, ObjType2::cntCp>& CpPos2, 
					const std::array<Vector3d, ObjType2::cntCp>& CpVel2,
					const std::string& filename="DoFs.txt"){
	std::ofstream f(filename);
	f<< std::fixed << std::setprecision(10);
	for(auto item : CpPos1)
		f<<item.transpose()<<"\n";
	for(auto item : CpVel1)
		f<<item.transpose()<<"\n";
	for(auto item : CpPos2)
		f<<item.transpose()<<"\n";
	for(auto item : CpVel2)
		f<<item.transpose()<<"\n";
	f.close();
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