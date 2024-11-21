# pragma once
#include"config.h"

template<typename ParamObj1, typename ParamObj2>
void setAxes(const std::array<Vector3d, ParamObj1::cntCp>& ptPos1, 
				const std::array<Vector3d, ParamObj2::cntCp>& ptPos2,
				std::vector<Vector3d>& axes,
				const BoundingBoxType& bb){	
	if(bb==BoundingBoxType::AABB){
		axes = {Vector3d::Unit(0), Vector3d::Unit(1), Vector3d::Unit(2)};
	}
	else if(bb==BoundingBoxType::OBB){
		Vector3d lu1 = ParamObj1::axisU(ptPos1);
		Vector3d lv1tmp = ParamObj1::axisV(ptPos1);
		Vector3d ln1 = lu1.cross(lv1tmp);
		Vector3d lv1 = ln1.cross(lu1);

		Vector3d lu2 = ParamObj2::axisU(ptPos2);
		Vector3d lv2tmp = ParamObj2::axisV(ptPos2);
		Vector3d ln2 = lu2.cross(lv2tmp);
		Vector3d lv2 = ln2.cross(lu2);

		axes = {lu1,lv1,ln1,lu2,lv2,ln2, 
			lu1.cross(lu2), lu1.cross(lv2), lu1.cross(ln2), 
			lv1.cross(lu2), lv1.cross(lv2), lv1.cross(ln2), 
			ln1.cross(lu2), ln1.cross(lv2), ln1.cross(ln2)};
	}
}

template<typename ParamObj1, typename ParamObj2>
void setAxes(const std::array<Vector3d, ParamObj1::cntCp>& ptPos1, 
				const std::array<Vector3d, ParamObj1::cntCp>& ptVel1, 
				const std::array<Vector3d, ParamObj2::cntCp>& ptPos2,
				const std::array<Vector3d, ParamObj2::cntCp>& ptVel2,
				std::vector<Vector3d>& axes,
				const BoundingBoxType& bb,
				const double t = 0){	
	if(bb==BoundingBoxType::AABB){
		axes = {Vector3d::Unit(0), Vector3d::Unit(1), Vector3d::Unit(2)};
	}
	else if(bb==BoundingBoxType::OBB){
		Vector3d lu1 = ParamObj1::axisU(ptPos1) + t*ParamObj1::axisU(ptVel1);
		Vector3d lv1tmp = ParamObj1::axisV(ptPos1) + t*ParamObj1::axisV(ptVel1);
		Vector3d ln1 = lu1.cross(lv1tmp);
		Vector3d lv1 = ln1.cross(lu1);

		Vector3d lu2 = ParamObj2::axisU(ptPos2) + t*ParamObj2::axisU(ptVel2);
		Vector3d lv2tmp = ParamObj2::axisV(ptPos2) + t*ParamObj2::axisU(ptVel2);
		Vector3d ln2 = lu2.cross(lv2tmp);
		Vector3d lv2 = ln2.cross(lu2);

		axes = {lu1,lv1,ln1,lu2,lv2,ln2, 
			lu1.cross(lu2), lu1.cross(lv2), lu1.cross(ln2), 
			lv1.cross(lu2), lv1.cross(lv2), lv1.cross(ln2), 
			ln1.cross(lu2), ln1.cross(lv2), ln1.cross(ln2)};
	}
}

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