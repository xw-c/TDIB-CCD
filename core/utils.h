#include "config.h"

static void readinDoFs(std::array<Vector3d, 16>& CpPos1, std::array<Vector3d, 16>& CpVel1, std::array<Vector3d, 16>& CpPos2, std::array<Vector3d, 16>& CpVel2){
	std::ifstream readin("DoFs.txt");
	for (int i = 0; i < 16; i++)
		for(int k=0;k<3;k++)
			readin>>CpPos1[i](k);
	for (int i = 0; i < 16; i++)
		for(int k=0;k<3;k++)
			readin>>CpPos2[i](k);
	for (int i = 0; i < 16; i++)
		for(int k=0;k<3;k++)
			readin>>CpVel1[i](k);
	for (int i = 0; i < 16; i++)
		for(int k=0;k<3;k++)
			readin>>CpVel2[i](k);
	readin.close();
}
static void saveDoFs(const std::array<Vector3d, 16>& CpPos1, const std::array<Vector3d, 16>& CpVel1, const std::array<Vector3d, 16>& CpPos2, const std::array<Vector3d, 16>& CpVel2){
	std::ofstream f("DoFs.txt");
	for(auto item : CpPos1)
		f<<item.transpose()<<"\n";
	for(auto item : CpPos2)
		f<<item.transpose()<<"\n";
	for(auto item : CpVel1)
		f<<item.transpose()<<"\n";
	for(auto item : CpVel2)
		f<<item.transpose()<<"\n";
	f.close();
}