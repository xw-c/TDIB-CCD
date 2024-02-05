#include"linearTriMeshCCD.h"
#include"config.h"
void testSegment(){
	double u1,u2;

    using steady_clock = std::chrono::steady_clock;
	using duration = std::chrono::duration<double>;
	const auto initialTime = steady_clock::now();

    for(int i=0;i<Kase;i++){
        Eigen::Vector3d dir=Eigen::Vector3d::Random().normalized()/3.,e10,e11,e20,e21,v10,v11,v20,v21;
        for(int dim=0; dim<3; dim++) {
            e10[dim] = randNormal(randGenerator);
            e11[dim] = randNormal(randGenerator);
            e20[dim] = randNormal(randGenerator);
            e21[dim] = randNormal(randGenerator);
            v10[dim] = randNormal(randGenerator);
            v11[dim] = randNormal(randGenerator);
            v20[dim] = randNormal(randGenerator);
            v21[dim] = randNormal(randGenerator);
        }
        const double t = modifiedEETest(e10+dir,e11+dir,v10-dir,v11-dir,e20-dir,e21-dir,v20+dir,v21+dir,u1,u2);
		std::cout<<t<<"\n";
    }
    const auto endTime = steady_clock::now();
	// std::cout<<"colTime: "<<t_max<<"\n";
	std::cout << "used ave seconds: " <<
		(duration(endTime - initialTime).count())/Kase
		<< std::endl;
}
int main(){
	// testBunny_EE_VF();
	testSegment();
}