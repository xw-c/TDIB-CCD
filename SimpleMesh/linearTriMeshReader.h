#include"linearTriMeshCCD.h"
class SimpleMesh{
public: 
    std::vector<Eigen::Vector3d> verts;
    std::vector<Eigen::Vector3i> faces;
    std::vector<Eigen::Vector2i> edges;
    int cntV, cntF, cntE;
    SimpleMesh(const std::string& filename){
        std::cout<<filename;
        std::ifstream in(filename);
        if(!in.is_open()){std::cout<<"open file failed!\n";exit(-1);}
        in >> cntV >> cntF;

        verts.clear();
        faces.clear();
        edges.clear();
        char c;
        for(int i=0;i<cntV;i++){
            Eigen::Vector3d vert;
            in >> c >> vert[0] >> vert[1] >> vert[2];
            verts.push_back(vert);
        }
        for(int i=0;i<cntF;i++){
            Eigen::Vector3i face;
            in >> c >> face[0] >> face[1] >> face[2];
            face[0]--,face[1]--,face[2]--;
            faces.emplace_back(face);
            for(int i=0;i<3;i++){
                // int a=std::min(face[i],face[(i+1)%3]),b=std::max(face[i],face[(i+1)%3]);
                int a=face[i],b=face[(i+1)%3];
                if(a<b)edges.emplace_back(a,b);
            }
        }
	    // edges.erase(std::unique(edges.begin(), edges.end()), edges.end()); // 去重
        cntE=edges.size();

        in.close();
    }
    void writeObj(const std::string& filename) const {
		std::ofstream out(filename);
		for(auto&vert:verts)
			out<<"v "<<vert[0]<<" "<<vert[1]<<" "<<vert[2]<<"\n";
		for(auto&face:faces)
			out<<"f "<<face[0]+1<<" "<<face[1]+1<<" "<<face[2]+1<<"\n";
		out.close();
	}
	void moveObj(const Eigen::Vector3d& dis){
		for(auto&vert:verts)
            vert+=dis;
	}
	void rotateObj(const double& angle, const Eigen::Vector3d& axis, const Eigen::Vector3d& origin = Eigen::Vector3d::Zero()){
		//似乎带仿射
		Eigen::AngleAxisd rot(angle, axis);
		for(auto&vert:verts)
            vert = (rot.matrix()*(vert-origin)+origin).eval();
	}
	void setObj(const Vector3d& vel){
		for(auto& v:verts)v=vel;
	}
};


// void testBunny_EE_VF(){
// 	SimpleMesh obj1("bunny292.obj");
// 	SimpleMesh obj2 = obj1, vel1=obj1, vel2=obj1;
// 	Vector3d displace(6,6,0);
// 	obj2.moveObj(displace*0.5);
// 	obj1.rotateObj(PI, Vector3d(0,1,0)+Vector3d::Random());
// 	vel1.setObj(displace);
// 	vel2.setObj(Vector3d::Zero());
// 	std::cout<<"preprocess done!\n";
// 	using steady_clock = std::chrono::steady_clock;
// 	using duration = std::chrono::duration<double>;
// 	const auto initialTime = steady_clock::now();
	
// 	double minT = DeltaT;
// 	{
//         // std::cin.get();
//         for(int e0=0;e0<obj1.cntV;e0++)
//             for(int e1=0;e1<obj2.cntF;e1++){
// 				Point CpPos1 = obj1.verts[e0];
// 				Point CpVel1 = vel1.verts[e0];
// 				Face CpPos2({obj2.verts[obj2.faces[e1][0]], obj2.verts[obj2.faces[e1][1]], obj2.verts[obj2.faces[e1][2]]});
// 				Face CpVel2({vel2.verts[vel2.faces[e1][0]], vel2.verts[vel2.faces[e1][1]], vel2.verts[vel2.faces[e1][2]]});
//                 Array2d pb;
// 				const double toi = VFTest(CpPos1, CpVel1, CpPos2, CpVel2, pb, minT);
//                 if(toi>0)minT=std::min(toi, minT);
//                 // std::cout<<res<<": "<<t_max<<"\n";
//             }
//         std::cout<<"VF done!\n";
// 		std::cout<<"colTime: "<<minT<<"\n";
//         for(int e0=0;e0<obj2.cntV;e0++)
//             for(int e1=0;e1<obj1.cntF;e1++){
// 				Point CpPos1 = obj2.verts[e0];
// 				Point CpVel1 = vel2.verts[e0];
// 				Face CpPos2({obj1.verts[obj1.faces[e1][0]], obj1.verts[obj1.faces[e1][1]], obj1.verts[obj1.faces[e1][2]]});
// 				Face CpVel2({vel1.verts[vel1.faces[e1][0]], vel1.verts[vel1.faces[e1][1]], vel1.verts[vel1.faces[e1][2]]});
//                 Array2d pb;
// 				const double toi = VFTest(CpPos1, CpVel1, CpPos2, CpVel2, pb, minT);
//                 if(toi>0)minT=std::min(toi, minT);
//                 // std::cout<<res<<": "<<t_max<<"\n";
//             }
//         std::cout<<"FV done!\n";
// 		std::cout<<"colTime: "<<minT<<"\n";
//         // std::cin.get();
//         for(int e0=0;e0<obj1.cntE;e0++)
//             for(int e1=0;e1<obj2.cntE;e1++){
// 				Edge CpPos1(obj1.verts[obj1.edges[e0][0]], obj1.verts[obj1.edges[e0][1]]);
// 				Edge CpVel1(vel1.verts[obj1.edges[e0][0]], vel1.verts[obj1.edges[e0][1]]);
// 				Edge CpPos2(obj2.verts[obj2.edges[e1][0]], obj2.verts[obj1.edges[e1][1]]);
// 				Edge CpVel2(vel2.verts[obj1.edges[e1][0]], vel2.verts[obj1.edges[e1][1]]);
//                 double u,v;
// 				const double toi = EETest(CpPos1, CpVel1, CpPos2, CpVel2, u,v, minT);
//                 if(toi>0)minT=std::min(toi, minT);
//             }
//         std::cout<<"EE done!\n";
// 		std::cout<<"colTime: "<<minT<<"\n";
//     }


// 	const auto endTime = steady_clock::now();
// 	std::cout<<"colTime: "<<minT<<"\n";
// 	std::cout << "used seconds: " <<
// 		duration(endTime - initialTime).count()
// 		<< std::endl;

// 	obj1.moveObj(displace*minT);
// 	// obj1.writeObj("bunnyCol1.obj");
// 	// obj2.writeObj("bunnyCol2.obj");
// }
