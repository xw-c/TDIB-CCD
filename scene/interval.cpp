#include"mathOps.h"
#include"recBezierMesh.h"
#include"recRatBezierMesh.h"
#include"triBezier.h"
#include"utils.h"


template<typename ParamObj1, typename ParamObj2, typename ParamBound1, typename ParamBound2>
static bool inclusionCheck(const ParamObj1 &CpPos1, const ParamObj1 &CpVel1, 
						const ParamObj2 &CpPos2, const ParamObj2 &CpVel2,
						const ParamBound1 &divUvB1, const ParamBound2 &divUvB2,
						const Array2d divTime = Array2d(0,DeltaT)){
	auto posStart1 = CpPos1.divideBezierPatch(divUvB1), posEnd1 = posStart1;
	auto ptVel1 = CpVel1.divideBezierPatch(divUvB1);
	auto posStart2 = CpPos2.divideBezierPatch(divUvB2), posEnd2 = posStart2;
	auto ptVel2 = CpVel2.divideBezierPatch(divUvB2);
	for(int i=0;i<ParamObj1::cntCp;i++){
		posStart1[i]+=ptVel1[i]*divTime[0],
		posEnd1[i]+=ptVel1[i]*divTime[1];
	}
	for(int i=0;i<ParamObj2::cntCp;i++){
		posStart2[i]+=ptVel2[i]*divTime[0],
		posEnd2[i]+=ptVel2[i]*divTime[1];
	}

	// std::cout<<"done!\n";
	std::vector<Vector3d> axes;
	// axes.clear();
	if(bbtype==BoundingBoxType::AABB){
		axes = {Vector3d::Unit(0), Vector3d::Unit(1), Vector3d::Unit(2)};
	}
	else if(bbtype==BoundingBoxType::DOP14){
		axes = {Vector3d::Unit(0), Vector3d::Unit(1), Vector3d::Unit(2), 
				Vector3d(1,1,1).normalized(), Vector3d(-1,1,1).normalized(), Vector3d(-1,-1,1).normalized()};
	}
	else if(bbtype==BoundingBoxType::OBB){
		Vector3d lu1 = ParamObj1::axisU(posStart1).normalized();//u延展的方向
		Vector3d lv1 = ParamObj1::axisV(posStart1);//v延展的方向
		lv1 = (lv1-lv1.dot(lu1)*lu1).eval();
		Vector3d ln1 = lu1.cross(lv1);

		Vector3d lu2 = ParamObj2::axisU(posStart2).normalized();//u延展的方向
		Vector3d lv2 = ParamObj2::axisV(posStart2);//v延展的方向
		lv2 = (lv2-lv2.dot(lu2)*lu2).eval();
		Vector3d ln2 = lu2.cross(lv2);

		axes = {lu1,lv1,ln1,lu2,lv2,ln2, 
			lu1.cross(lu2), lu1.cross(lv2), lu1.cross(ln2), 
			lv1.cross(lu2), lv1.cross(lv2), lv1.cross(ln2), 
			ln1.cross(lu2), ln1.cross(lv2), ln1.cross(ln2)};
	}

	for(auto& axis:axes){
		double maxProj1 = -std::numeric_limits<double>::infinity(), minProj1 = std::numeric_limits<double>::infinity();
		for(const auto&p:posStart1){
			maxProj1 = std::max(maxProj1, p.dot(axis));
			minProj1 = std::min(minProj1, p.dot(axis));
		}
		for(const auto&p:posEnd1){
			maxProj1 = std::max(maxProj1, p.dot(axis));
			minProj1 = std::min(minProj1, p.dot(axis));
		}
		double maxProj2 = -std::numeric_limits<double>::infinity(), minProj2 = std::numeric_limits<double>::infinity();
		for(const auto&p:posStart2){
			maxProj2 = std::max(maxProj2, p.dot(axis));
			minProj2 = std::min(minProj2, p.dot(axis));
		}
		for(const auto&p:posEnd2){
			maxProj2 = std::max(maxProj2, p.dot(axis));
			minProj2 = std::min(minProj2, p.dot(axis));
		}
		if(maxProj2<minProj1 || maxProj1<minProj2) return 0;
	}
	return 1;
}

template<typename ParamObj1, typename ParamObj2, typename ParamBound1, typename ParamBound2>
static double inclusionCCD(const ParamObj1 &CpPos1, const ParamObj1 &CpVel1, 
						const ParamObj2 &CpPos2, const ParamObj2 &CpVel2,
						Array2d& uv1, Array2d& uv2, 
						const double upperTime = DeltaT) {
	struct PatchPair{
		ParamBound1 pb1;
		ParamBound2 pb2;
		Array2d tIntv;
		PatchPair(const ParamBound1& c1, const ParamBound2& c2, 
				Array2d t = Array2d(0,DeltaT)): pb1(c1), pb2(c2), tIntv(t) {}
		bool operator<(PatchPair const &o) const { return tIntv[0] > o.tIntv[0]; }

		double calcDist(const ParamObj1 &CpPos1, const ParamObj1 &CpVel1, 
						const ParamObj2 &CpPos2, const ParamObj2 &CpVel2) const{
			Vector3d const p1 = CpPos1.evaluatePatchPoint(pb1.centerParam());
			Vector3d const v1 = CpVel1.evaluatePatchPoint(pb1.centerParam());
			Vector3d const p2 = CpPos2.evaluatePatchPoint(pb2.centerParam());
			Vector3d const v2 = CpVel2.evaluatePatchPoint(pb2.centerParam());
			double t = (tIntv[0] + tIntv[1]) * 0.5;
        	Vector3d const pt1=(v1*t+p1), pt2=(v2*t+p2);
			return (pt2-pt1).norm();
		}
		double calcL1Dist(const ParamObj1 &CpPos1, const ParamObj1 &CpVel1, 
						const ParamObj2 &CpPos2, const ParamObj2 &CpVel2) const{
			auto const ptPos1 = CpPos1.divideBezierPatch(pb1);
			auto const ptPos2 = CpPos2.divideBezierPatch(pb2);
			double d=0;
			for(int axis=0;axis<3;axis++){
				double maxv = ptPos1[0][axis], minv=maxv;
				for(int i = 1; i < ParamObj1::cntCp; i++) {
					maxv=std::max(maxv, ptPos1[i][axis]);
					minv=std::min(minv, ptPos1[i][axis]);
				}
				d=std::max(d,maxv-minv);
			}
			for(int axis=0;axis<3;axis++){
				double maxv = ptPos2[0][axis], minv=maxv;
				for(int i = 1; i < ParamObj2::cntCp; i++) {
					maxv=std::max(maxv, ptPos2[i][axis]);
					minv=std::min(minv, ptPos2[i][axis]);
				}
				d=std::max(d,maxv-minv);
			}
			return d;
		}
	};

	using steady_clock = std::chrono::steady_clock;
	using duration = std::chrono::duration<double>;
	const auto initialTime = steady_clock::now();

	std::priority_queue<PatchPair> heap;
	ParamBound1 initParam1;
	ParamBound2 initParam2;
	double colTime = inclusionCheck(CpPos1, CpVel1, CpPos2, CpVel2, initParam1, initParam2);
	if (inclusionCheck(CpPos1, CpVel1, CpPos2, CpVel2, initParam1, initParam2))
		heap.emplace(initParam1, initParam2);
	// std::cout<<"done!\n";
	cnt=1;
	while (!heap.empty()) {
		auto const cur = heap.top();
		// std::cout << "patch1 : {" << cur.pb1.pMin.transpose()<<"; "<< cur.pb1.pMax.transpose()<<"; " <<"}\n" 
		// 	<< " patch2 : {" << cur.pb2.pMin.transpose()<<"; "<< cur.pb2.pMax.transpose()<<"; "<<"}\n";
		// std::cin.get();
		heap.pop();
		cnt++;
		if(DEBUG) std::cout<<cnt<<"\n";
			// std::cout << cur.pb1.centerParam().transpose() << std::endl;
			// std::cout << cur.pb2.centerParam().transpose() << std::endl;
			// std::cout << cur.tIntv.transpose() << std::endl;
		// Decide whether the algorithm converges
		// calcSquaredDist(cur) < MinSquaredDist || 
		if (cur.calcL1Dist(CpPos1, CpVel1, CpPos2, CpVel2) < MinL1Dist) {
		// if (calcSquaredDist(cur) < MinSquaredDist) {
			// std::cout << cur.pb1.centerParam() << std::endl;
			// std::cout << cur.pb2.centerParam() << std::endl;
			uv1 = cur.pb1.centerParam();
			uv2 = cur.pb2.centerParam();
			const auto endTime = steady_clock::now();
			// std::cout << "min time: "<<  cur.tIntv[0] 
			// 	<< "\nused seconds: " << duration(endTime - initialTime).count()
			// 	<< std::endl;
			return cur.tIntv[0];
		}

		// Divide the current patch into four-to-four pieces
		double tMid = (cur.tIntv[0]+cur.tIntv[1])*0.5;
		Array2d divTime1(cur.tIntv[0],tMid), divTime2(tMid, cur.tIntv[1]);
		for (int i = 0; i < 4; i++) {
			ParamBound1 divUvB1(cur.pb1.interpSubpatchParam(i));
			for (int j = 0; j < 4; j++) {
				ParamBound2 divUvB2(cur.pb2.interpSubpatchParam(j));
				if (inclusionCheck(CpPos1, CpVel1, CpPos2, CpVel2, divUvB1, divUvB2, divTime1)){
					heap.emplace(divUvB1, divUvB2, divTime1);
				}
				if (inclusionCheck(CpPos1, CpVel1, CpPos2, CpVel2, divUvB1, divUvB2, divTime2)){
					heap.emplace(divUvB1, divUvB2, divTime2);
				}
			}
		}
	}

	const auto endTime = steady_clock::now();
	// std::cout << "used seconds: " <<
	// 	duration(endTime - initialTime).count()
	// 	<< std::endl;
	return -1;
}

template<typename ObjType, typename ParamType>
static void singleTest(){
	ObjType obj1, obj2, vel1, vel2;

	using steady_clock = std::chrono::steady_clock;
	using duration = std::chrono::duration<double>;
	// kase 19
	Array2d uv1,uv2;
	readinDoFs<ObjType>(obj1.ctrlp, vel1.ctrlp, obj2.ctrlp, vel2.ctrlp);
	const auto initialTime = steady_clock::now();
	double t = inclusionCCD<ObjType,ObjType,ParamType,ParamType>
						(obj1,vel1,obj2,vel2,uv1,uv2, DeltaT);
	const auto endTime = steady_clock::now();
	std::cout << "used seconds: " <<
		duration(endTime - initialTime).count()
		<< std::endl;
	Vector3d const p1 = obj1.evaluatePatchPoint(uv1);
	Vector3d const v1 = vel1.evaluatePatchPoint(uv1);
	Vector3d const p2 = obj2.evaluatePatchPoint(uv2);
	Vector3d const v2 = vel2.evaluatePatchPoint(uv2);
	Vector3d const pt1=(v1*t+p1), pt2=(v2*t+p2);
	std::cout<<"delta: "<<(pt2-pt1).norm()<<"\n";

	// RecBezierMesh obj(2);
	// for(int i=0;i<16;i++)obj1.ctrlp[i]+=t*vel1.ctrlp[i];
	// for(int i=0;i<16;i++)obj2.ctrlp[i]+=t*vel2.ctrlp[i];
	// obj.patches[0]=obj1;
	// obj.patches[1]=obj2;
	// obj.writeObj("check-intv.obj");
}
template<typename ObjType, typename ParamType>
static void randomTest(const double denom){
	ObjType obj1, obj2, vel1, vel2;
	std::srand(0);
	using steady_clock = std::chrono::steady_clock;
	using duration = std::chrono::duration<double>;
	const auto initialTime = steady_clock::now();
	int hasCol = 0;
	for(int kase = 0;kase<Kase;kase++){
		generatePatchPair<ObjType>(obj1.ctrlp, vel1.ctrlp, obj2.ctrlp, vel2.ctrlp, denom);
		Array2d uv1,uv2;
		double t = inclusionCCD<ObjType,ObjType,ParamType,ParamType>
							(obj1,vel1,obj2,vel2,uv1,uv2, DeltaT);
		// std::cout<<cnt<<"\n";
		if(t>=0)hasCol++;
		// std::cout<<kase<<": "<<duration(steady_clock::now() - initialTime).count()<<"s\n";
	}
	const auto endTime = steady_clock::now();
	std::cout << hasCol<<" pairs have collided.\n";
	std::cout << "used seconds: " <<
		duration(endTime - initialTime).count()/Kase
		<< std::endl;
}
int main(){
	randomTest<TriLinearBezier, TriParamBound>(3);
	randomTest<TriQuadBezier, TriParamBound>(2);
	randomTest<TriCubicBezier, TriParamBound>(1.625);
	randomTest<RecLinearBezier, RecParamBound>(3);
	randomTest<RecQuadBezier, RecParamBound>(2);
	randomTest<RecCubicBezier, RecParamBound>(1.625);
	randomTest<RecQuadRatBezier, RecParamBound>(2);
	randomTest<RecCubicRatBezier, RecParamBound>(1.625);
}