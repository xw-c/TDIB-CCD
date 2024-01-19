#include"config.h"
#include"recBezierMesh.h"
#include"triBezier.h"


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
		Vector3d lu1 = (posStart1[ParamObj1::cornerId(2)]-posStart1[ParamObj1::cornerId(0)]
				+posStart1[ParamObj1::cornerId(3)]-posStart1[ParamObj1::cornerId(1)]).normalized();//u延展的方向
		Vector3d lv1 = (posStart1[ParamObj1::cornerId(1)]-posStart1[ParamObj1::cornerId(0)]
				+posStart1[ParamObj1::cornerId(3)]-posStart1[ParamObj1::cornerId(2)]);//v延展的方向
		lv1 = (lv1-lv1.dot(lu1)*lu1).eval();
		Vector3d ln1 = lu1.cross(lv1);

		Vector3d lu2 = (posStart2[ParamObj2::cornerId(2)]-posStart2[ParamObj2::cornerId(0)]
				+posStart2[ParamObj2::cornerId(3)]-posStart2[ParamObj2::cornerId(1)]).normalized();//u延展的方向
		Vector3d lv2 = (posStart2[ParamObj2::cornerId(1)]-posStart2[ParamObj2::cornerId(0)]
				+posStart2[ParamObj2::cornerId(3)]-posStart2[ParamObj2::cornerId(2)]);//v延展的方向
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

	while (!heap.empty()) {
		auto const cur = heap.top();
		// std::cout << "patch1 : {" << cur.pb1.pMin.transpose()<<"; "<< cur.pb1.pMax.transpose()<<"; " <<"}\n" 
		// 	<< " patch2 : {" << cur.pb2.pMin.transpose()<<"; "<< cur.pb2.pMax.transpose()<<"; "<<"}\n";
		// std::cin.get();
		heap.pop();
		cnt++;
		if(DEBUG) std::cout<<cnt<<"\n";

		// Decide whether the algorithm converges
		// calcSquaredDist(cur) < MinSquaredDist || 
		if (cur.calcL1Dist(CpPos1, CpVel1, CpPos2, CpVel2) < MinL1Dist) {
		// if (calcSquaredDist(cur) < MinSquaredDist) {
			// std::cout << cur.pb1.centerParam() << std::endl;
			// std::cout << cur.pb2.centerParam() << std::endl;
			uv1 = cur.pb1.centerParam();
			uv2 = cur.pb2.centerParam();
			const auto endTime = steady_clock::now();
			std::cout << "min time: "<<  cur.tIntv[0] << "\nused seconds: " <<
				duration(endTime - initialTime).count()
				<< std::endl;
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
	std::cout << "used seconds: " <<
		duration(endTime - initialTime).count()
		<< std::endl;
	return -1;
}

static void singleTest(){
	RecCubicBezier obj1, obj2, vel1, vel2;

	const int Kase = 100;
	std::srand(0);
	using steady_clock = std::chrono::steady_clock;
	using duration = std::chrono::duration<double>;
	const auto initialTime = steady_clock::now();
	// kase 19
	for(int kase = 0;kase<Kase;kase++){
		for(auto &p:obj1.ctrlp)p=Vector3d::Random();
		for(auto &p:obj1.ctrlp)p[0]-=3;
		for(auto &p:obj2.ctrlp)p=Vector3d::Random();
		for(auto &p:obj2.ctrlp)p[0]+=3;
		for(auto &p:vel1.ctrlp)p=Vector3d(3,0,0);
		for(auto &p:vel2.ctrlp)p=Vector3d(-3,0,0);

		Array2d uv1,uv2;
		double t = inclusionCCD<RecCubicBezier,RecCubicBezier,RecParamBound,RecParamBound>
							(obj1,vel1,obj2,vel2,uv1,uv2, DeltaT);
		std::cout<<kase<<": "<<duration(steady_clock::now() - initialTime).count()<<"s\n";
	}
	const auto endTime = steady_clock::now();
	std::cout << "used seconds: " <<
		duration(endTime - initialTime).count()/Kase
		<< std::endl;
}
static void randomTest(){
	RecCubicBezier obj1, obj2, vel1, vel2;
	std::srand(0);
	std::default_random_engine randGenerator(0);
	const int Kase = 100;
	using steady_clock = std::chrono::steady_clock;
	using duration = std::chrono::duration<double>;
	const auto initialTime = steady_clock::now();
	for(int kase = 0;kase<Kase;kase++){
		Vector3d dir=Vector3d::Random().normalized()/1.625;
		for (int i = 0; i < RecCubicBezier::cntCp; i++) {
			for(int dim=0; dim<3; dim++) obj1.ctrlp[i][dim] = randNormal(randGenerator);
			for(int dim=0; dim<3; dim++) vel1.ctrlp[i][dim] = randNormal(randGenerator);
			for(int dim=0; dim<3; dim++) obj2.ctrlp[i][dim] = randNormal(randGenerator);
			for(int dim=0; dim<3; dim++) vel2.ctrlp[i][dim] = randNormal(randGenerator);
			obj1.ctrlp[i]+=dir;
			vel1.ctrlp[i]-=dir;
			obj2.ctrlp[i]-=dir;
			vel2.ctrlp[i]+=dir;
		}
		Array2d uv1,uv2;
		double t = inclusionCCD<RecCubicBezier,RecCubicBezier,RecParamBound,RecParamBound>
							(obj1,vel1,obj2,vel2,uv1,uv2, DeltaT);
		std::cout<<kase<<": "<<duration(steady_clock::now() - initialTime).count()<<"s\n";
	}
	const auto endTime = steady_clock::now();
	std::cout << "used seconds: " <<
		duration(endTime - initialTime).count()/Kase
		<< std::endl;
}
int main(){
	randomTest();
}