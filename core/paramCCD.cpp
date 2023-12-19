#include"paramCCD.h"
#include"triBezier.h"
// #include"recBezier.h"

template<const PatchType cntCp>
void ParamCCD<cntCp>::ParamCCD<cntCp>(const BoundingBoxType _bbtype):
	bbtype(_bbtype)
{
	switch(cntCp):{
		case PatchType::TriBezier:
			CpPos1 = new TriBezierObj(), 
			CpPos2 = new TriBezierObj(), 
			CpVel1 = new TriBezierObj(), 
			CpVel2 = new TriBezierObj();
			break;
		default:
			std::cerr<<"no such case!";
	}
	std::default_random_engine randGenerator(0);
	// std::default_random_engine randGenerator(std::random_device());
    Vector3d dir;
    for(int dim=0; dim<3; dim++) dir[dim] = randNormal(randGenerator);
    dir.normalize();
    dir/=1.625;
    // std::cout<<dir;

    for (int i = 0; i < cntCp; i++) {
        for(int dim=0; dim<3; dim++) CpPos1->ctrlp[i][dim] = randNormal(randGenerator);
        for(int dim=0; dim<3; dim++) CpVel1->ctrlp[i][dim] = randNormal(randGenerator);
        for(int dim=0; dim<3; dim++) CpPos2->ctrlp[i][dim] = randNormal(randGenerator);
        for(int dim=0; dim<3; dim++) CpVel2->ctrlp[i][dim] = randNormal(randGenerator);
        CpPos1->ctrlp[i]+=dir;
        CpVel1->ctrlp[i]-=dir;
        CpPos2->ctrlp[i]-=dir;
        CpVel2->ctrlp[i]+=dir;
	}
}

template<const PatchType cntCp>
double ParamCCD<cntCp>::PrimitiveCheck(ParamBound const &divUvB1, ParamBound const &divUvB2){
	auto const ptPos1 = CpPos1->divideBezierPatch(divUvB1);
	auto const ptVel1 = CpVel1->divideBezierPatch(divUvB1);
	auto const ptPos2 = CpPos2->divideBezierPatch(divUvB2);
	auto const ptVel2 = CpVel2->divideBezierPatch(divUvB2);

	auto setAxes = [&] (std::vector<Vector3d>& axes) {
		if(bbtype==BoundingBoxType::AABB){
			axes = {Vector3d::Unit(0), Vector3d::Unit(1), Vector3d::Unit(2)};
		}
		else if(bbtype==BoundingBoxType::OBB) {
			Vector3d 
			lu1 = (ptPos1[9]-ptPos1[0]),//v-w，即u的对边
			lv1 = (ptPos1[9]-ptPos1[3]+ptPos1[0]-ptPos1[3]),//中线
			ln1 = lu1.cross(lv1);
			Vector3d 
			lu2 = (ptPos2[9]-ptPos2[0]),
			lv2 = (ptPos2[9]-ptPos2[3]+ptPos2[0]-ptPos2[3]),
			ln2 = lu2.cross(lv2);
			axes = {lu1,lv1,ln1,lu2,lv2,ln2, 
				lu1.cross(lu2), lu1.cross(lv2), lu1.cross(ln2), 
				lv1.cross(lu2), lv1.cross(lv2), lv1.cross(ln2), 
				ln1.cross(lu2), ln1.cross(lv2), ln1.cross(ln2)};
		}
	};

	// std::cout<<"done!\n";
	std::vector<Vector3d> axes;
	axes.clear();
	setAxes(axes);
	std::vector<Array2d> feasibleIntvs;
	feasibleIntvs.clear();

	auto AxisCheck=[&](const std::array<Vector3d, CpPos1.cntCp>& p1, const std::array<Vector3d, CpVel1.cntCp>& v1, 
						const std::array<Vector3d, CpPos2.cntCp>& p2, const std::array<Vector3d, CpVel1.cntCp>& v2, 
						const Vector3d& axis){
		std::vector<Line> lines1, lines2;
		std::vector<Line> ch1, ch2;
		std::vector<double> pts1, pts2;
		lines1.clear(); lines2.clear();
		ch1.clear(); ch2.clear();
		pts1.clear(); pts2.clear();
		for(int i = 0; i < CpPos1.cntCp; i++) lines1.emplace_back(v1[i].dot(axis), p1[i].dot(axis));
		for(int i = 0; i < CpPos2.cntCp; i++) lines2.emplace_back(-v2[i].dot(axis), -p2[i].dot(axis));
		std::sort(lines1.begin(), lines1.end());
		std::sort(lines2.begin(), lines2.end());
		getCH(lines1, ch1, pts2);
		getCH(lines2, ch2, pts2);
		// std::cout<<"getCHOK!\n";
		for(auto & l:ch2)
			l.k = -l.k, l.b = -l.b;
		const auto intvT = linearCHIntersect(ch1, ch2, pts1, pts2);
		if(SHOWANS) std::cout<<intvT.transpose()<<"\n";
		if(DEBUG) std::cin.get();
		if(intvT[0]!=-1)feasibleIntvs.push_back(intvT);
	};

	for(const auto& axis:axes)
		AxisCheck(ptPos1, ptVel1, ptPos2, ptVel2, axis);
	for(const auto& axis:axes)
		AxisCheck(ptPos2, ptVel2, ptPos1, ptVel1, axis);

	if(SHOWANS) std::cout<<"done!\n";

	if (feasibleIntvs.size()==0) return 0; //这意味着整段时间都有碰撞

	//无碰撞发生的并，剩下的就是有碰撞发生的
	std::sort(feasibleIntvs.begin(), feasibleIntvs.end(), 
		[](const Array2d& intv1, const Array2d& intv2){
			return (intv1(0)<intv2(0));
		});
	if(feasibleIntvs.size()==0 || (feasibleIntvs[0](0)>0)) return 0;
	double minT = feasibleIntvs[0](1);
	for(int i=1;i<feasibleIntvs.size();i++)
		if(feasibleIntvs[i](0)<minT) //不能加等，因为无碰撞给的是开区间，如果有),(的情况加等号会把这个情况漏掉
			minT=std::max(minT, feasibleIntvs[i](1));
		else break;
	if(minT<DeltaT)return minT;
	else return -1;
}

template<const PatchType cntCp>
double ParamCCD<cntCp>::CCD() {
	using steady_clock = std::chrono::steady_clock;
	using duration = std::chrono::duration<double>;
	const auto initialTime = steady_clock::now();

	std::priority_queue<ParamPatchPair> heap;
	ParamBound initParam();
	initParam.initialize();
	double colTime = PrimitiveCheck(initParam, initParam);
	if (colTime>=0 && colTime<=DeltaT)
		heap.emplace(&initParam, &initParam, colTime);
	std::cout<<"done!\n";

	while (!heap.empty()) {
		auto const cur = heap.top();
		// std::cout << "patch1 : {" << cur.tpb1.nodes[0]<<"; "<< cur.tpb1.nodes[1]<<"; "<< cur.tpb1.nodes[2]<<"; " <<"}\n" 
		// 	<< " patch2 : {" << cur.tpb2.nodes[0]<<"; "<< cur.tpb2.nodes[1]<<"; "<< cur.tpb2.nodes[2]<<"; "<<"}\n";
		// std::cin.get();
		heap.pop();
		cnt++;
		if(DEBUG) std::cout<<cnt<<"\n";

		// Decide whether the algorithm converges
		// calcSquaredDist(cur) < MinSquaredDist || 
		if (std::max(cur.pb1.diameter(), cur.pb2.diameter()) < MinDeltaUV) {
		// if (calcSquaredDist(cur) < MinSquaredDist) {
			std::cout << cur.pb1.centerParam() << std::endl;
			std::cout << cur.pb2.centerParam() << std::endl;
			uv1 = cur.tpb1.centerParam();
			uv2 = cur.tpb2.centerParam();
			const auto endTime = steady_clock::now();
			std::cout << "min time: "<<  cur.tLower << "used seconds: " <<
				duration(endTime - initialTime).count()
				<< std::endl;
			return cur.tLower;
		}

		// Divide the current patch into four-to-four pieces
		for (int i = 0; i < 4; i++) {
			ParamBound divUvB1(cur.tpb1.interpSubpatchParam(i));
			for (int j = 0; j < 4; j++) {
				ParamBound divUvB2(cur.tpb2.interpSubpatchParam(j));
				colTime = PrimitiveCheck(divUvB1, divUvB2, bbtype);//maybe also need timeLB?
				if (colTime>=0 && colTime<=DeltaT){
					heap.emplace(divUvB1, divUvB2, colTime);
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