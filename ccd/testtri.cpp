#include "triBezier.h"
#include "config.h"

#include <algorithm>
#include <array>
#include <string>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <limits>
#include <queue>
#include <span>
#include <chrono>

static constexpr int cntCp = 10;

TriBezier CpPos1, CpPos2, CpVel1, CpVel2;
BaryCoord uv1,uv2;

struct Line
{
	double k, b;
	Line(const double& k,const double& b): k(k), b(b) {}
	bool operator<(const Line &l) const { 
		return k < l.k || (k == l.k && b > l.b); // 相同斜率的直线中只有截距最大的被留下来
	}
	bool operator==(const Line &l) const {return k == l.k;}

	double lineIntersect_x(const Line &l) const {
		if(k==l.k){
			std::cout<<"parallel lines do not intersect at a single point!\n";
			exit(-1);
		}
		return -(b-l.b)/(k-l.k);
	}

};

static void getCH(std::vector<Line>& lines, std::vector<Line>& ch, std::vector<double>& pts) {
	lines.erase(std::unique(lines.begin(), lines.end()), lines.end()); // 去重
	ch.clear();
	pts.clear();
	pts.push_back(0);
	ch.push_back(lines[0]);
	if(DEBUG) for(const auto&l:lines)std::cout<<"lines:\t"<<l.k<<" "<<l.b<<"\n";
	double id = 1, intsctX = 0;
	while(id < lines.size()){
		// std::cout<<id<<"  "<<pts.size()<<"\n";
		while(!ch.empty()){
			intsctX = lines[id].lineIntersect_x(ch.back());
			if(intsctX<=pts.back()){
			// if(ch.back().k*pts.back()+ch.back().b<=lines[id].k*pts.back()+lines[id].b){
				pts.pop_back();
				ch.pop_back();
			}
			else break;
		}
		ch.push_back(lines[id]);
		pts.push_back(std::max(0.,intsctX));
		id++;
	}
	while(pts.back()>=DeltaT){
		pts.pop_back();
		ch.pop_back();
	}
	pts.push_back(DeltaT);
	if(DEBUG) for(const auto&l:ch)std::cout<<"ch:\t"<<l.k<<" "<<l.b<<"\n";
	if(DEBUG) for(const auto&pt:pts)std::cout<<"pt:\t"<<pt<<"\n";
	if(ch.empty()){
		std::cout<<"empty CH!\n";
		exit(-1);
	}
	if(ch.size()+1!=pts.size()){
		std::cout<<"segments and inflections are not compatible!\n";
		exit(-1);
	}
}

static Array2d linearCHIntersect(const std::vector<Line>& ch1, const std::vector<Line>& ch2, 
								const std::vector<double>& pts1, const std::vector<double>& pts2) {
	int id1=1, id2=1;
	double intvL=-1, intvR=-1;
	double sweep=0, lastsweep=0;
	bool stopInAdv = false; // 还没写上，但是比如：intvL==-1&&ch1[id1-1].k>ch1[id1-1].k

	auto checkSweepLine = [&] (const int id1, const int id2) {
		double y1=ch1[id1-1].k*sweep+ch1[id1-1].b;
		double y2=ch2[id2-1].k*sweep+ch2[id2-1].b;
		if (y1>y2){
			if(intvL!=-1)
				intvR = ch1[id1-1].lineIntersect_x(ch2[id2-1]);
		} // 如果k1==k2，那么必然前一个节点必然已经满足y1<y2了
		else if(y1<y2){
			if(intvL==-1)
				intvL = ch1[id1-1].lineIntersect_x(ch2[id2-1]); 
		}// 如果k1==k2，那么必然前一个节点必然已经满足y1<y2了
		if (DEBUG) std::cout<<id1<<"  "<<id2<<" /  "<<y1<<"  "<<y2<<" /  "<<intvL<<" "<<intvR<<"\n";
		lastsweep = sweep;
	};

	if(ch1[0].b<ch2[0].b)
		intvL = 0;
	// else if(ch1[0].b==ch2[0].b)
	// 	intvL = intvR = 0;
	while(id1<pts1.size()&&id2<pts1.size()){
		sweep = std::min(pts1[id1], pts2[id2]);
		if(sweep!=lastsweep)//并不知道这样跳过能不能更快
			checkSweepLine(id1, id2);
		if(pts1[id1] < pts2[id2]) id1++;
		else id2++;
	}
	while(id1<pts1.size()){
		sweep = pts1[id1];
		if(sweep!=lastsweep)//并不知道这样跳过能不能更快
			checkSweepLine(id1, id2);
		id1++;
	}
	while(id2<pts2.size()){
		sweep = pts2[id2];
		if(sweep!=lastsweep)//并不知道这样跳过能不能更快axesOBB
			checkSweepLine(id1, id2);
		id2++;
	}
	if(intvL!=-1 && intvR==-1)intvR=DeltaT;
	return Array2d(intvL,intvR);
}

// maybe no need to contruct from scratch?
static double PrimitiveCheck(TriParamBound const &divUvB1, TriParamBound const &divUvB2, const BB bbtype){
	auto const ptPos1 = CpPos1.divideBezierPatch(divUvB1);
	auto const ptVel1 = CpVel1.divideBezierPatch(divUvB1);
	auto const ptPos2 = CpPos2.divideBezierPatch(divUvB2);
	auto const ptVel2 = CpVel2.divideBezierPatch(divUvB2);

	auto setAxes = [&] (std::vector<Vector3d>& axes) {
		if(bbtype==BB::AABB){
			axes = {Vector3d::Unit(0), Vector3d::Unit(1), Vector3d::Unit(2)};
		}
		else if(bbtype==BB::OBB) {
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

	auto AxisCheck=[&](const std::array<Vector3d, cntCp>& p1, const std::array<Vector3d, cntCp>& v1, 
						const std::array<Vector3d, cntCp>& p2, const std::array<Vector3d, cntCp>& v2, 
						const Vector3d& axis){
		std::vector<Line> lines1, lines2;
        std::vector<Line> ch1, ch2;
        std::vector<double> pts1, pts2;
		lines1.clear(); lines2.clear();
		ch1.clear(); ch2.clear();
		pts1.clear(); pts2.clear();
		for(int i = 0; i < cntCp; i++) lines1.emplace_back(v1[i].dot(axis), p1[i].dot(axis));
		for(int i = 0; i < cntCp; i++) lines2.emplace_back(-v2[i].dot(axis), -p2[i].dot(axis));
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

double calcSquaredDist(const TriPatchPair& tpp){
    Vector3d const p1 = CpPos1.blossomBicubicBezier(tpp.tpb1.centerParam());
    Vector3d const v1 = CpVel1.blossomBicubicBezier(tpp.tpb1.centerParam());
    Vector3d const p2 = CpPos2.blossomBicubicBezier(tpp.tpb2.centerParam());
    Vector3d const v2 = CpVel2.blossomBicubicBezier(tpp.tpb2.centerParam());
    Vector3d const pt1=(v1*tpp.tLower+p1), pt2=(v2*tpp.tLower+p2);
    return (p1-p2).squaredNorm();
}

static double ccd(const BB bbtype) {
	using steady_clock = std::chrono::steady_clock;
	using duration = std::chrono::duration<double>;
	const auto initialTime = steady_clock::now();

	std::priority_queue<TriPatchPair> heap;
	TriParamBound initParam(BaryCoord(1,0,0), BaryCoord(0,1,0), BaryCoord(0,0,1));
    double colTime = PrimitiveCheck(initParam, initParam, bbtype);
	if (colTime>=0 && colTime<=DeltaT)
		heap.emplace(initParam, initParam, colTime);
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
		if (calcSquaredDist(cur) < MinSquaredDist || std::max(cur.tpb1.maxEdgeDist(), cur.tpb2.maxEdgeDist()) < MinDeltaUV) {
        // if (calcSquaredDist(cur) < MinSquaredDist) {
			std::cout << cur.tpb1.centerParam() << std::endl;
			std::cout << cur.tpb2.centerParam() << std::endl;
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
			TriParamBound divUvB1(cur.tpb1.interpSubpatchParam(i));
			for (int j = 0; j < 4; j++) {
			    TriParamBound divUvB2(cur.tpb2.interpSubpatchParam(j));
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

void saveDoFs(){
	std::ofstream f("DoFs.txt");
	for(auto item : CpPos1.ctrlp)
		f<<item.transpose()<<"\n";
	for(auto item : CpPos2.ctrlp)
		f<<item.transpose()<<"\n";
	for(auto item : CpVel1.ctrlp)
		f<<item.transpose()<<"\n";
	for(auto item : CpVel2.ctrlp)
		f<<item.transpose()<<"\n";
	f.close();
}

void generateTriBeizer(){
    // std::srand(0);
    std::random_device rd;
    std::default_random_engine randGenerator(rd());
    Vector3d dir;
    for(int dim=0; dim<3; dim++) dir[dim] = randNormal(randGenerator);
    dir.normalize();
    dir/=1.625;
    // std::cout<<dir;

    for (int i = 0; i < cntCp; i++) {
        for(int dim=0; dim<3; dim++) CpPos1.ctrlp[i][dim] = randNormal(randGenerator);
        for(int dim=0; dim<3; dim++) CpVel1.ctrlp[i][dim] = randNormal(randGenerator);
        for(int dim=0; dim<3; dim++) CpPos2.ctrlp[i][dim] = randNormal(randGenerator);
        for(int dim=0; dim<3; dim++) CpVel2.ctrlp[i][dim] = randNormal(randGenerator);
        CpPos1.ctrlp[i]+=dir;
        CpVel1.ctrlp[i]-=dir;
        CpPos2.ctrlp[i]-=dir;
        CpVel2.ctrlp[i]+=dir;
	}

    // saveDoFs();
    // for (int i = 0; i < 10; i++) {
	// 	CpPos1.ctrlp[i] = Vector3d::Random() - Vector3d::Constant(.6);
	// 	CpVel1.ctrlp[i] = Vector3d::Random()*0.3 + Vector3d::Constant(.6);
	// 	CpPos2.ctrlp[i] = Vector3d::Random() + Vector3d::Constant(.6);
	// 	CpVel2.ctrlp[i] = Vector3d::Random()*0.3 - Vector3d::Constant(.6);
	// }
}
void randomTest(){
	using steady_clock = std::chrono::steady_clock;
	using duration = std::chrono::duration<double>;
	int cntAABB=0, cntSample=0;
	const int Kase = 100;
	double ans[2][Kase];

	const auto initOBB = steady_clock::now();
	std::srand(0);
	for(int kase=0;kase<Kase;kase++){
		generateTriBeizer();
		ans[0][kase]=ccd(BB::OBB);
	}
	const auto endOBB = steady_clock::now();
	std::cout<<"OBB used seconds: "<<duration(endOBB - initOBB).count()/Kase<<"\n";

	const auto initAABB = steady_clock::now();
	// std::srand(0);
	for(int kase=0;kase<Kase;kase++){
		generateTriBeizer();
		ans[0][kase]=ccd(BB::AABB);
	}
	const auto endAABB = steady_clock::now();
	std::cout<<"AABB used seconds: "<<duration(endAABB - initAABB).count()/Kase<<"\n";
	std::cout<<"OBB used seconds: "<<duration(endOBB - initOBB).count()/Kase<<"\n";
}
void certificate(){
	std::srand(std::time(nullptr));
    generateTriBeizer();

    {
        double t = ccd(BB::OBB);

        Vector3d const p1 = CpPos1.blossomBicubicBezier(uv1);
        Vector3d const v1 = CpVel1.blossomBicubicBezier(uv1);
        Vector3d const p2 = CpPos2.blossomBicubicBezier(uv2);
        Vector3d const v2 = CpVel2.blossomBicubicBezier(uv2);
        Vector3d const pt1=(v1*t+p1), pt2=(v2*t+p2);
        std::cout<<"pos at: "<<pt1.transpose()<<"     "<<pt2.transpose()<<"\n";
        std::cout<<"delta: "<<(pt2-pt1).norm()<<"\n";
    }

    // {
    //     double t = ccd(BB::AABB);

    //     Vector3d const p1 = CpPos1.blossomBicubicBezier(uv1);
    //     Vector3d const v1 = CpVel1.blossomBicubicBezier(uv1);
    //     Vector3d const p2 = CpPos2.blossomBicubicBezier(uv2);
    //     Vector3d const v2 = CpVel2.blossomBicubicBezier(uv2);
    //     Vector3d const pt1=(v1*t+p1), pt2=(v2*t+p2);
    //     std::cout<<"pos at: "<<pt1.transpose()<<"     "<<pt2.transpose()<<"\n";
    //     std::cout<<"delta: "<<(pt2-pt1).norm()<<"\n";
    // }
}
int main(){
    randomTest();

	// generateTriBeizer();
    // std::cout<<ccd(BB::AABB);

    // int x;
    // std::cin>>x;
    // certificate();
}