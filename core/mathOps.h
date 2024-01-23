# pragma once
#include"config.h"

struct Line
{
	double k, b;
	Line():k(0),b(0){}
	void set(const double& _k, const double &_b){k=_k,b=_b;}
	Line(const double& k,const double& b): k(k), b(b) {}
	bool operator<(const Line &l) const { 
		return k < l.k || (k == l.k && b > l.b); // 相同斜率的直线中只有截距最大的被留下来
	}
	bool operator==(const Line &l) const {return k == l.k;}
};
inline double lineIntersect_x(const Line& l1, const Line &l2) {
	if(l1.k==l2.k){
				std::cout<<"parallel lines do not intersect at a single point!\n";
		exit(-1);
	}
	return (l2.b-l1.b)/(l1.k-l2.k);
}
void getCH(std::vector<Line>& lines, std::vector<Line>& ch, std::vector<double>& pts,
			 const bool getUpperCH = true, const double& upperT = DeltaT) {
	if(!getUpperCH)std::reverse(lines.begin(),lines.end());
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
			intsctX = lineIntersect_x(lines[id], ch.back());
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
	while(pts.back()>=upperT){
		pts.pop_back();
		ch.pop_back();
	}
	pts.push_back(upperT);
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
// void getCH(std::vector<Line>& lines, std::vector<Line>& ch, std::vector<double>& pts,
// 			 const bool getUpperCH = true, const double& upperT = DeltaT) {
// 	if(!getUpperCH)std::reverse(lines.begin(),lines.end());
// 	lines.erase(std::unique(lines.begin(), lines.end()), lines.end()); // 去重
// 	ch.clear();
// 	pts.clear();
// 	pts.push_back(0);
// 	ch.push_back(lines[0]);
// 	std::vector<double> heights; heights.clear();heights.push_back(lines[0].b);
// 	if(DEBUG) for(const auto&l:lines)std::cout<<"lines:\t"<<l.k<<" "<<l.b<<"\n";
// 	double id = 1, height = 0;
// 	while(id < lines.size()){
// 		// std::cout<<id<<"  "<<pts.size()<<"\n";
// 		while(!ch.empty()){
// 			height=lines[id].k*pts.back()+lines[id].b;
// 			if((getUpperCH&&height>heights.back())||(!getUpperCH&&height<heights.back())){
// 			// if(ch.back().k*pts.back()+ch.back().b<=lines[id].k*pts.back()+lines[id].b){
// 				pts.pop_back();
// 				ch.pop_back();
// 				heights.pop_back();
// 			}
// 			else break;
// 		}
// 		// intsctX = lineIntersect_x(lines[id], ch.back());
// 		heights.push_back(height);
// 		pts.push_back(ch.empty()?0:std::max(0.,lineIntersect_x(lines[id], ch.back(), false)));
// 		ch.push_back(lines[id]);
// 		id++;
// 	}
// 	while(pts.back()>=upperT){
// 		pts.pop_back();
// 		ch.pop_back();
// 	}
// 	pts.push_back(upperT);
// 	if(DEBUG) for(const auto&l:ch)std::cout<<"ch:\t"<<l.k<<" "<<l.b<<"\n";
// 	if(DEBUG) for(const auto&pt:pts)std::cout<<"pt:\t"<<pt<<"\n";
// 	if(ch.empty()){
// 		std::cout<<"empty CH!\n";
// 		exit(-1);
// 	}
// 	if(ch.size()+1!=pts.size()){
// 		std::cout<<"segments and inflections are not compatible!\n";
// 		exit(-1);
// 	}
// }

Array2d linearCHIntersect(const std::vector<Line>& ch1, const std::vector<Line>& ch2, 
							const std::vector<double>& pts1, const std::vector<double>& pts2,
							const double& upperT = DeltaT) {
	int id1=1, id2=1;
	double intvL=-1, intvR=-1;
	double sweep=0, lastsweep=0;
	bool stopInAdv = false; // 还没写上，但是比如：intvL==-1&&ch1[id1-1].k>ch1[id1-1].k

	auto checkSweepLine = [&] (const int id1, const int id2) {
		double y1=ch1[id1].k*sweep+ch1[id1].b;
		double y2=ch2[id2].k*sweep+ch2[id2].b;
		if (y1>y2){
			if(intvL!=-1)
				intvR = lineIntersect_x(ch1[id1], ch2[id2]);
		} // 如果k1==k2，那么必然前一个节点必然已经满足y1<y2了
		else if(y1<y2){
			if(intvL==-1)
				intvL = lineIntersect_x(ch1[id1], ch2[id2]); 
		}// 如果k1==k2，那么必然前一个节点必然已经满足y1<y2了
		// else{
		// 	//y1==y2
		// 	// if(intvL==-1) intvR = intvL = sweep;//这个不对，如果intv一直持续到deltaT就会变成只有一个点
		// 	if(intvL!=-1) 
		// 		intvR = sweep;
		// }
		if (DEBUG) std::cout<<id1<<"  "<<id2<<" /  "<<y1<<"  "<<y2<<" /  "<<intvL<<" "<<intvR<<"\n";
		lastsweep = sweep;
	};

	if(ch1[0].b<ch2[0].b)
		intvL = 0;
	// else if(ch1[0].b==ch2[0].b)
	// 	intvL = intvR = 0;
	while(id1<pts1.size()&&id2<pts2.size()&&intvR==-1){
		sweep = std::min(pts1[id1], pts2[id2]);
		if(sweep!=lastsweep)//并不知道这样跳过能不能更快
			checkSweepLine(id1-1, id2-1);
		if(pts1[id1] < pts2[id2]) id1++;
		else id2++;
	}
	while(id1<pts1.size()&&intvR==-1){
		sweep = pts1[id1];
		if(sweep!=lastsweep)//并不知道这样跳过能不能更快
			checkSweepLine(id1-1, id2-1);
		id1++;
	}
	while(id2<pts2.size()&&intvR==-1){
		sweep = pts2[id2];
		if(sweep!=lastsweep)//并不知道这样跳过能不能更快axesOBB
			checkSweepLine(id1-1, id2-1);
		id2++;
	}
	if(intvL!=-1 && intvR==-1)intvR=upperT;
	return Array2d(intvL,intvR);
}

template<typename ObjType>
static void generatePatchPair(std::array<Vector3d, ObjType::cntCp> &CpPos1, std::array<Vector3d, ObjType::cntCp> &CpVel1,
							 std::array<Vector3d, ObjType::cntCp> &CpPos2, std::array<Vector3d, ObjType::cntCp> &CpVel2, const double& denom){
	Vector3d dir=Vector3d::Random().normalized()/denom;
	for (int i = 0; i < ObjType::cntCp; i++) {
		for(int dim=0; dim<3; dim++) CpPos1[i][dim] = randNormal(randGenerator);
		for(int dim=0; dim<3; dim++) CpVel1[i][dim] = randNormal(randGenerator);
		for(int dim=0; dim<3; dim++) CpPos2[i][dim] = randNormal(randGenerator);
		for(int dim=0; dim<3; dim++) CpVel2[i][dim] = randNormal(randGenerator);
		CpPos1[i]+=dir;
		CpVel1[i]-=dir;
		CpPos2[i]-=dir;
		CpVel2[i]+=dir;
	}
}

template<typename ObjType>
static void generatePatchPair(std::array<Vector4d, ObjType::cntCp> &CpPos1, std::array<Vector4d, ObjType::cntCp> &CpVel1,
							 std::array<Vector4d, ObjType::cntCp> &CpPos2, std::array<Vector4d, ObjType::cntCp> &CpVel2, const double& denom){
	Vector3d dir=Vector3d::Random().normalized()/denom;
	for (int i = 0; i < ObjType::cntCp; i++) {
		for(int dim=0; dim<4; dim++) CpPos1[i][dim] = randNormal(randGenerator);
		for(int dim=0; dim<4; dim++) CpVel1[i][dim] = randNormal(randGenerator);
		for(int dim=0; dim<4; dim++) CpPos2[i][dim] = randNormal(randGenerator);
		for(int dim=0; dim<4; dim++) CpVel2[i][dim] = randNormal(randGenerator);
		CpPos1[i].segment(0,3)+=dir*CpPos1[i][3];
		CpVel1[i].segment(0,3)+=dir*CpVel1[i][3];
		CpPos2[i].segment(0,3)+=dir*CpPos2[i][3];
		CpVel2[i].segment(0,3)+=dir*CpVel2[i][3];
	}
}