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

template<typename ObjType1, typename ObjType2>
static double calcDist(const ObjType1& CpPos1, const ObjType1& CpVel1,
					const ObjType2& CpPos2, const ObjType2& CpVel2,
					const Array2d& uv1, const Array2d& uv2, const double& t){
	Vector3d const p1 = CpPos1.evaluatePatchPoint(uv1);
	Vector3d const v1 = CpVel1.evaluatePatchPoint(uv1);
	Vector3d const p2 = CpPos2.evaluatePatchPoint(uv2);
	Vector3d const v2 = CpVel2.evaluatePatchPoint(uv2);
	Vector3d const pt1=(v1*t+p1), pt2=(v2*t+p2);
	return (pt2-pt1).norm();
}
template<typename ObjType>
static double calcAAExtent(const std::array<Vector3d, ObjType::cntCp>& ptPos) {
	double d=0;
	for(int axis=0;axis<3;axis++){
		double maxv = ptPos[0][axis], minv=maxv;
		for(int i = 1; i < ObjType::cntCp; i++) {
			maxv=std::max(maxv, ptPos[i][axis]);
			minv=std::min(minv, ptPos[i][axis]);
		}
		d=std::max(d,maxv-minv);
	}
	return d;
}

template<typename ObjType>
static double calcAAExtent_vec(const std::array<Vector3d, ObjType::cntCp>& ptPos) {
	Vector3d maxv=Vector3d::Constant(-INFT), minv=Vector3d::Constant(INFT);
	for(const auto& p:ptPos){
		maxv.cwiseMax(p);
		minv.cwiseMax(p);
	}
	return (maxv-minv).maxCoeff();
}

template<typename ParamObj1, typename ParamObj2>
void setAxes(const std::array<Vector3d, ParamObj1::cntCp>& ptPos1, 
					const std::array<Vector3d, ParamObj2::cntCp>& ptPos2,
					std::vector<Vector3d>& axes){	
	if(bbType==BoundingBoxType::AABB){
		axes = {Vector3d::Unit(0), Vector3d::Unit(1), Vector3d::Unit(2)};
	}
	else if(bbType==BoundingBoxType::DOP14){
		axes = {Vector3d::Unit(0), Vector3d::Unit(1), Vector3d::Unit(2), 
				Vector3d(1,1,1).normalized(), Vector3d(-1,1,1).normalized(), Vector3d(-1,-1,1).normalized()};
	}
	else if(bbType==BoundingBoxType::OBB){
		Vector3d lu1 = ParamObj1::axisU(ptPos1).normalized();//u延展的方向
		Vector3d lv1 = ParamObj1::axisV(ptPos1);//v延展的方向
		lv1 = (lv1-lv1.dot(lu1)*lu1).eval();
		Vector3d ln1 = lu1.cross(lv1);

		Vector3d lu2 = ParamObj2::axisU(ptPos2).normalized();//u延展的方向
		Vector3d lv2 = ParamObj2::axisV(ptPos2);//v延展的方向
		lv2 = (lv2-lv2.dot(lu2)*lu2).eval();
		Vector3d ln2 = lu2.cross(lv2);

		axes = {lu1,lv1,ln1,lu2,lv2,ln2, 
			lu1.cross(lu2), lu1.cross(lv2), lu1.cross(ln2), 
			lv1.cross(lu2), lv1.cross(lv2), lv1.cross(ln2), 
			ln1.cross(lu2), ln1.cross(lv2), ln1.cross(ln2)};
	}
}
