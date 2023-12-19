# pragma once
#include"config.h"

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

void getCH(std::vector<Line>& lines, std::vector<Line>& ch, std::vector<double>& pts) {
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

Array2d linearCHIntersect(const std::vector<Line>& ch1, const std::vector<Line>& ch2, 
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
