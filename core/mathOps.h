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
void getCH(std::vector<Line> lines, std::vector<Line>& ch, std::vector<double>& pts,
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

// void getCH(const std::map<double,double>& lines, std::map<double,double>& ch, std::vector<double>& pts,
// 			 const bool getUpperCH = true, const double& upperT = DeltaT) {
// 	ch.clear();
// 	pts.clear();
// 	pts.push_back(0);
// 	auto it=lines.begin();
// 	ch.insert(*it);
// 	double intsctX = 0;
// 	it++;
// 	while(it != lines.end()){
// 		// std::cout<<id<<"  "<<pts.size()<<"\n";
// 		while(!ch.empty()){
// 			intsctX = -(it->second-ch.rbegin()->second)/(it->first-ch.rbegin()->first);
// 			if(intsctX<=pts.back()){
// 			// if(ch.back().k*pts.back()+ch.back().b<=lines[id].k*pts.back()+lines[id].b){
// 				pts.pop_back();
// 				ch.erase(std::prev(ch.end()));
// 			}
// 			else break;
// 		}
// 		ch.insert(*it);
// 		pts.push_back(std::max(0.,intsctX));
// 		it++;
// 	}
// 	while(pts.back()>=upperT){
// 		pts.pop_back();
// 		ch.erase(std::prev(ch.end()));
// 	}
// 	pts.push_back(upperT);
// 	if(ch.empty()){
// 		std::cout<<"empty CH!\n";
// 		exit(-1);
// 	}
// 	if(ch.size()+1!=pts.size()){
// 		std::cout<<"segments and inflections are not compatible!\n"<<ch.size()<<" "<<pts.size();
// 		exit(-1);
// 	}
// }

Array2d linearCHIntersect(const std::vector<Line>& ch1, const std::vector<Line>& ch2, 
							const std::vector<double>& pts1, const std::vector<double>& pts2,
							const double& upperT = DeltaT) {
	int id1=0, id2=0;
	double intvL=-1, intvR=-1;
	double sweep=0, lastsweep=0;
	// bool stopInAdv = false; // 还没写上，但是比如：intvL==-1&&ch1[id1-1].k>ch1[id1-1].k

	auto checkSweepLine = [&] (const int id1, const int id2) {
		double y1=ch1[id1].k*sweep+ch1[id1].b;
		double y2=ch2[id2].k*sweep+ch2[id2].b;
		if (intvL!=-1&&y1>y2){
			intvR = lineIntersect_x(ch1[id1], ch2[id2]);
		} // 如果k1==k2，那么必然前一个节点必然已经满足y1<y2了
		else if(intvL==-1&&y1<y2){
			intvL = lineIntersect_x(ch1[id1], ch2[id2]);
		}// 如果k1==k2，那么必然前一个节点必然已经满足y1<y2了
		lastsweep = sweep;
	};

	if(ch1[0].b<ch2[0].b)
		intvL = 0;
	// else if(ch1[0].b==ch2[0].b)
	// 	intvL = intvR = 0;
	while(id1<ch1.size()&&id2<ch2.size()){
		sweep = std::min(pts1[id1+1], pts2[id2+1]);
		// if(ch1[id1].k!=ch2[id2].k&&sweep!=lastsweep){
		if(sweep!=lastsweep){
			// double y1=ch1[id1].k*sweep+ch1[id1].b;
			// double y2=ch2[id2].k*sweep+ch2[id2].b;
			// if (intvL!=-1&&y1>y2){
			// 	intvR = lineIntersect_x(ch1[id1], ch2[id2]);
			// 	return Array2d(intvL,intvR);
			// } // 如果k1==k2，那么必然前一个节点必然已经满足y1<y2了
			// else if(intvL==-1&&y1<y2){
			// 	intvL = lineIntersect_x(ch1[id1], ch2[id2]);
			// }// 如果k1==k2，那么必然前一个节点必然已经满足y1<y2了
			// lastsweep = sweep;
			checkSweepLine(id1,id2);
		}
		if(pts1[id1] < pts2[id2]) id1++;
		else id2++;
	}
	while(id1<ch1.size()){
		sweep = pts1[id1+1];
		// if(ch1[id1].k!=ch2[id2].k&&sweep!=lastsweep){
		if(sweep!=lastsweep){
			// double y1=ch1[id1].k*sweep+ch1[id1].b;
			// double y2=ch2[id2].k*sweep+ch2[id2].b;
			// if (intvL!=-1&&y1>y2){
			// 	intvR = lineIntersect_x(ch1[id1], ch2[id2]);
			// 	return Array2d(intvL,intvR);
			// } // 如果k1==k2，那么必然前一个节点必然已经满足y1<y2了
			// else if(intvL==-1&&y1<y2){
			// 	intvL = lineIntersect_x(ch1[id1], ch2[id2]);
			// }// 如果k1==k2，那么必然前一个节点必然已经满足y1<y2了
			// lastsweep = sweep;
			checkSweepLine(id1,id2);
		}
		id1++;
	}
	while(id2<ch2.size()){
		sweep = pts2[id2+1];
		// if(ch1[id1].k!=ch2[id2].k&&sweep!=lastsweep){
		if(sweep!=lastsweep){
			// double y1=ch1[id1].k*sweep+ch1[id1].b;
			// double y2=ch2[id2].k*sweep+ch2[id2].b;
			// if (intvL!=-1&&y1>y2){
			// 	intvR = lineIntersect_x(ch1[id1], ch2[id2]);
			// 	return Array2d(intvL,intvR);
			// } // 如果k1==k2，那么必然前一个节点必然已经满足y1<y2了
			// else if(intvL==-1&&y1<y2){
			// 	intvL = lineIntersect_x(ch1[id1], ch2[id2]);
			// }// 如果k1==k2，那么必然前一个节点必然已经满足y1<y2了
			// lastsweep = sweep;
			checkSweepLine(id1,id2);
		}
		id2++;
	}
	if(intvL!=-1 && intvR==-1)intvR=upperT;
	// std::cout<<intvL<<" "<<intvR<<"\n";
	return Array2d(intvL,intvR);
}

// Array2d linearCHIntersect(const std::vector<Line>& ch1, const std::vector<Line>& ch2, 
// 							const std::vector<double>& pts1, const std::vector<double>& pts2,
// 							const double& upperT = DeltaT) {
// 	int id1=0, id2=0;
// 	double intvL=-1, intvR=-1;
// 	std::vector<double> pts(pts1.size()+pts2.size());
// 	std::merge(pts1.begin(), pts1.end(), pts2.begin(), pts2.end(), pts.begin());
// 	pts.erase(std::unique(pts.begin(), pts.end()), pts.end()); // 去重
// 	pts.erase(pts.begin());
// 	// for(const auto& pt: pts){std::cout<<pt<<"\n";}

// 	if(ch1[0].b<ch2[0].b)
// 		intvL = 0;

// 	for(const auto& pt: pts){
// 		// std::cout<<pt<<" "<<id1<<" "<<id2<<"\n";
// 		if(intvL==-1&&ch1[id1].k>ch1[id1].k)break;
// 		if(pt>pts1[id1+1]&&id1<ch1.size()-1)id1++;
// 		if(pt>pts2[id1+1]&&id2<ch2.size()-1)id2++;
// 		double y1=ch1[id1].k*pt+ch1[id1].b;
// 		double y2=ch2[id2].k*pt+ch2[id2].b;
// 		if(y1==y2){
// 			if(intvL==-1)intvL=pt;
// 			else intvR=std::max(intvL,pt);
// 			// std::cout<<ch1[id1].k<<" "<<pt<<" "<<ch1[id1].b<<" "<<y1<<"\n";
// 			// std::cout<<ch2[id2].k<<" "<<pt<<" "<<ch2[id2].b<<" "<<y2<<"\n";
// 			// std::cerr<<"not implemented!\n";
// 			// exit(-1);
// 		}
// 		else if (intvL!=-1&&y1>y2){
// 			intvR = std::max(intvL,lineIntersect_x(ch1[id1], ch2[id2]));
// 			break;
// 		} // 如果k1==k2，那么必然前一个节点必然已经满足y1<y2了
// 		else if(intvL==-1&&y1<y2){
// 			intvL = lineIntersect_x(ch1[id1], ch2[id2]);
// 		}// 如果k1==k2，那么必然前一个节点必然已经满足y1<y2了
// 	}

// 	if(intvL!=-1 && intvR==-1)intvR=upperT;
// 	// std::cout<<intvL<<" "<<intvR<<"\n";
// 	return Array2d(intvL,intvR);
// }
