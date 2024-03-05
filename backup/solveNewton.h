# pragma once
#include "mathOps.h"
#include "triBezier.h"
#include "recBezier.h"
#include "recRatBezier.h"
template<typename ParamObj1, typename ParamObj2, typename ParamBound1, typename ParamBound2>
class SolverNewton{
	struct Interval1d{
		Array2d intv;
		Interval1d(): intv(0,1) {}
		Interval1d(const Array2d& a): intv(a) {}
		Interval1d(const double& a1, const double& a2): intv(a1,a2) {}
		double operator[](int i) const { return i == 0 ? intv[0] : intv[1]; }
		double mid() const { return 0.5 * (intv[0] + intv[1]); }
		Interval1d operator*(const Interval1d& i) const {
			Vector4d mltp(intv[0] * i[0], intv[0] * i[1], intv[1] * i[0], intv[1] * i[1]);
			return Interval1d(mltp.minCoeff(), mltp.maxCoeff());
		}
		Interval1d operator+(const Interval1d& i) const {
			return Interval1d(intv[0] + i[0], intv[1] + i[1]);
		}
		Interval1d operator-(const Interval1d& i) const {
			return Interval1d(intv[0] - i[1], intv[1] - i[0]);
		}
		Interval1d operator/(const Interval1d& i) const {
			if(i[0] <= 0 && i[1] >= 0){
				std::cerr<<"denominator intv includes 0!\n";
				exit(-1);
			}
			Interval1d intvInverse(1/i[1], 1/i[0]);
			return (*this) * intvInverse;
		}
	};
	template<int dim>
	struct IntervalXd{
		std::array<Interval1d, dim> intvs;
		Interval1d operator[](const int i) const { 
			if(i < 0 || i >= dim){
				std::cerr<<"IntervalXd out of range!\n";
				exit(-1);
			}
			return intvs[i];
		}
	};
	struct PatchPair{
		ParamBound1 pb1;
		ParamBound2 pb2;
		Array2d tIntv;
		PatchPair(const ParamBound1& c1, const ParamBound2& c2, 
				Array2d t = Array2d(0,DeltaT)): pb1(c1), pb2(c2), tIntv(t) {}
		bool operator<(PatchPair const &o) const { return tIntv[0] > o.tIntv[0]; }
		double calcL1Dist(const ParamObj1 &CpPos1, const ParamObj1 &CpVel1, 
						const ParamObj2 &CpPos2, const ParamObj2 &CpVel2) const{
			auto const ptPos1 = CpPos1.divideBezierPatch(pb1);
			auto const ptPos2 = CpPos2.divideBezierPatch(pb2);
			double d1=calcAAExtent<ParamObj1>(ptPos1);
			double d2=calcAAExtent<ParamObj1>(ptPos2);
			return std::max(d1, d2);
		}
		void bisectInterval(const ParamObj1 &CpPos1, const ParamObj1 &CpVel1, 
						const ParamObj2 &CpPos2, const ParamObj2 &CpVel2, 
						PatchPair& subp1, PatchPair& subp2) const{
			//算一下
			
		}
	};
public:
	static bool primitiveCheck(const ParamObj1 &CpPos1, const ParamObj1 &CpVel1, 
							const ParamObj2 &CpPos2, const ParamObj2 &CpVel2,
							const PatchPair& patchPair) {
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

		std::vector<Vector3d> axes;
		axes.clear();
		setAxes<ParamObj1, ParamObj2>(posStart1, posStart2, axes);

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
			if(maxProj2<minProj1 || maxProj1<minProj2) return false;
		}
		return true;
	}
	static double solveCCD(const ParamObj1 &CpPos1, const ParamObj1 &CpVel1, 
						const ParamObj2 &CpPos2, const ParamObj2 &CpVel2,
						Array2d& uv1, Array2d& uv2, 
						const double upperTime = DeltaT) {
		using steady_clock = std::chrono::steady_clock;
		using duration = std::chrono::duration<double>;
		const auto initialTime = steady_clock::now();

		std::priority_queue<PatchPair> heap;
		ParamBound1 initParam1;
		ParamBound2 initParam2;
		if (primitiveCheck(CpPos1, CpVel1, CpPos2, CpVel2, PatchPair(initParam1, initParam2)))
			heap.emplace(initParam1, initParam2);
		// cnt=1;
		while (!heap.empty()) {
			auto const cur = heap.top();
			heap.pop();
			// cnt++;
			// if(SHOWANS) std::cout<<cnt<<"\n";
			if (cur.calcL1Dist(CpPos1, CpVel1, CpPos2, CpVel2) < MinL1Dist) {
				uv1 = cur.pb1.centerParam();
				uv2 = cur.pb2.centerParam();
				const auto endTime = steady_clock::now();
				if(SHOWANS)
					std::cout << "min time: "<<  cur.tIntv[0] 
						<< "\nused seconds: " << duration(endTime - initialTime).count()
						<< std::endl;
				return cur.tIntv[0];
			}

			// Bisect the current param interval
			PatchPair subp1(cur), subp2(cur);
			cur.bisectInterval(CpPos1, CpVel1, CpPos2, CpVel2, subp1, subp2);
			if (primitiveCheck(CpPos1, CpVel1, CpPos2, CpVel2, subp1)){
				heap.push(subp1);
			}
			if (primitiveCheck(CpPos1, CpVel1, CpPos2, CpVel2, subp2)){
				heap.push(subp2);
			}
		}

		const auto endTime = steady_clock::now();
		if(SHOWANS)
			std::cout << "used seconds: " <<
				duration(endTime - initialTime).count()
				<< std::endl;
		return -1;
	}
};