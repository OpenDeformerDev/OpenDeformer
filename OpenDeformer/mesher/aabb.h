#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_MESHER_AABB_H
#define ODER_MESHER_AABB_H

#include "oder.h"
#include "latool.h"
#include <numeric>

namespace ODER{
	template<class FT> struct AABB{
		using Point = VectorBase<FT>;
		AABB(){
			constexpr FT max = std::numeric_limits<FT>::max();
			pMin = Point(max, max, max);
			pMax = Point(-max, -max, -max);
		}
		AABB(const Point& p) :pMin(p), pMax(p){}
		AABB(const Point& p1, const Point& p2){
			pMin = Point(std::min(p1.x, p2.x), std::min(p1.y, p2.y), std::min(p1.z, p2.z));
			pMax = Point(std::max(p1.x, p2.x), std::max(p1.y, p2.y), std::max(p1.z, p2.z));
		}

		template<class T> void Insert(const VectorBase<T> &p){
			pMin.x = std::min(pMin.x, FT(p.x));
			pMin.y = std::min(pMin.y, FT(p.y));
			pMin.z = std::min(pMin.z, FT(p.z));
			pMax.x = std::max(pMax.x, FT(p.x));
			pMax.y = std::max(pMax.y, FT(p.y));
			pMax.z = std::max(pMax.z, FT(p.z));
		}

		template<> void Insert<FT>(const Point &p){
			pMin.x = std::min(pMin.x, p.x);
			pMin.y = std::min(pMin.y, p.y);
			pMin.z = std::min(pMin.z, p.z);
			pMax.x = std::max(pMax.x, p.x);
			pMax.y = std::max(pMax.y, p.y);
			pMax.z = std::max(pMax.z, p.z);
		}

		bool Overlap(const AABB& b){
			if (pMax.x<b.pMin.x || pMin.x>b.pMax.x) return false;
			if (pMax.y<b.pMin.y || pMin.y>b.pMax.y) return false;
			if (pMax.z<b.pMin.z || pMin.z>b.pMax.z) return false;

			return true;
		}

		bool Inside(const Point& p) const{
			return (p.x >= pMin.x && p.x <= pMax.x &&
				p.y >= pMin.y && p.y <= pMax.y &&
				p.z >= pMin.z && p.z <= pMax.z);
		}

		void Expand(FT delta){
			VectorBase<FT> v = VectorBase<FT>(delta, delta, delta);
			pMin -= v;
			pMax += v;
		}

		void Scale(FT factor) {
			VectorBase<FT> d = (factor * FT(0.5)) * (pMax - pMin);
			VectorBase<FT> mid = FT(0.5) * (pMax + pMin);
			pMin = mid - d;
			pMax = mid + d;
		}

		int maxExtent() const{
			VectorBase<FT> d = pMax - pMin;
			return d.x > d.y ? (d.x > d.z ? 0 : 2) : (d.y > d.z ? 1 : 2);
		}
	private:
		Point pMin, pMax;
	};
}

#endif