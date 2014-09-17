#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_MESHER_AABB_H
#define ODER_MESHER_AABB_H

#include "oder.h"
#include "latool.h"

namespace ODER{
	template<class FT> struct AABB{
		typedef VectorBase<FT> Point;
		AABB(){
			pMin = Point(FLT_MAX, FLT_MAX, FLT_MAX);
			pMax = Point(-FLT_MAX, -FLT_MAX, -FLT_MAX);
		}
		AABB(const Point& p) :pMin(p), pMax(p){}
		AABB(const Point& p1, const Point& p2){
			pMin = Point(min(p1.x, p2.x), min(p1.y, p2.y), min(p1.z, p2.z));
			pMax = Point(max(p1.x, p2.x), max(p1.y, p2.y), max(p1.z, p2.z));
		}

		void Insert(const Point &p){
			pMin.x = min(pMin.x, p.x);
			pMin.y = min(pMin.y, p.y);
			pMin.z = min(pMin.z, p.z);
			pMax.x = max(pMax.x, p.x);
			pMax.y = max(pMax.y, p.y);
			pMax.z = max(pMax.z, p.z);
		}

		void Insert(const Vector& p){
			pMin.x = min(pMin.x, FT(p.x));
			pMin.y = min(pMin.y, FT(p.y));
			pMin.z = min(pMin.z, FT(p.z));
			pMax.x = max(pMax.x, FT(p.x));
			pMax.y = max(pMax.y, FT(p.y));
			pMax.z = max(pMax.z, FT(p.z));
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
		void Expand(float delta){
			VectorBase<FT> v = VectorBase<FT>(delta, delta, delta);
			pMin -= v;
			pMax += v;
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