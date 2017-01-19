#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_MESHER_GEOMETER_H
#define ODER_MESHER_GEOMETER_H

#include "oder.h"
#include "latool.h"
#include "numerMethod.h"
#include "predicate.h"
#include <random>

namespace ODER {
	namespace Geometer {
		template<class FT> void Circumcircle(const VectorBase<FT>& a, const VectorBase<FT>& b, const VectorBase<FT>& c, 
			VectorBase<FT> *center, FT *r = NULL){
			VectorBase<FT> ca = a - c;
			VectorBase<FT> cb = b - c;
			VectorBase<FT> n = ca%cb;

			FT A[9];
			A[0] = ca[0], A[1] = ca[1], A[2] = ca[2],
			A[3] = cb[0], A[4] = cb[1], A[5] = cb[2],
			A[6] = n[0], A[7] = n[1], A[8] = n[2];

			FT rhs[3];
			rhs[0] = (ca.length2())*FT(0.5);
			rhs[1] = (cb.length2())*FT(0.5);
			rhs[2] = FT(0);

			VectorBase<FT> rr;
			gaussianElimination3x3(A, rhs, &rr[0]);

			if (center)
				*center = c + rr;

			if (r)
				*r = sqrt(rr.length2());
		}

		template<class FT> void Circumsphere(const VectorBase<FT>& a, const VectorBase<FT>& b, const VectorBase<FT>& c, const VectorBase<FT>& d, 
			VectorBase<FT> *center, FT *r = NULL){
			VectorBase<FT> t = a - d;
			VectorBase<FT> u = b - d;
			VectorBase<FT> v = c - d;

			FT A[9];
			A[0] = t[0], A[1] = t[1], A[2] = t[2],
			A[3] = u[0], A[4] = u[1], A[5] = u[2],
			A[6] = v[0], A[7] = v[1], A[8] = v[2];

			FT rhs[3];
			rhs[0] = t.length2()*FT(0.5);
			rhs[1] = u.length2()*FT(0.5);
			rhs[2] = v.length2()*FT(0.5);

			VectorBase<FT> rr;
			gaussianElimination3x3(A, rhs, &rr[0]);

			if (center)
				*center = d + rr;

			if (r)
				*r = sqrt(rr.length2());
		}

		template<class FT> void Orthocircle(const VectorBase<FT>& a, FT aWeight, const VectorBase<FT>& b, FT bWeight,
			const VectorBase<FT>& c, FT cWeight, VectorBase<FT> *center, FT *r = NULL) {
			VectorBase<FT> ca = a - c;
			VectorBase<FT> cb = b - c;
			VectorBase<FT> n = ca%cb;

			FT A[9];
			A[0] = ca[0], A[1] = ca[1], A[2] = ca[2],
			A[3] = cb[0], A[4] = cb[1], A[5] = cb[2],
			A[6] = n[0], A[7] = n[1], A[8] = n[2];

			FT rhs[3];
			rhs[0] = (ca.length2() + cWeight - aWeight)*FT(0.5);
			rhs[1] = (cb.length2() + cWeight - bWeight)*FT(0.5);
			rhs[2] = FT(0);

			VectorBase<FT> rr;
			gaussianElimination3x3(A, rhs, &rr[0]);

			if (center)
				*center = c + rr;

			if (r)
				*r = sqrt(rr.length2() - cWeight);
		}

		template<class FT> void Orthosphere(const VectorBase<FT>& a, FT aWeight, const VectorBase<FT>& b, FT bWeight,
			const VectorBase<FT>& c, FT cWeight, const VectorBase<FT>& d, FT dWeight,
			VectorBase<FT> *center, FT *r = NULL) {

			VectorBase<FT> t = a - d;
			VectorBase<FT> u = b - d;
			VectorBase<FT> v = c - d;

			FT A[9];
			A[0] = t[0], A[1] = t[1], A[2] = t[2],
			A[3] = u[0], A[4] = u[1], A[5] = u[2],
			A[6] = v[0], A[7] = v[1], A[8] = v[2];

			FT rhs[3];
			rhs[0] = (t.length2() + dWeight - aWeight)*FT(0.5);
			rhs[1] = (u.length2() + dWeight - bWeight)*FT(0.5);
			rhs[2] = (v.length2() + dWeight - cWeight)*FT(0.5);

			VectorBase<FT> rr;
			gaussianElimination3x3(A, rhs, &rr[0]);

			if (center)
				*center = d + rr;

			if (r)
				*r = sqrt(rr.length2() - dWeight);
		}

		template<class FT> void Orthocenter(const VectorBase<FT>& a, FT aWeight, const VectorBase<FT>& b, FT bWeight, 
			VectorBase<FT> *center, FT *r = NULL) {
			VectorBase<FT> ab = b - a;
			FT abLen = ab.length();
			FT extraDis = (aWeight - bWeight) / (FT(2.0) * abLen);

			if (center)
				*center = a + (FT(0.5) + extraDis / abLen)*ab;

			if (r) {
				FT distance = abLen*FT(0.5) + extraDis;
				*r = sqrt(distance*distance - aWeight);
			}
		}


		template<class FT> VectorBase<FT> triangleNormal(const VectorBase<FT>& a,
			const VectorBase<FT>& b, const VectorBase<FT>& c) {
			VectorBase<FT> ab = b - a;
			VectorBase<FT> ac = c - a;
			VectorBase<FT> bc = c - b;
			FT abLen2 = ab.length2();
			FT acLen2 = ac.length2();
			FT bcLen2 = bc.length2();

			VectorBase<FT> *u = &ab;
			VectorBase<FT> *v = &ac;

			if (bcLen2 < abLen2) {
				if (abLen2 < acLen2) v = &bc;
				else {
					u = &ac; v = &bc;
				}
			}
			else if (bcLen2 < acLen2) v = &bc;

			return *u % *v;
		}

		template<class FT> VectorBase<FT> calculateAbovePoint(const VectorBase<FT>& a,
			const VectorBase<FT>& b, const VectorBase<FT>& c) {
			VectorBase<FT> ab = b - a;
			VectorBase<FT> ac = c - a;
			VectorBase<FT> bc = c - b;
			FT abLen2 = ab.length2();
			FT acLen2 = ac.length2();
			FT bcLen2 = bc.length2();
			FT maxLen2 = bcLen2;

			VectorBase<FT> *u = &ab;
			VectorBase<FT> *v = &ac;
			VectorBase<FT> corner = a;

			if (bcLen2 < abLen2) {
				if (abLen2 < acLen2) {
					v = &bc; 
					maxLen2 = acLen2;
				}
				else {
					u = &ac; v = &bc; 
					maxLen2 = abLen2;
				}
			}
			else if (bcLen2 < acLen2) {
				v = &bc;
				maxLen2 = acLen2;
			}
			else corner = b;

			VectorBase<FT> n = *u % *v;
			FT nLen2 = n.length2();
			if (nLen2 == FT(0)) {
				constexpr Predicator<FT> predicator;
				n.x = predicator.orient2d(a.y, a.z, b.y, b.z, c.y, c.z);
				n.y = predicator.orient2d(a.z, a.x, b.z, b.x, c.z, c.x);
				n.z = predicator.orient2d(a.x, a.y, b.x, b.y, c.x, c.y);
				nLen2 = n.length2();
			}
			FT scale = sqrt(maxLen2) / sqrt(nLen2);
			return corner + (scale * n);
		}

		template<class FT> inline FT interiorAngle(const VectorBase<FT>& o, const VectorBase<FT>& a, const VectorBase<FT>& b) {
			constexpr FT pi = M_PI;
			VectorBase<FT> oa = a - o;
			VectorBase<FT> ob = b - o;

			FT len = oa.length() * ob.length();
			Assert(len != FT(0));

			return acos(Clamp((oa * ob) / len, FT(-1), FT(1)));
		}

		//a and b are common edge vertices
		template<class FT> FT dihedralAngle(const VectorBase<FT>& a, const VectorBase<FT>& b,
			const VectorBase<FT>& c, const VectorBase<FT>& d) {
			constexpr FT pi = M_PI;
			VectorBase<FT> nc = triangleNormal(a, b, d);
			VectorBase<FT> nd = triangleNormal(b, a, c);
			FT abLen = (b - a).length();

			constexpr Predicator<FT> predicator;
			FT ori = predicator.orient3d(d, a, b, c);

			FT theta = pi * FT(0.5);
			FT dot = nc * nd;
			if (dot != 0) {
				theta = atan(-ori * abLen / dot);
				if (dot > 0) theta += (ori > 0 ? pi : -pi);
				theta = fabs(theta);
			}

			return theta;
		}

		template<class FT> VectorBase<FT> getProjectedPoint(const VectorBase<FT>& a, const VectorBase<FT>& b,
			const VectorBase<FT>& c, const VectorBase<FT>& p) {
			VectorBase<FT> normal = Normalize(triangleNormal(a, b, c));
			VectorBase<FT> pa = a - p;
			return p + (pa * normal) * normal;
		}

		namespace HilbertSortInternal {
			template<int x, bool xFlip, bool yFlip, bool zFlip, size_t threshold,
				class RandomAccessIterator, class CoordTrait>
			void HilbertSort3dRecursive(RandomAccessIterator begin, RandomAccessIterator end, CoordTrait ct) {
				constexpr int y = (x + 1) % 3, z = (x + 2) % 3;
				if (end - begin <= threshold) return;

				auto CompXLess = [ct](auto p, auto q){ return xFlip ? ct(q).x < ct(p).x : ct(p).x < ct(q).x; };
				auto CompYLess = [ct](auto p, auto q){ return yFlip ? ct(q).y < ct(p).y : ct(p).y < ct(q).y; };
				auto CompZLess = [ct](auto p, auto q){ return zFlip ? ct(q).z < ct(p).z : ct(p).z < ct(q).z; };
				auto CompXGreat = [ct](auto p, auto q){ return xFlip ? ct(q).x > ct(p).x : ct(p).x > ct(q).x; };
				auto CompYGreat = [ct](auto p, auto q){ return yFlip ? ct(q).y > ct(p).y : ct(p).y > ct(q).y; };
				auto CompZGreat = [ct](auto p, auto q){ return zFlip ? ct(q).z > ct(p).z : ct(p).z > ct(q).z; };

				auto HilbertSplit = [](auto begin, auto end, auto comp) {
					auto mid = begin + (end - begin) / 2;
					std::nth_element(begin, mid, end, comp);
					return mid;
				};

				RandomAccessIterator iter0 = begin, iter8 = end;
				RandomAccessIterator iter4 = HilbertSplit(iter0, iter8, CompXLess);
				RandomAccessIterator iter2 = HilbertSplit(iter0, iter4, CompYLess);
				RandomAccessIterator iter1 = HilbertSplit(iter0, iter2, CompZLess);
				RandomAccessIterator iter3 = HilbertSplit(iter2, iter4, CompZGreat);
				RandomAccessIterator iter6 = HilbertSplit(iter4, iter8, CompYGreat);
				RandomAccessIterator iter5 = HilbertSplit(iter4, iter6, CompZLess);
				RandomAccessIterator iter7 = HilbertSplit(iter6, iter8, CompZGreat);

				HilbertSort3dRecursive<z, zFlip, xFlip, yFlip, threshold>(iter0, iter1, ct);
				HilbertSort3dRecursive<y, yFlip, zFlip, xFlip, threshold>(iter1, iter2, ct);
				HilbertSort3dRecursive<y, yFlip, zFlip, xFlip, threshold>(iter2, iter3, ct);
				HilbertSort3dRecursive<x, xFlip, !yFlip, !zFlip, threshold>(iter3, iter4, ct);
				HilbertSort3dRecursive<x, xFlip, !yFlip, !zFlip, threshold>(iter4, iter5, ct);
				HilbertSort3dRecursive<y, !yFlip, zFlip, !xFlip, threshold>(iter5, iter6, ct);
				HilbertSort3dRecursive<y, !yFlip, zFlip, !xFlip, threshold>(iter6, iter7, ct);
				HilbertSort3dRecursive<z, !zFlip, !xFlip, yFlip, threshold>(iter7, iter8, ct);
			}
		}

		template<bool reverse = false, size_t threshold = 1u, class RandomAccessIterator, class CoordTrait>
		void HilbertSort3d(RandomAccessIterator begin, RandomAccessIterator end, CoordTrait ct) {
			static_assert(std::is_same<std::iterator_traits<RandomAccessIterator>::iterator_category,
				std::random_access_iterator_tag>::value, "ODER::HilbertSort3d support random access iterators only");
			HilbertSortInternal::HilbertSort3dRecursive<0, reverse, false, false, threshold>(begin, end, ct);
		}
		template<bool reverse = false, size_t threshold = 1u, class RandomAccessIterator>
		void HilbertSort3d(RandomAccessIterator begin, RandomAccessIterator end) {
			return HilbertSort3d<reverse, threshold>(begin, end, [](auto p){ return p; });
		}

		template<class RandomAccessIterator, class CoordTrait> 
		void brioSort3d(RandomAccessIterator begin, RandomAccessIterator end, CoordTrait ct) {
			static_assert(std::is_same<std::iterator_traits<RandomAccessIterator>::iterator_category,
				std::random_access_iterator_tag>::value, "ODER::brioSort3d support random access iterators only");

			using difference_type = std::iterator_traits<RandomAccessIterator>::difference_type;
			std::random_device rng;
			std::default_random_engine randomEngine(rng());
			std::shuffle(begin, end, randomEngine);

			constexpr double ratio = 0.5;
			constexpr difference_type brioThreshold = 64u;
			constexpr size_t HilbertThreshold = 8u;
			difference_type range = end - begin;
			size_t roundCount = size_t(ceil(log2(range) * log2(1.0 / ratio)));

			RandomAccessIterator curEnd = end;
			size_t round = 0u;
			for (; round < roundCount && range >= brioThreshold; round++) {
				std::binomial_distribution<difference_type> distrbute(range, 1.0 - ratio);
				range = distrbute(randomEngine);
				RandomAccessIterator curMid = begin + range;

				if (round % 2u == 0u) HilbertSort3d<false, HilbertThreshold>(curMid, curEnd, ct);
				else HilbertSort3d<true, HilbertThreshold>(curMid, curEnd, ct);

				curEnd = curMid;
				std::shuffle(begin, curEnd, randomEngine);
			}
			if (round % 2u == 0u) HilbertSort3d<false, HilbertThreshold>(begin, curEnd, ct);
			else HilbertSort3d<false, HilbertThreshold>(begin, curEnd, ct);
		}

		template<class RandomAccessIterator>
		void brioSort3d(RandomAccessIterator begin, RandomAccessIterator end) {
			return brioSort3d(begin, end, [](auto p){ return p; });
		}
	}
}

#endif
