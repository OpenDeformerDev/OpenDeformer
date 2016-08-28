#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_CORE_ARITHMETRIC_H
#define ODER_CORE_ARITHMETRIC_H

#include "oder.h"
#include <cstdint>

#define ABOSOLUTE_GREATER(x, y) ((x)>(y)) == ((x)>(-y)) //returns ture if |x| > |y|

#if defined(_MSC_VER)
#pragma float_control(precise, on, push)
#pragma fp_contract(off)
#endif

namespace ODER{
	template<class FT> class ExactArthmeticer{
	public:
		FT Estimate(FT *e, int n) const{
			FT ret = e[0];
			for (int i = 1; i < n; i++){
				ret += e[i];
			}
			return ret;
		}

		//|a| >= |b| must be confirmed
		void fastTwoSum(FT a, FT b, FT &esti, FT& err) const {
			esti = a + b;
			err = b - (esti - a);
		}

		//|a| >= |b| must be confirmed
		FT fastTwoSumError(FT a, FT b, FT esti) const{
			return b - (esti - a);
		}

		//|a| >= |b| must be confirmed
		void fastTwoDiff(FT a, FT b, FT& esti, FT& err) const {
			esti = a - b;
			err = (a - esti) - b;
		}

		//|a| >= |b| must be confirmed
		void fastTwoDiffError(FT a, FT b, FT esti) const {
			return (a - esti) - b;
		}

		void twoSum(FT a, FT b, FT &esti, FT& err) const {
			esti = a + b;
			FT bvirt = esti - a;
			FT avirt = esti - bvirt;
			err = (b - bvirt) + (a - avirt);
		}

		FT twoSumError(FT a, FT b, FT esti) const{
			FT bvirt = esti - a;
			FT avirt = esti - bvirt;
			return (b - bvirt) + (a - avirt);
		}

		void twoDiff(FT a, FT b, FT &esti, FT& err) const{
			esti = a - b;
			FT bvirt = a - esti;
			FT avirt = esti + bvirt;
			err = (bvirt - b) + (a - avirt);
		}

		FT twoDiffError(FT a, FT b, FT esti) const{
			FT bvirt = a - esti;
			FT avirt = esti + bvirt;
			return (bvirt - b) + (a - avirt);
		}

		void Split(FT a, FT &high, FT& low) const {
			FT c = splitter*a;
			FT big = c - a;
			high = c - big;
			low = a - high;
		}

		void twoProduct(FT a, FT b, FT& esti, FT& err) const {
			esti = a * b;
			FT aHigh, bHigh;
			FT aLow, bLow;
			Split(a, aHigh, aLow);
			Split(b, bHigh, bLow);
			FT e0 = esti - (aHigh*bHigh);
			FT e1 = e0 - (aLow*bHigh);
			FT e2 = e1 - (aHigh*bLow);
			err = (aLow*bLow) - e2;
		}

		void twoProductPreSplit(FT a, FT b, FT bHigh, FT bLow, FT& esti, FT& err) const{
			esti = a * b;
			FT aHigh, aLow;
			Split(a, aHigh, aLow);
			FT e0 = esti - (aHigh*bHigh);
			FT e1 = e0 - (aLow*bHigh);
			FT e2 = e1 - (aHigh*bLow);
			err = (aLow*bLow) - e2;
		}

		void twoOneSum(FT a1, FT a0, FT b, FT& x2, FT& x1, FT& x0) const {
			FT iq;
			//a1 + a0 + b
			twoSum(a0, b, iq, x0);
			twoSum(a1, iq, x2, x1);
		}

		void twoTwoSum(FT a1, FT a0, FT b1, FT b0, FT& x3, FT& x2, FT& x1, FT& x0) const {
			FT i, j;
			twoOneSum(a1, a0, b0, i, j, x0);
			twoOneSum(i, j, b1, x3, x2, x1);
		}

		void fourOneSum(FT a3, FT a2, FT a1, FT a0, FT b, FT& x4, FT& x3, FT& x2, FT& x1, FT& x0) const {
			FT iq;
			twoOneSum(a1, a0, b, iq, x1, x0);
			twoOneSum(a3, a2, iq, x4, x3, x2);
		}

		void fourTwoSum(FT a3, FT a2, FT a1, FT a0, FT b1, FT b0,
			FT& x5, FT& x4, FT& x3, FT& x2, FT& x1, FT& x0) const {
			FT i, j, k, l;
			fourOneSum(a3, a2, a1, a0, b0, i, j, k, l, x0);
			fourOneSum(i, j, k, l, b1, x5, x4, x3, x2, x1);
		}

		void twoTwoTwoSum(FT a1, FT a0, FT b1, FT b0, FT c1, FT c0,
			FT& x5, FT& x4, FT& x3, FT& x2, FT& x1, FT& x0) const {
			FT i, j, k, l;
			twoTwoSum(a1, a0, b1, b0, i, j, k, l);
			fourTwoSum(i, j, k, l, c1, c0, x5, x4, x3, x2, x1, x0);
		}

		void twoOneDiff(FT a1, FT a0, FT b, FT& x2, FT& x1, FT& x0) const{
			FT iq;
			//a1 + a0 - b
			twoDiff(a0, b, iq, x0);
			twoSum(a1, iq, x2, x1);
		}

		void twoTwoDiff(FT a1, FT a0, FT b1, FT b0, FT& x3, FT& x2, FT& x1, FT& x0) const{
			FT i, j;
			twoOneDiff(a1, a0, b0, i, j, x0);
			twoOneDiff(i, j, b1, x3, x2, x1);
		}

		void twoOneProduct(FT a1, FT a0, FT b, FT& x3, FT& x2, FT& x1, FT& x0) const{
			FT bHigh, bLow, i, j, k, l;
			Split(b, bHigh, bLow);
			twoProductPreSplit(a0, b, bHigh, bLow, i, x0);
			twoProductPreSplit(a1, b, bHigh, bLow, j, k);
			twoSum(i, k, l, x1);
			fastTwoSum(j, l, x3, x2);
		}

		int growExpansion(FT *e, FT b, FT *h, int n) const{
			FT q, qnext;
			q = b;
			for (int i = 0; i < n; i++){
				twoSum(q, e[i], qnext, h[i]);
				q = qnext;
			}
			h[n] = q;
			return n + 1;
		}

		int fastExpansionSum(FT *e, FT *f, FT *h, int em, int fn) const{
			FT q, qnext, hh;
			FT enext = e[0], fnext = f[0];
			int eIndex = 0, fIndex = 0, hIndex = 0;

			if (ABOSOLUTE_GREATER(fnext, enext)){
				q = enext;
				enext = e[++eIndex];
			}
			else{
				q = fnext;
				fnext = f[++fIndex];
			}
			if ((eIndex < em) && (fIndex < fn)){
				if (ABOSOLUTE_GREATER(fnext, enext)){
					fastTwoSum(enext, q, qnext, hh);
					enext = e[++eIndex];
				}
				else{
					fastTwoSum(fnext, q, qnext, hh);
					fnext = f[++fIndex];
				}
				q = qnext;
				if (hh != FT(0))
					h[hIndex++] = hh;

				while ((eIndex < em) && (fIndex < fn)){
					if (ABOSOLUTE_GREATER(fnext, enext)){
						twoSum(enext, q, qnext, hh);
						enext = e[++eIndex];
					}
					else{
						twoSum(fnext, q, qnext, hh);
						fnext = f[++fIndex];
					}
					q = qnext;
					if (hh != FT(0))
						h[hIndex++] = hh;
				}
			}
			while (eIndex < em){
				twoSum(enext, q, qnext, hh);
				enext = e[++eIndex];
				q = qnext;
				if (hh != FT(0))
					h[hIndex++] = hh;
			}
			while (fIndex < fn){
				twoSum(fnext, q, qnext, hh);
				fnext = f[++fIndex];
				q = qnext;
				if (hh != FT(0))
					h[hIndex++] = hh;
			}
			if ((q != FT(0)) || (hIndex == 0)){
				h[hIndex++] = q;
			}
			return hIndex;
		}

		int scaleExpansion(FT *e, FT b, FT *h, int n) const{
			FT q, qnext;
			FT hh;
			FT t0, t1;
			FT bLow, bHigh;
			int hi = 0;
			Split(b, bHigh, bLow);
			twoProductPreSplit(e[0], b, bHigh, bLow, q, hh);
			if (hh != 0.f)
				h[hi++] = hh;
			for (int i = 1; i < n; i++){
				twoProductPreSplit(e[i], b, bHigh, bLow, t1, t0);
				twoSum(q, t0, qnext, hh);
				if (hh != 0.f)
					h[hi++] = hh;
				fastTwoSum(t1, qnext, q, hh);
				if (hh != 0.f)
					h[hi++] = hh;
			}
			if ((q != 0.0) || (hi == 0))
				h[hi++] = q;
			return hi;
		}

		private:
			template<class FT> static constexpr FT getEpsilon() noexcept{
				static_assert(false, "ODER::ExactArthmeticer support IEEE 754-1985 floating point only");
				return 0;
			}
			template<> static constexpr float getEpsilon<float>() noexcept{
				return const_pow(2.f, -24); 
				//uint32_t native = 0x33800000;
				//return *(float *)&native;
			}
			template<> static constexpr double getEpsilon<double>() noexcept{
				return const_pow(2.0, -53); 
				//uint64_t native = 0x3ca0000000000000;
				//return *(double *)&native;
			}

			template<class FT> static constexpr FT getSplitter() noexcept{
				static_assert(false, "ODER::ExactArthmeticer support IEEE 754-1985 floating point only");
				return 0;
			}
			template<> static constexpr float getSplitter<float>() noexcept{
				return const_pow(2.f, 12) + 1.f;
				//return 4097.f;
			}
			template<> static constexpr double getSplitter<double>() noexcept{
				return const_pow(2.0, 27) + 1.0;
				//return 134217729.0; 
			}

			public:
				static constexpr FT epsilon = getEpsilon<FT>();
				static constexpr FT splitter = getSplitter<FT>();
	};

}

#if defined(_MSC_VER)
#pragma float_control(pop)
#pragma fp_contract(on)
#endif

#endif