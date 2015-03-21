#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_CORE_ARITHMETRIC_H
#define ODER_CORE_ARITHMETRIC_H

#include "oder.h"
#include <cstdint>

#define ABOSOLUTE_GREATER(x, y) ((x)>(y)) == ((x)>(-y)) //returns ture if |x| > |y|

namespace ODER{
	template<class FT> class ExactArthmeticer{
	public:
		inline FT Estimate(FT *e, int n) const{
			FT ret = e[0];
			for (int i = 1; i < n; i++){
				ret += e[i];
			}
			return ret;
		}

		inline void fastTwoSum(FT a, FT b, FT &esti, FT& err) const {
			//|a| >= |b| must be confirmed
			esti = a + b;
			err = b - (esti - a);
		}

		inline FT fastTwoSumError(FT a, FT b, FT esti) const{
			//|a| >= |b| must be confirmed
			return b - (esti - a);
		}

		inline void twoSum(FT a, FT b, FT &esti, FT& err) const {
			esti = a + b;
			FT bvirt = esti - a;
			FT avirt = esti - bvirt;
			err = (b - bvirt) + (a - avirt);
		}

		inline FT twoSumError(FT a, FT b, FT esti) const{
			FT bvirt = esti - a;
			FT avirt = esti - bvirt;
			return (b - bvirt) + (a - avirt);
		}

		inline void twoDiff(FT a, FT b, FT &esti, FT& err) const{
			esti = a - b;
			FT bvirt = a - esti;
			FT avirt = esti + bvirt;
			err = (bvirt - b) + (a - avirt);
		}

		inline FT twoDiffError(FT a, FT b, FT esti) const{
			FT bvirt = a - esti;
			FT avirt = esti + bvirt;
			return (bvirt - b) + (a - avirt);
		}

		inline void Split(FT a, FT &high, FT& low) const {
			FT c = splitter*a;
			FT big = c - a;
			high = c - big;
			low = a - high;
		}

		inline void twoProduct(FT a, FT b, FT& esti, FT& err) const {
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

		inline void twoProductPreSplit(FT a, FT b, FT bHigh, FT bLow, FT& esti, FT& err) const{
			esti = a * b;
			FT aHigh, aLow;
			Split(a, aHigh, aLow);
			FT e0 = esti - (aHigh*bHigh);
			FT e1 = e0 - (aLow*bHigh);
			FT e2 = e1 - (aHigh*bLow);
			err = (aLow*bLow) - e2;
		}

		inline void twoOneDiff(FT a1, FT a0, FT b, FT& x2, FT& x1, FT& x0) const{
			FT iq;
			//a1 + a0 - b
			twoDiff(a0, b, iq, x0);
			twoSum(a1, iq, x2, x1);
		}

		inline void twoTwoDiff(FT a1, FT a0, FT b1, FT b0, FT& x3, FT& x2, FT& x1, FT& x0) const{
			FT i, j;
			twoOneDiff(a1, a0, b0, i, j, x0);
			twoOneDiff(i, j, b1, x3, x2, x1);
		}

		inline void twoOneProduct(FT a1, FT a0, FT b, FT& x3, FT& x2, FT& x1, FT& x0) const{
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
				if (hh != 0.f)
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
					if (hh != 0.f)
						h[hIndex++] = hh;
				}
			}
			while (eIndex < em){
				twoSum(enext, q, qnext, hh);
				enext = e[++eIndex];
				q = qnext;
				if (hh != 0.f)
					h[hIndex++] = hh;
			}
			while (fIndex < fn){
				twoSum(fnext, q, qnext, hh);
				fnext = f[++fIndex];
				q = qnext;
				if (hh != 0.f)
					h[hIndex++] = hh;
			}
			if ((q != 0.f) || (hIndex == 0)){
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

		static const FT epsilon;
		static const FT splitter;

		private:
			template<class FT> static FT getEpsilon(){
				static_assert(false, "ODER::ExactArthmeticer support IEEE 754-1985 floating point only");
				return 0;
			}
			template<> static float getEpsilon<float>(){
				//return pow(2.f, -24); 
				uint32_t native = 0x33800000;
				return *(float *)&native;
			}
			template<> static double getEpsilon<double>(){ 
				//return pow(2.0, -53); 
				uint64_t native = 0x3ca0000000000000;
				return *(double *)&native;
			}

			template<class FT> static constexpr FT getSplitter(){
				static_assert(false, "ODER::ExactArthmeticer support IEEE 754-1985 floating point only");
				return 0;
			}
			template<> static constexpr float getSplitter<float>(){
				//return pow(2.f, 12) + 1.f; 
				return 4097.f;
			}
			template<> static constexpr double getSplitter<double>(){
				//return pow(2.0, 27) + 1.0;
				return 134217729.0; 
			}
	};

	template<class FT> const FT ExactArthmeticer<FT>::epsilon = ExactArthmeticer<FT>::getEpsilon<FT>();
	template<class FT> const FT ExactArthmeticer<FT>::splitter = ExactArthmeticer<FT>::getSplitter<FT>();

}

#endif