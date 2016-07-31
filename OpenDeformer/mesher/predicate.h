/*****************************************************************************/
/*                                                                           */
/*    This file contains implementation of algorithms for exact addition     */
/*    and multiplication of floating-point numbers, and predicates for       */
/*    robustly performing the orientation and incircle tests used in         */
/*    computational geometry.  The algorithms and underlying theory are      */
/*    described in Jonathan Richard Shewchuk. "Adaptive Precision FTing-     */
/*    Point Arithmetic and Fast Robust Geometric Predicates."  Technical     */
/*    Report CMU-CS-96-140, School of Computer Science, Carnegie Mellon      */
/*    University, Pittsburgh, Pennsylvania, May 1996.  (Submitted to         */
/*    Discrete & Computational Geometry.)                                    */
/*                                                                           */
/*****************************************************************************/

#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_MESHER_PREDICATE_H
#define ODER_MESHER_PREDICATE_H

#include "oder.h"
#include "arithmetic.h"

#if defined(_MSC_VER)
#pragma float_control(precise, on, push)
#pragma fp_contract(off)
#endif

#include "latool.h"
#include "geometer.h"

namespace ODER{
template<class FT> class Predicator{
public:
	static_assert(std::is_same<FT, float>::value || std::is_same<FT, double>::value, "ODER::Predicator support IEEE 754-1985 floating point only");

	FT orient2d(FT ax, FT ay, FT bx, FT by, FT cx, FT cy) const;

	FT orient2d(const VectorBase<FT>& a, const VectorBase<FT>& b, const VectorBase<FT>& c, const VectorBase<FT>& above) const;

	FT orient3d(const VectorBase<FT>& a, const VectorBase<FT>& b, const VectorBase<FT>& c, const VectorBase<FT>& d) const;

	inline FT orientCoplane(const VectorBase<FT>& a, const VectorBase<FT>& b, const VectorBase<FT>& c) const;

	FT inCircle(const VectorBase<FT>& a, const VectorBase<FT>& b, const VectorBase<FT>& c, const VectorBase<FT>& d, const VectorBase<FT>& n) const;

	FT inSphere(const VectorBase<FT>& a, const VectorBase<FT>& b, const VectorBase<FT>& c, const VectorBase<FT>& d, const VectorBase<FT>& e) const;

	FT inOrthoCircle(const VectorBase<FT> &a, FT aWeight, const VectorBase<FT> &b, FT bWeight,
		                const VectorBase<FT> &c, FT cWeight, const VectorBase<FT> &d, FT dWeight, const VectorBase<FT> &above) const;

	FT inOrthoCirclePerturbed(const VectorBase<FT> &a, FT aWeight, const VectorBase<FT> &b, FT bWeight,
		                        const VectorBase<FT> &c, FT cWeight, const VectorBase<FT> &d, FT dWeight, const VectorBase<FT> &above) const;

	FT inOrthoSphere(const VectorBase<FT> &a, FT aWeight, const VectorBase<FT> &b, FT bWeight,
		                const VectorBase<FT> &c, FT cWeight, const VectorBase<FT> &d, FT dWeight, const VectorBase<FT> &e, FT eWeight) const;

	FT inOrthoSpherePerturbed(const VectorBase<FT> &a, FT aWeight, const VectorBase<FT> &b, FT bWeight,
		                         const VectorBase<FT> &c, FT cWeight, const VectorBase<FT> &d, FT dWeight, const VectorBase<FT> &e, FT eWeight) const;

	FT inOrthoSphereExact(const VectorBase<FT> &a, FT aWeight, const VectorBase<FT> &b, FT bWeight,
		const VectorBase<FT> &c, FT cWeight, const VectorBase<FT> &d, FT dWeight, const VectorBase<FT> &e, FT eWeight) const;

	bool fastCoPlane(const VectorBase<FT> &a, const VectorBase<FT> &b, const VectorBase<FT>&c, const VectorBase<FT> &d) const{
		VectorBase<FT> n = (b - a) % (c - a);
		return fabs(n*(d - a)) <= arthemetricer.epsilon*std::max(fabs(n.x), std::max(fabs(n.y), fabs(n.z)));
	}

	bool fastCoLine(const VectorBase<FT> &a, const VectorBase<FT> &b, const VectorBase<FT>&c) const{
		VectorBase<FT> n = (b - a) % (c - a);
		return n.length2() <= arthemetricer.epsilon*arthemetricer.epsilon*std::max(fabs(n.x), std::max(fabs(n.y), fabs(n.z)));
	}

	bool inHalfSpace3D(const VectorBase<FT> &u, const VectorBase<FT> &a, const VectorBase<FT>& b, const VectorBase<FT> &c) const;

	bool inHalfSpace2D(const VectorBase<FT>& u, const VectorBase<FT>& a, const VectorBase<FT>& b, const VectorBase<FT>& above) const;

	bool inOrthoHalfSpace3D(const VectorBase<FT> &u, FT uWeight, const VectorBase<FT> &a, FT aWeight, const VectorBase<FT>& b, FT bWeight, const VectorBase<FT> &c, FT cWeight) const;

	bool Intersection(const VectorBase<FT>& p, const VectorBase<FT>& q, const VectorBase<FT>& r,
		const VectorBase<FT>& a, const VectorBase<FT>& b, const VectorBase<FT>& c) const;

	bool Intersection(const VectorBase<FT>& a, const VectorBase<FT>& b, const VectorBase<FT>& c, 
		const VectorBase<FT>& p, const VectorBase<FT>& q) const;

private:

	enum Plane { Plane_XY, Plane_YZ, Plane_XZ, Plane_Arbitary };

	FT orient2dAdaptive(FT ax, FT ay, FT bx, FT by, FT cx, FT cy, FT norm) const;

	FT orient3dAdaptive(const VectorBase<FT>& a, const VectorBase<FT>& b, const VectorBase<FT>& c, const VectorBase<FT>& d, FT norm) const;

	inline FT orientCoplane(const VectorBase<FT>& a, const VectorBase<FT>& b, const VectorBase<FT>& c, Plane hint) const;

	FT inOrthoSphereAdaptive(const VectorBase<FT> &a, FT aWeight, const VectorBase<FT> &b, FT bWeight,
		                   const VectorBase<FT> &c, FT cWeight, const VectorBase<FT> &d, FT dWeight, const VectorBase<FT> &e, FT eWeight, FT norm) const;

	FT inSphereAdaptive(const VectorBase<FT>& a, const VectorBase<FT>& b, const VectorBase<FT>& c, const VectorBase<FT>& d, const VectorBase<FT>& e, FT norm) const;

	FT inSphereExact(const VectorBase<FT>& a, const VectorBase<FT>& b, const VectorBase<FT>& c, const VectorBase<FT>& d, const VectorBase<FT>& e) const;

	bool inSegmentRange(const VectorBase<FT>& u, const VectorBase<FT>& a, const VectorBase<FT>& b) const;

	FT inSegmentRangeHalfAdaptive(const VectorBase<FT>& u, const VectorBase<FT>& a, const VectorBase<FT>& b, 
		const VectorBase<FT>& ab, const VectorBase<FT>& au, FT norm) const;

	bool intersectionTestEdge(const VectorBase<FT>& p, const VectorBase<FT>& q, const VectorBase<FT>& r,
		const VectorBase<FT>& a, const VectorBase<FT>& b, Plane hint) const;

	bool intersectionTestVertex(const VectorBase<FT>& p, const VectorBase<FT>& q, const VectorBase<FT>& r,
		const VectorBase<FT>& a, const VectorBase<FT>& b, const VectorBase<FT>& c, Plane hint) const;

	bool intersectionTriSegCoplane(const VectorBase<FT>& a, const VectorBase<FT>& b, const VectorBase<FT>& c,
		const VectorBase<FT>& p, const VectorBase<FT>& q, Plane hint) const;

	static const ExactArthmeticer<FT> arthemetricer;
	static const FT epsilon;
	static const FT commonBound;
	static const FT o2dErrorBoundA;
	static const FT o2dErrorBoundB;
	static const FT o2dErrorBoundC;
	static const FT o3dErrorBoundA;
	static const FT o3dErrorBoundB;
	static const FT o3dErrorBoundC;
	static const FT inSpeErrorBoundA;
	static const FT inSpeErrorBoundB;
	static const FT inSpeErrorBoundC;
	static const FT inSegRangeErrorBoundA;
	static const FT inSegRangeErrorBoundB;
	static const FT inSegRangeErrorBoundC;
};

template<class FT> const ExactArthmeticer<FT> Predicator<FT>::arthemetricer;
template<class FT> const FT Predicator<FT>::epsilon = FT(1e-8);
template<class FT> const FT Predicator<FT>::commonBound = (FT(3.0) + FT(8.0) * ExactArthmeticer<FT>::epsilon) * ExactArthmeticer<FT>::epsilon;
template<class FT> const FT Predicator<FT>::o2dErrorBoundA = (FT(3.0) + FT(16.0) * ExactArthmeticer<FT>::epsilon) * ExactArthmeticer<FT>::epsilon;
template<class FT> const FT Predicator<FT>::o2dErrorBoundB = (FT(2.0) + FT(12.0) * ExactArthmeticer<FT>::epsilon) * ExactArthmeticer<FT>::epsilon;
template<class FT> const FT Predicator<FT>::o2dErrorBoundC = (FT(9.0) + FT(64.0) * ExactArthmeticer<FT>::epsilon) * ExactArthmeticer<FT>::epsilon * ExactArthmeticer<FT>::epsilon;
template<class FT> const FT Predicator<FT>::o3dErrorBoundA = (FT(7.0) + FT(56.0) * ExactArthmeticer<FT>::epsilon) * ExactArthmeticer<FT>::epsilon;
template<class FT> const FT Predicator<FT>::o3dErrorBoundB = (FT(3.0) + FT(28.0) * ExactArthmeticer<FT>::epsilon) * ExactArthmeticer<FT>::epsilon;
template<class FT> const FT Predicator<FT>::o3dErrorBoundC = (FT(26.0) + FT(288.0) * ExactArthmeticer<FT>::epsilon) * ExactArthmeticer<FT>::epsilon * ExactArthmeticer<FT>::epsilon;
template<class FT> const FT Predicator<FT>::inSpeErrorBoundA = (FT(16.0) + FT(224.0) * ExactArthmeticer<FT>::epsilon) * ExactArthmeticer<FT>::epsilon;
template<class FT> const FT Predicator<FT>::inSpeErrorBoundB = (FT(5.0) + FT(72.0) * ExactArthmeticer<FT>::epsilon) * ExactArthmeticer<FT>::epsilon;
template<class FT> const FT Predicator<FT>::inSpeErrorBoundC = (FT(71.0) + FT(1408.0) * ExactArthmeticer<FT>::epsilon) * ExactArthmeticer<FT>::epsilon * ExactArthmeticer<FT>::epsilon;
template<class FT> const FT Predicator<FT>::inSegRangeErrorBoundA = (FT(5.0) + FT(24.0) * ExactArthmeticer<FT>::epsilon) * ExactArthmeticer<FT>::epsilon;
template<class FT> const FT Predicator<FT>::inSegRangeErrorBoundB = (FT(2.0) + FT(16.0) * ExactArthmeticer<FT>::epsilon) * ExactArthmeticer<FT>::epsilon;
template<class FT> const FT Predicator<FT>::inSegRangeErrorBoundC = (FT(11.0) + FT(80.0) * ExactArthmeticer<FT>::epsilon) * ExactArthmeticer<FT>::epsilon * ExactArthmeticer<FT>::epsilon;

template<class FT> FT Predicator<FT>::orient2d(FT ax, FT ay, FT bx, FT by, FT cx, FT cy) const{
	FT detLeft = (ax - cx)*(by - cy);
	FT detRight = (bx - cx)*(ay - cy);
	FT det = detLeft - detRight;
	FT norm;
	FT errorBound;

	if (detLeft > 0.0){
		if (detRight <= 0.0)
			return det;
		else
			norm = detLeft + detRight;
	}
	else if (detLeft < 0.0){
		if (detRight >= 0.0)
			return det;
		else
			norm = -detLeft - detRight;
	}
	else
		return det;

	errorBound = o2dErrorBoundA*norm;
	if (fabs(det) > errorBound)
		return det;

	return orient2dAdaptive(ax, ay, bx, by, cx, cy, norm);
}

template<class FT> FT Predicator<FT>::orient2dAdaptive(FT ax, FT ay, FT bx, FT by, FT cx, FT cy, FT norm) const{
	FT det, errorBound;
	FT cax = ax - cx;
	FT cay = ay - cy;
	FT cbx = bx - cx;
	FT cby = by - cy;

	FT detLeft, detLeftError, detRight, detRightError;
	arthemetricer.twoProduct(cax, cby, detLeft, detLeftError);
	arthemetricer.twoProduct(cbx, cay, detRight, detRightError);

	FT detEsti[4];
	arthemetricer.twoTwoDiff(detLeft, detLeftError, detRight, detRightError, detEsti[3], detEsti[2], detEsti[1], detEsti[0]);

	det = arthemetricer.Estimate(detEsti, 4);
	errorBound = o2dErrorBoundB*norm;
	if (fabs(det) >= errorBound)
		return det;

	FT caxError = arthemetricer.twoDiffError(ax, cx, cax);
	FT cayError = arthemetricer.twoDiffError(ay, cy, cay);
	FT cbxError = arthemetricer.twoDiffError(bx, cx, cbx);
	FT cbyError = arthemetricer.twoDiffError(by, cy, cby);

	if (caxError == 0.0 && cayError == 0.0 && cbxError == 0.0 && cbyError == 0.0)
		return det;

	errorBound = o2dErrorBoundC*norm + commonBound*fabs(det);
	det += ((cax*cbyError + cby*caxError) - (cbx*cayError + cay*cbxError));
	if (fabs(det) >= errorBound)
		return det;

	FT u[4], c0[8], c1[12], fin[16];
	int c0Len, c1Len, finLength;
	FT s1, s0, t1, t0;

	arthemetricer.twoProduct(caxError, cby, s1, s0);
	arthemetricer.twoProduct(cayError, cbx, t1, t0);
	arthemetricer.twoTwoDiff(s1, s0, t1, t0, u[3], u[2], u[1], u[0]);
	c0Len = arthemetricer.fastExpansionSum(detEsti, u, c0, 4, 4);

	arthemetricer.twoProduct(cbyError, cax, s1, s0);
	arthemetricer.twoProduct(cbxError, cay, t1, t0);
	arthemetricer.twoTwoDiff(s1, s0, t1, t0, u[3], u[2], u[1], u[0]);
	c1Len = arthemetricer.fastExpansionSum(c0, u, c1, c0Len, 4);

	arthemetricer.twoProduct(caxError, cbyError, s1, s0);
	arthemetricer.twoProduct(cayError, cbxError, t1, t0);
	arthemetricer.twoTwoDiff(s1, s0, t1, t0, u[3], u[2], u[1], u[0]);
	finLength = arthemetricer.fastExpansionSum(c1, u, fin, c1Len, 4);

	return fin[finLength - 1];
}

template<class FT> FT Predicator<FT>::orient2d(const VectorBase<FT>& a, const VectorBase<FT>& b, const VectorBase<FT>& c, const VectorBase<FT>& above) const{
	return orient3d(above, a, b, c);
}


template<class FT> FT Predicator<FT>::inCircle(const VectorBase<FT>& a, const VectorBase<FT>& b, const VectorBase<FT>& c, const VectorBase<FT>& d, const VectorBase<FT>& above) const{
	return inSphere(above, a, b, c, d);
}

template<class FT> FT Predicator<FT>::inOrthoCircle(const VectorBase<FT> &a, FT aWeight, const VectorBase<FT> &b, FT bWeight,
	const VectorBase<FT> &c, FT cWeight, const VectorBase<FT> &d, FT dWeight, const VectorBase<FT> &above) const{
	return inOrthoSphere(above, 0.0, a, aWeight, b, bWeight, c, cWeight, d, dWeight);
}

template<class FT> FT Predicator<FT>::orient3d(const VectorBase<FT>& a, const VectorBase<FT>& b, const VectorBase<FT>& c, const VectorBase<FT>& d) const{
	VectorBase<FT> da = a - d;
	VectorBase<FT> db = b - d;
	VectorBase<FT> dc = c - d;

	FT daxdby = da.x*db.y;
	FT daydbx = da.y*db.x;

	FT dbxdcy = db.x*dc.y;
	FT dbydcx = db.y*dc.x;

	FT dcxday = dc.x*da.y;
	FT dcycax = dc.y*da.x;

	FT det = (daxdby - daydbx) * dc.z + (dbxdcy - dbydcx) * da.z + (dcxday - dcycax) * db.z;

	FT norm = (fabs(daxdby) + fabs(daydbx)) * fabs(dc.z)
		+ (fabs(dbxdcy) + fabs(dbydcx)) * fabs(da.z)
		+ (fabs(dcxday) + fabs(dcycax)) * fabs(db.z);

	FT errorBound = o3dErrorBoundA*norm;
	if (fabs(det) >= errorBound)
		return det;

	return orient3dAdaptive(a, b, c, d, norm);
}


template<class FT> FT Predicator<FT>::orient3dAdaptive(const VectorBase<FT>& a, const VectorBase<FT>& b, const VectorBase<FT>& c, const VectorBase<FT>& d, FT norm) const{
	int finLength;
	FT det, errorBound;
	FT fin0[196];
	VectorBase<FT> da = a - d;
	VectorBase<FT> db = b - d;
	VectorBase<FT> dc = c - d;

	//Adaptivly approach
	FT dbxdcy1, dbxdcy0, dbydcx1, dbydcx0;
	arthemetricer.twoProduct(db.x, dc.y, dbxdcy1, dbxdcy0);
	arthemetricer.twoProduct(db.y, dc.x, dbydcx1, dbydcx0);

	FT bc[4], bdet[8];
	arthemetricer.twoTwoDiff(dbxdcy1, dbxdcy0, dbydcx1, dbydcx0, bc[3], bc[2], bc[1], bc[0]);
	int blen = arthemetricer.scaleExpansion(bc, da.z, bdet, 4);

	FT dcxday1, dcxday0, dcydax1, dcydax0;
	arthemetricer.twoProduct(dc.x, da.y, dcxday1, dcxday0);
	arthemetricer.twoProduct(dc.y, da.x, dcydax1, dcydax0);

	FT ca[4], cdet[8];
	arthemetricer.twoTwoDiff(dcxday1, dcxday0, dcydax1, dcydax0, ca[3], ca[2], ca[1], ca[0]);
	int clen = arthemetricer.scaleExpansion(ca, db.z, cdet, 4);

	FT daxdby1, daxdby0, daydbx1, daydbx0;
	arthemetricer.twoProduct(da.x, db.y, daxdby1, daxdby0);
	arthemetricer.twoProduct(da.y, db.x, daydbx1, daydbx0);

	FT ab[4], adet[8];
	arthemetricer.twoTwoDiff(daxdby1, daxdby0, daydbx1, daydbx0, ab[3], ab[2], ab[1], ab[0]);
	int alen = arthemetricer.scaleExpansion(ab, dc.z, adet, 4);

	FT abdet[16];
	int ablen = arthemetricer.fastExpansionSum(adet, bdet, abdet, alen, blen);
	finLength = arthemetricer.fastExpansionSum(abdet, cdet, fin0, ablen, clen);

	det = arthemetricer.Estimate(fin0, finLength);
	errorBound = o3dErrorBoundB*norm;

	if (fabs(det) >= errorBound)
		return det;

	FT daxError, dayError, dazError;
	FT dbxError, dbyError, dbzError;
	FT dcxError, dcyError, dczError;

	daxError = arthemetricer.twoDiffError(a.x, d.x, da.x);
	dayError = arthemetricer.twoDiffError(a.y, d.y, da.y);
	dazError = arthemetricer.twoDiffError(a.z, d.z, da.z);
	dbxError = arthemetricer.twoDiffError(b.x, d.x, db.x);
	dbyError = arthemetricer.twoDiffError(b.y, d.y, db.y);
	dbzError = arthemetricer.twoDiffError(b.z, d.z, db.z);
	dcxError = arthemetricer.twoDiffError(c.x, d.x, dc.x);
	dcyError = arthemetricer.twoDiffError(c.y, d.y, dc.y);
	dczError = arthemetricer.twoDiffError(c.z, d.z, dc.z);

	if (daxError == 0.0 && dayError == 0.0 && dazError == 0.0
		&& dbxError == 0.0 && dbyError == 0.0 && dbzError == 0.0
		&& dcxError == 0.0 && dcyError == 0.0 && dczError == 0.0)
		return det;

	errorBound = commonBound * fabs(det) + o3dErrorBoundC * norm;

	FT adetEstimate = da.z*((db.x*dcyError + dc.y*dbxError) - (db.y*dcxError + dc.x*dbyError))
		+ dazError*(db.x*dc.y - db.y*dc.x);
	FT bdetEstimate = db.z*((dc.x*dayError + da.y*dcxError) - (dc.y*daxError + da.x*dcyError))
		+ dbzError*(dc.x*da.y - dc.y*da.x);
	FT cdetEstimate = dc.z*((da.x*dbyError + db.y*daxError) - (da.y*dbxError + db.x*dayError))
		+ dczError*(da.x*db.y - da.y*db.x);

	det += (adetEstimate + bdetEstimate + cdetEstimate);

	if (fabs(det) >= errorBound)
		return det;

	//Exact Arithmetic
	FT fin1[196];
	FT *finNow = fin0;
	FT *finOther = fin1;

	bool daxErrorZeroTest = (daxError == 0.0);
	bool dayErrorZeroTest = (dayError == 0.0);
	bool dazErrorZeroTest = (dazError == 0.0);
	bool dbxErrorZeroTest = (dbxError == 0.0);
	bool dbyErrorZeroTest = (dbyError == 0.0);
	bool dbzErrorZeroTest = (dbzError == 0.0);
	bool dcxErrorZeroTest = (dcxError == 0.0);
	bool dcyErrorZeroTest = (dcyError == 0.0);
	bool dczErrorZeroTest = (dczError == 0.0);

	FT a_tb[4], a_tc[4], b_ta[4], b_tc[4], c_ta[4], c_tb[4];
	int a_tbLen, a_tcLen, b_taLen, b_tcLen, c_taLen, c_tbLen;
	if (daxErrorZeroTest){
		if (dayErrorZeroTest){
			a_tb[0] = 0.0;
			a_tbLen = 1;
			a_tc[0] = 0.0;
			a_tcLen = 1;
		}
		else{
			arthemetricer.twoProduct(-dayError, db.x, a_tb[1], a_tb[0]);
			a_tbLen = 2;
			arthemetricer.twoProduct(dc.x, dayError, a_tc[1], a_tc[0]);
			a_tcLen = 2;
		}
	}
	else if (dayErrorZeroTest){
		arthemetricer.twoProduct(daxError, db.y, a_tb[1], a_tb[0]);
		a_tbLen = 2;
		arthemetricer.twoProduct(dc.y, -daxError, a_tc[1], a_tc[0]);
		a_tcLen = 2;
	}
	else{
		FT dax_tdby1, dax_tdby0, day_tdbx1, day_tdbx0;
		arthemetricer.twoProduct(daxError, db.y, dax_tdby1, dax_tdby0);
		arthemetricer.twoProduct(dayError, db.x, day_tdbx1, day_tdbx0);
		arthemetricer.twoTwoDiff(dax_tdby1, dax_tdby0, day_tdbx1, day_tdbx0,
			a_tb[3], a_tb[2], a_tb[1], a_tb[0]);
		a_tbLen = 4;

		FT dax_tdcy1, dax_tdcy0, day_tdcx1, day_tdcx0;
		arthemetricer.twoProduct(dayError, dc.x, day_tdcx1, day_tdcx0);
		arthemetricer.twoProduct(daxError, dc.y, dax_tdcy1, dax_tdcy0);
		arthemetricer.twoTwoDiff(day_tdcx1, day_tdcx0, dax_tdcy1, dax_tdcy0,
			a_tc[3], a_tc[2], a_tc[1], a_tc[0]);
		a_tcLen = 4;
	}

	if (dbxErrorZeroTest){
		if (dbyErrorZeroTest){
			b_ta[0] = 0.0;
			b_taLen = 1;
			b_tc[0] = 0.0;
			b_tcLen = 1;
		}
		else{
			arthemetricer.twoProduct(-dbyError, dc.x, b_tc[1], b_tc[0]);
			b_tcLen = 2;
			arthemetricer.twoProduct(da.x, dbyError, b_ta[1], b_ta[0]);
			b_taLen = 2;
		}
	}
	else if (dbyErrorZeroTest){
		arthemetricer.twoProduct(dbxError, dc.y, b_tc[1], b_tc[0]);
		b_tcLen = 2;
		arthemetricer.twoProduct(da.y, -dbxError, b_ta[1], b_ta[0]);
		b_taLen = 2;
	}
	else{
		FT dbx_tdcy1, dbx_tdcy0, dby_tdcx1, dby_tdcx0;
		arthemetricer.twoProduct(dbxError, dc.y, dbx_tdcy1, dbx_tdcy0);
		arthemetricer.twoProduct(dbyError, dc.x, dby_tdcx1, dby_tdcx0);
		arthemetricer.twoTwoDiff(dbx_tdcy1, dbx_tdcy0, dby_tdcx1, dby_tdcx0,
			b_tc[3], b_tc[2], b_tc[1], b_tc[0]);
		b_tcLen = 4;

		FT dbx_tday1, dbx_tday0, dby_tdax1, dby_tdax0;
		arthemetricer.twoProduct(dbyError, da.x, dby_tdax1, dby_tdax0);
		arthemetricer.twoProduct(dbxError, da.y, dbx_tday1, dbx_tday0);
		arthemetricer.twoTwoDiff(dby_tdax1, dby_tdax0, dbx_tday1, dbx_tday0,
			b_ta[3], b_ta[2], b_ta[1], b_ta[0]);
		b_taLen = 4;
	}

	if (dcxErrorZeroTest){
		if (dcyErrorZeroTest){
			c_ta[0] = 0.0;
			c_taLen = 1;
			c_tb[0] = 0.0;
			c_tbLen = 1;
		}
		else{
			arthemetricer.twoProduct(-dcyError, da.x, c_ta[1], c_ta[0]);
			c_taLen = 2;
			arthemetricer.twoProduct(db.x, dcyError, c_tb[1], c_tb[0]);
			c_tbLen = 2;
		}
	}
	else if (dcyErrorZeroTest){
		arthemetricer.twoProduct(dcxError, da.y, c_ta[1], c_ta[0]);
		c_taLen = 2;
		arthemetricer.twoProduct(db.y, -dcxError, c_tb[1], c_tb[0]);
		c_tbLen = 2;
	}
	else{
		FT dcx_tday1, dcx_tday0, dcy_tdax1, dcy_tdax0;
		arthemetricer.twoProduct(dcxError, da.y, dcx_tday1, dcx_tday0);
		arthemetricer.twoProduct(dcyError, da.x, dcy_tdax1, dcy_tdax0);
		arthemetricer.twoTwoDiff(dcx_tday1, dcx_tday0, dcy_tdax1, dcy_tdax0,
			c_ta[3], c_ta[2], c_ta[1], c_ta[0]);
		c_taLen = 4;

		FT dcx_tdby1, dcx_tdby0, dcy_tdbx1, dcy_tdbx0;
		arthemetricer.twoProduct(dcyError, db.x, dcy_tdbx1, dcy_tdbx0);
		arthemetricer.twoProduct(dcxError, db.y, dcx_tdby1, dcx_tdby0);
		arthemetricer.twoTwoDiff(dcy_tdbx1, dcy_tdbx0, dcx_tdby1, dcx_tdby0,
			c_tb[3], c_tb[2], c_tb[1], c_tb[0]);
		c_tbLen = 4;
	}

	FT bct[8], w[16];
	int bctLen = arthemetricer.fastExpansionSum(b_tc, c_tb, bct, b_tcLen, c_tbLen);
	int wLen = arthemetricer.scaleExpansion(bct, da.z, w, bctLen);
	finLength = arthemetricer.fastExpansionSum(finNow, w, finOther, finLength, wLen);
	std::swap(finNow, finOther);

	FT cat[8];
	int catLen = arthemetricer.fastExpansionSum(c_ta, a_tc, cat, c_taLen, a_tcLen);
	wLen = arthemetricer.scaleExpansion(cat, db.z, w, catLen);
	finLength = arthemetricer.fastExpansionSum(finNow, w, finOther, finLength, wLen);
	std::swap(finNow, finOther);

	FT abt[8];
	int abtLen = arthemetricer.fastExpansionSum(a_tb, b_ta, abt, a_tbLen, b_taLen);
	wLen = arthemetricer.scaleExpansion(abt, dc.z, w, abtLen);
	finLength = arthemetricer.fastExpansionSum(finNow, w, finOther, finLength, wLen);
	std::swap(finNow, finOther);

	FT v[8];
	int vLen;
	if (!dazErrorZeroTest){
		vLen = arthemetricer.scaleExpansion(bc, dazError, v, 4);
		finLength = arthemetricer.fastExpansionSum(finNow, v, finOther, finLength, vLen);
		std::swap(finNow, finOther);
	}
	if (!dbzErrorZeroTest){
		vLen = arthemetricer.scaleExpansion(ca, dbzError, v, 4);
		finLength = arthemetricer.fastExpansionSum(finNow, v, finOther, finLength, vLen);
		std::swap(finNow, finOther);
	}
	if (!dczErrorZeroTest){
		vLen = arthemetricer.scaleExpansion(ab, dczError, v, 4);
		finLength = arthemetricer.fastExpansionSum(finNow, v, finOther, finLength, vLen);
		std::swap(finNow, finOther);
	}

	FT u[4];
	if (!daxErrorZeroTest){
		if (!dbyErrorZeroTest){
			FT daxe_dbye1, daxe_dbye0;
			arthemetricer.twoProduct(daxError, dbyError, daxe_dbye1, daxe_dbye0);
			arthemetricer.twoOneProduct(daxe_dbye1, daxe_dbye0, dc.z, u[3], u[2], u[1], u[0]);
			finLength = arthemetricer.fastExpansionSum(finNow, u, finOther, finLength, 4);
			std::swap(finNow, finOther);
			if (!dczErrorZeroTest){
				arthemetricer.twoOneProduct(daxe_dbye1, daxe_dbye0, dczError, u[3], u[2], u[1], u[0]);
				finLength = arthemetricer.fastExpansionSum(finNow, u, finOther, finLength, 4);
				std::swap(finNow, finOther);
			}
		}
		if (!dcyErrorZeroTest){
			FT daxe_dcye1, dcxe_dbye0;
			arthemetricer.twoProduct(-daxError, dbyError, daxe_dcye1, dcxe_dbye0);
			arthemetricer.twoOneProduct(daxe_dcye1, dcxe_dbye0, db.z, u[3], u[2], u[1], u[0]);
			finLength = arthemetricer.fastExpansionSum(finNow, u, finOther, finLength, 4);
			std::swap(finNow, finOther);
			if (!dbzErrorZeroTest){
				arthemetricer.twoOneProduct(daxe_dcye1, dcxe_dbye0, dbzError, u[3], u[2], u[1], u[0]);
				finLength = arthemetricer.fastExpansionSum(finNow, u, finOther, finLength, 4);
				std::swap(finNow, finOther);
			}
		}
	}

	if (!dbxErrorZeroTest){
		if (!dcyErrorZeroTest){
			FT dbxe_dcye1, dbxe_dcye0;
			arthemetricer.twoProduct(dbxError, dcyError, dbxe_dcye1, dbxe_dcye0);
			arthemetricer.twoOneProduct(dbxe_dcye1, dbxe_dcye0, da.z, u[3], u[2], u[1], u[0]);
			finLength = arthemetricer.fastExpansionSum(finNow, u, finOther, finLength, 4);
			std::swap(finNow, finOther);
			if (!dazErrorZeroTest){
				arthemetricer.twoOneProduct(dbxe_dcye1, dbxe_dcye0, dazError, u[3], u[2], u[1], u[0]);
				finLength = arthemetricer.fastExpansionSum(finNow, u, finOther, finLength, 4);
				std::swap(finNow, finOther);
			}
		}
		if (!dayErrorZeroTest){
			FT dbxe_daye1, dbxe_dcye0;
			arthemetricer.twoProduct(-dbxError, dayError, dbxe_daye1, dbxe_dcye0);
			arthemetricer.twoOneProduct(dbxe_daye1, dbxe_dcye0, dc.z, u[3], u[2], u[1], u[0]);
			finLength = arthemetricer.fastExpansionSum(finNow, u, finOther, finLength, 4);
			std::swap(finNow, finOther);
			if (!dczErrorZeroTest){
				arthemetricer.twoOneProduct(dbxe_daye1, dbxe_dcye0, dczError, u[3], u[2], u[1], u[0]);
				finLength = arthemetricer.fastExpansionSum(finNow, u, finOther, finLength, 4);
				std::swap(finNow, finOther);
			}
		}
	}

	if (!dcxErrorZeroTest){
		if (!dayErrorZeroTest){
			FT dcxe_daye1, dcxe_daye0;
			arthemetricer.twoProduct(dcxError, dayError, dcxe_daye1, dcxe_daye0);
			arthemetricer.twoOneProduct(dcxe_daye1, dcxe_daye0, db.z, u[3], u[2], u[1], u[0]);
			finLength = arthemetricer.fastExpansionSum(finNow, u, finOther, finLength, 4);
			std::swap(finNow, finOther);
			if (!dbzErrorZeroTest){
				arthemetricer.twoOneProduct(dcxe_daye1, dcxe_daye0, dbzError, u[3], u[2], u[1], u[0]);
				finLength = arthemetricer.fastExpansionSum(finNow, u, finOther, finLength, 4);
				std::swap(finNow, finOther);
			}
		}
		if (!dbyErrorZeroTest){
			FT dcxe_dcye1, dcxe_dbye0;
			arthemetricer.twoProduct(-dcxError, dbyError, dcxe_dcye1, dcxe_dbye0);
			arthemetricer.twoOneProduct(dcxe_dcye1, dcxe_dbye0, dc.z, u[3], u[2], u[1], u[0]);
			finLength = arthemetricer.fastExpansionSum(finNow, u, finOther, finLength, 4);
			std::swap(finNow, finOther);
			if (!dazErrorZeroTest){
				arthemetricer.twoOneProduct(dcxe_dcye1, dcxe_dbye0, dczError, u[3], u[2], u[1], u[0]);
				finLength = arthemetricer.fastExpansionSum(finNow, u, finOther, finLength, 4);
				std::swap(finNow, finOther);
			}
		}
	}

	if (!dazErrorZeroTest){
		wLen = arthemetricer.scaleExpansion(bct, dazError, w, bctLen);
		finLength = arthemetricer.fastExpansionSum(finNow, w, finOther, finLength, wLen);
		std::swap(finNow, finOther);
	}
	if (!dbzErrorZeroTest){
		wLen = arthemetricer.scaleExpansion(cat, dbzError, w, catLen);
		finLength = arthemetricer.fastExpansionSum(finNow, w, finOther, finLength, wLen);
		std::swap(finNow, finOther);
	}
	if (!dczErrorZeroTest){
		wLen = arthemetricer.scaleExpansion(abt, dazError, w, abtLen);
		finLength = arthemetricer.fastExpansionSum(finNow, w, finOther, finLength, wLen);
		std::swap(finNow, finOther);
	}
	return finNow[finLength - 1];
}

template<class FT> FT Predicator<FT>::orientCoplane(const VectorBase<FT>& a, const VectorBase<FT>& b, const VectorBase<FT>& c) const {
	FT xy = orient2d(a.x, a.y, b.x, b.y, c.x, c.y);
	if (xy != 0) return xy;
	FT xz = orient2d(a.x, a.z, b.x, b.z, c.x, c.z);
	if (xz != 0) return xz;
	return orient2d(a.y, a.z, b.y, b.z, c.y, c.z);
}

template<class FT> FT Predicator<FT>::orientCoplane(const VectorBase<FT>& a, const VectorBase<FT>& b, const VectorBase<FT>& c, Plane hint) const {
	switch (hint) {
	case Plane::Plane_XY:
		return orient2d(a.x, a.y, b.x, b.y, c.x, c.y);
	case Plane::Plane_XZ:
		return orient2d(a.x, a.z, b.x, b.z, c.x, c.z);
	case Plane::Plane_YZ:
		return orient2d(a.y, a.z, b.y, b.z, c.y, c.z);
	default:
		return orientCoplane(a, b, c);
	}

	return orientCoplane(a, b, c); // never here
}

template<class FT> FT Predicator<FT>::inSphere(const VectorBase<FT>& a, const VectorBase<FT>& b, const VectorBase<FT>& c, const VectorBase<FT>& d, const VectorBase<FT>& e) const{
	VectorBase<FT> ea = a - e;
	VectorBase<FT> eb = b - e;
	VectorBase<FT> ec = c - e;
	VectorBase<FT> ed = d - e;

	FT eaxeby = ea.x*eb.y;
	FT eayebx = ea.y*eb.x;
	FT ab = eaxeby - eayebx;

	FT ebxecy = eb.x*ec.y;
	FT ebyecx = eb.y*ec.x;
	FT bc = ebxecy - ebyecx;

	FT ecxedy = ec.x*ed.y;
	FT ecyedx = ec.y*ed.x;
	FT cd = ecxedy - ecyedx;

	FT edxeay = ed.x*ea.y;
	FT edyeax = ed.y*ea.x;
	FT da = edxeay - edyeax;

	FT eaxecy = ea.x*ec.y;
	FT eayecx = ea.y*ec.x;
	FT ac = eaxecy - eayecx;

	FT ebxedy = eb.x*ed.y;
	FT ebyedx = eb.y*ed.x;
	FT bd = ebxedy - ebyedx;

	FT abc = ea.z*bc - eb.z*ac + ec.z*ab;
	FT bcd = eb.z*cd - ec.z*bd + ed.z*bc;
	FT cda = ec.z*da + ed.z*ac + ea.z*cd;
	FT dab = ed.z*ab + ea.z*bd + eb.z*da;

	FT aLift = ea.length2();
	FT bLift = eb.length2();
	FT cLift = ec.length2();
	FT dLift = ed.length2();

	FT det = -aLift*bcd + bLift*cda - cLift*dab + dLift*abc;

	FT eaxebyAbsolute = fabs(eaxeby);
	FT eayebxAbsolute = fabs(eayebx);
	FT ebxecyAbsolute = fabs(ebxecy);
	FT ebyecxAbsolute = fabs(ebyecx);
	FT ecxedyAbsolute = fabs(ecxedy);
	FT ecyedxAbsolute = fabs(ecyedx);
	FT edxeayAbsolute = fabs(edxeay);
	FT edyeaxAbsolute = fabs(edyeax);
	FT eaxecyAbsolute = fabs(eaxecy);
	FT eayecxAbsolute = fabs(eayecx);
	FT ebxedyAbsolute = fabs(ebxedy);
	FT ebyedxAbsolute = fabs(ebyedx);
	FT eazAbsolute = fabs(ea.z);
	FT ebzAbsolute = fabs(eb.z);
	FT eczAbsolute = fabs(ec.z);
	FT edzAbsolute = fabs(ed.z);

	abc = eazAbsolute*(ebxecyAbsolute + ebyecxAbsolute) + ebzAbsolute*(ecxedyAbsolute + eayecxAbsolute) + eczAbsolute*(eaxebyAbsolute + eayebxAbsolute);
	bcd = ebzAbsolute*(ecxedyAbsolute + ecyedxAbsolute) + ebzAbsolute*(ebxedyAbsolute + ebyedxAbsolute) + edzAbsolute*(ebxecyAbsolute + ebyecxAbsolute);
	cda = eczAbsolute*(edxeayAbsolute + edyeaxAbsolute) + edzAbsolute*(eaxecyAbsolute + eayecxAbsolute) + eazAbsolute*(ecxedyAbsolute + ecyedxAbsolute);
	dab = edzAbsolute*(eaxebyAbsolute + eayebxAbsolute) + eazAbsolute*(ebxedyAbsolute + ebyedxAbsolute) + ebzAbsolute*(edxeayAbsolute + edyeaxAbsolute);

	FT norm = aLift*bcd + bLift*cda + cLift*dab + dLift*abc;
	FT errorBound = inSpeErrorBoundA*norm;
	if (fabs(det) >= errorBound)
		return det;

	return inSphereAdaptive(a, b, c, d, e, norm);
}

template<class FT> FT Predicator<FT>::inSphereAdaptive(const VectorBase<FT>& a, const VectorBase<FT>& b, const VectorBase<FT>& c, const VectorBase<FT>& d, const VectorBase<FT>& e, FT norm) const{
	FT fin[1152];
	int finLength;
	FT det, errorBound;

	VectorBase<FT> ea = a - e;
	VectorBase<FT> eb = b - e;
	VectorBase<FT> ec = c - e;
	VectorBase<FT> ed = d - e;

	FT eaxeby1, eaxeby0, eayebx1, eayebx0, ab[4];
	arthemetricer.twoProduct(ea.x, eb.y, eaxeby1, eaxeby0);
	arthemetricer.twoProduct(ea.y, eb.x, eayebx1, eayebx0);
	arthemetricer.twoTwoDiff(eaxeby1, eaxeby0, eayebx1, eayebx0, ab[3], ab[2], ab[1], ab[0]);

	FT ebxecy1, ebxecy0, ebyecx1, ebyecx0, bc[4];
	arthemetricer.twoProduct(eb.x, ec.y, ebxecy1, ebxecy0);
	arthemetricer.twoProduct(eb.y, ec.x, ebyecx1, ebyecx0);
	arthemetricer.twoTwoDiff(ebxecy1, ebxecy0, ebyecx1, ebyecx0, bc[3], bc[2], bc[1], bc[0]);

	FT ecxedy1, ecxedy0, ecyedx1, ecyedx0, cd[4];
	arthemetricer.twoProduct(ec.x, ed.y, ecxedy1, ecxedy0);
	arthemetricer.twoProduct(ec.y, ed.x, ecyedx1, ecyedx0);
	arthemetricer.twoTwoDiff(ecxedy1, ecxedy0, ecyedx1, ecyedx0, cd[3], cd[2], cd[1], cd[0]);

	FT edxeay1, edxeay0, edyeax1, edyeax0, da[4];
	arthemetricer.twoProduct(ed.x, ea.y, edxeay1, edxeay0);
	arthemetricer.twoProduct(ed.y, ea.x, edyeax1, edyeax0);
	arthemetricer.twoTwoDiff(edxeay1, edxeay0, edyeax1, edyeax0, da[3], da[2], da[1], da[0]);

	FT eaxecy1, eaxecy0, eayecx1, eayecx0, ac[4];
	arthemetricer.twoProduct(ea.x, ec.y, eaxecy1, eaxecy0);
	arthemetricer.twoProduct(ea.y, ec.x, eayecx1, eayecx0);
	arthemetricer.twoTwoDiff(eaxecy1, eaxecy0, eayecx1, eayecx0, ac[3], ac[2], ac[1], ac[0]);

	FT ebxedy1, ebxedy0, ebyedx1, ebyedx0, bd[4];
	arthemetricer.twoProduct(eb.x, ed.y, ebxedy1, ebxedy0);
	arthemetricer.twoProduct(eb.y, ed.x, ebyedx1, ebyedx0);
	arthemetricer.twoTwoDiff(ebxedy1, ebxedy0, ebyedx1, ebyedx0, bd[3], bd[2], bd[1], bd[0]);

	FT temp8a[8], temp8b[8], temp8c[8], temp16[16], temp24[24], temp48[48];
	int temp8aLen, temp8bLen, temp8cLen, temp16Len, temp24Len, temp48Len;
	FT xdet[96], ydet[96], zdet[96], xydet[192];
	int xdetLen, ydetLen, zdetLen, xydetLen;
	FT adet[288], bdet[288], cdet[288], ddet[288];
	int adetLen, bdetLen, cdetLen, ddetLen;

	//abc = ea.z*bc - eb.z*ac + ec.z*ab;
	temp8aLen = arthemetricer.scaleExpansion(bc, ea.z, temp8a, 4);
	temp8bLen = arthemetricer.scaleExpansion(ac, -eb.z, temp8b, 4);
	temp8cLen = arthemetricer.scaleExpansion(ab, ec.z, temp8c, 4);
	temp16Len = arthemetricer.fastExpansionSum(temp8a, temp8b, temp16, temp8aLen, temp8bLen);
	temp24Len = arthemetricer.fastExpansionSum(temp8c, temp16, temp24, temp8cLen, temp16Len);

	//ed.leng2()*abc
	temp48Len = arthemetricer.scaleExpansion(temp24, ed.x, temp48, temp24Len);
	xdetLen = arthemetricer.scaleExpansion(temp48, ed.x, xdet, temp48Len);
	temp48Len = arthemetricer.scaleExpansion(temp24, ed.y, temp48, temp24Len);
	ydetLen = arthemetricer.scaleExpansion(temp48, ed.y, ydet, temp48Len);
	temp48Len = arthemetricer.scaleExpansion(temp24, ed.z, temp48, temp24Len);
	zdetLen = arthemetricer.scaleExpansion(temp48, ed.z, zdet, temp48Len);
	xydetLen = arthemetricer.fastExpansionSum(xdet, ydet, xydet, xdetLen, ydetLen);
	ddetLen = arthemetricer.fastExpansionSum(xydet, zdet, ddet, xydetLen, zdetLen);

	//bcd = eb.z*cd - ec.z*bd + ed.z*bc
	temp8aLen = arthemetricer.scaleExpansion(cd, eb.z, temp8a, 4);
	temp8bLen = arthemetricer.scaleExpansion(bd, -ec.z, temp8b, 4);
	temp8cLen = arthemetricer.scaleExpansion(bc, ed.z, temp8c, 4);
	temp16Len = arthemetricer.fastExpansionSum(temp8a, temp8b, temp16, temp8aLen, temp8bLen);
	temp24Len = arthemetricer.fastExpansionSum(temp8c, temp16, temp24, temp8cLen, temp16Len);

	//-ea.length2()*bcd
	temp48Len = arthemetricer.scaleExpansion(temp24, ea.x, temp48, temp24Len);
	xdetLen = arthemetricer.scaleExpansion(temp48, -ea.x, xdet, temp48Len);
	temp48Len = arthemetricer.scaleExpansion(temp24, ea.y, temp48, temp24Len);
	ydetLen = arthemetricer.scaleExpansion(temp48, -ea.y, ydet, temp48Len);
	temp48Len = arthemetricer.scaleExpansion(temp24, ea.z, temp48, temp24Len);
	zdetLen = arthemetricer.scaleExpansion(temp48, -ea.z, zdet, temp48Len);
	xydetLen = arthemetricer.fastExpansionSum(xdet, ydet, xydet, xdetLen, ydetLen);
	adetLen = arthemetricer.fastExpansionSum(xydet, zdet, adet, xydetLen, zdetLen);

	//cda = ec.z*da + ed.z*ac + ea.z*cd
	temp8aLen = arthemetricer.scaleExpansion(da, ec.z, temp8a, 4);
	temp8bLen = arthemetricer.scaleExpansion(ac, ed.z, temp8b, 4);
	temp8cLen = arthemetricer.scaleExpansion(cd, ea.z, temp8c, 4);
	temp16Len = arthemetricer.fastExpansionSum(temp8a, temp8b, temp16, temp8aLen, temp8bLen);
	temp24Len = arthemetricer.fastExpansionSum(temp8c, temp16, temp24, temp8cLen, temp16Len);

	//eb.length2()*cda
	temp48Len = arthemetricer.scaleExpansion(temp24, eb.x, temp48, temp24Len);
	xdetLen = arthemetricer.scaleExpansion(temp48, eb.x, xdet, temp48Len);
	temp48Len = arthemetricer.scaleExpansion(temp24, eb.y, temp48, temp24Len);
	ydetLen = arthemetricer.scaleExpansion(temp48, eb.y, ydet, temp48Len);
	temp48Len = arthemetricer.scaleExpansion(temp24, eb.z, temp48, temp24Len);
	zdetLen = arthemetricer.scaleExpansion(temp48, eb.z, zdet, temp48Len);
	xydetLen = arthemetricer.fastExpansionSum(xdet, ydet, xydet, xdetLen, ydetLen);
	bdetLen = arthemetricer.fastExpansionSum(xydet, zdet, bdet, xydetLen, zdetLen);

	//FT dab = ed.z*ab + ea.z*bd + eb.z*da
	temp8aLen = arthemetricer.scaleExpansion(ab, ed.z, temp8a, 4);
	temp8bLen = arthemetricer.scaleExpansion(bd, ea.z, temp8b, 4);
	temp8cLen = arthemetricer.scaleExpansion(da, eb.z, temp8c, 4);
	temp16Len = arthemetricer.fastExpansionSum(temp8a, temp8b, temp16, temp8aLen, temp8bLen);
	temp24Len = arthemetricer.fastExpansionSum(temp8c, temp16, temp24, temp8cLen, temp16Len);

	//-ec.length2()*dab
	temp48Len = arthemetricer.scaleExpansion(temp24, ec.x, temp48, temp24Len);
	xdetLen = arthemetricer.scaleExpansion(temp48, -ec.x, xdet, temp48Len);
	temp48Len = arthemetricer.scaleExpansion(temp24, ec.y, temp48, temp24Len);
	ydetLen = arthemetricer.scaleExpansion(temp48, -ec.y, ydet, temp48Len);
	temp48Len = arthemetricer.scaleExpansion(temp24, ec.z, temp48, temp24Len);
	zdetLen = arthemetricer.scaleExpansion(temp48, -ec.z, zdet, temp48Len);
	xydetLen = arthemetricer.fastExpansionSum(xdet, ydet, xydet, xdetLen, ydetLen);
	cdetLen = arthemetricer.fastExpansionSum(xydet, zdet, cdet, xydetLen, zdetLen);

	FT abdet[576], cddet[576];
	int abdetLen, cddetLen;
	abdetLen = arthemetricer.fastExpansionSum(adet, bdet, abdet, adetLen, bdetLen);
	cddetLen = arthemetricer.fastExpansionSum(cdet, ddet, cddet, cdetLen, ddetLen);
	finLength = arthemetricer.fastExpansionSum(abdet, cddet, fin, abdetLen, cddetLen);

	det = arthemetricer.Estimate(fin, finLength);
	errorBound = inSpeErrorBoundB*norm;
	if (fabs(det) >= errorBound)
		return det;

	FT eaxError = arthemetricer.twoDiffError(a.x, e.x, ea.x);
	FT eayError = arthemetricer.twoDiffError(a.y, e.y, ea.y);
	FT eazError = arthemetricer.twoDiffError(a.z, e.z, ea.z);

	FT ebxError = arthemetricer.twoDiffError(b.x, e.x, eb.x);
	FT ebyError = arthemetricer.twoDiffError(b.y, e.y, eb.y);
	FT ebzError = arthemetricer.twoDiffError(b.z, e.z, eb.z);

	FT ecxError = arthemetricer.twoDiffError(c.x, e.x, ec.x);
	FT ecyError = arthemetricer.twoDiffError(c.y, e.y, ec.y);
	FT eczError = arthemetricer.twoDiffError(c.z, e.z, ec.z);

	FT edxError = arthemetricer.twoDiffError(d.x, e.x, ed.x);
	FT edyError = arthemetricer.twoDiffError(d.y, e.y, ed.y);
	FT edzError = arthemetricer.twoDiffError(d.z, e.z, ed.z);

	if (eaxError == 0.0 && eayError == 0.0 && eazError == 0.0
		&& ebxError == 0.0 && ebyError == 0.0 && ebzError == 0.0
		&& ecxError == 0.0 && ecyError == 0.0 && eczError == 0.0
		&& edxError == 0.0 && edyError == 0.0 && edzError == 0.0)
		return det;

	FT abError = (ea.x*ebyError + eb.y*eaxError) - (ea.y*ebxError + eb.x*eayError);
	FT bcError = (eb.x*ecyError + ec.y*ebxError) - (eb.y*ecxError + ec.x*ebyError);
	FT cdError = (ec.x*edyError + ed.y*ecxError) - (ec.y*edxError + ed.x*ecyError);
	FT daError = (ed.x*eayError + ea.y*edxError) - (ed.y*eaxError + ea.x*edyError);
	FT acError = (ea.x*ecyError + ec.y*eaxError) - (ea.y*ecxError + ec.x*eayError);
	FT bdError = (eb.x*edyError + ed.y*ebxError) - (eb.y*edxError + ed.x*ebyError);

	FT abc = ea.z*bc[3] - eb.z*ac[3] + ec.z*ab[3];
	FT bcd = eb.z*cd[3] - ec.z*bd[3] + ed.z*bc[3];
	FT cda = ec.z*da[3] + ed.z*ac[3] + ea.z*cd[3];
	FT dab = ed.z*ab[3] + ea.z*bd[3] + eb.z*da[3];

	FT abcError = (ea.z*bcError - eb.z*acError + ec.z*abError) + (eazError*bc[3] - ebzError*ac[3] + eczError*ab[3]);
	FT bcdError = (eb.z*cdError - ec.z*bdError + ed.z*bcError) + (ebzError*cd[3] - eczError*bd[3] + edzError*bc[3]);
	FT cdaError = (ec.z*daError + ed.z*acError + ea.z*cdError) + (eczError*da[3] + edzError*ac[3] + eazError*cd[3]);
	FT dabError = (ed.z*abError + ea.z*bdError + eb.z*daError) + (edzError*ab[3] + eazError*bd[3] + ebzError*da[3]);

	FT eaLength2ErrorOne = ea.x*eaxError + ea.y*eayError + ea.z*eazError;
	FT ebLength2ErrorOne = eb.x*ebxError + eb.y*ebyError + eb.z*ebzError;
	FT ecLength2ErrorOne = ec.x*ecxError + ec.y*ecyError + ec.z*eczError;
	FT edLength2ErrorOne = ed.x*edxError + ed.y*edyError + ed.z*edzError;


	det += (-ea.length2()*bcdError + eb.length2()*cdaError - ec.length2()*dabError + ed.length2()*abcError
		+ FT(2)*(-eaLength2ErrorOne*bcd + ebLength2ErrorOne*cda - ecLength2ErrorOne*dab + edLength2ErrorOne*abc));

	errorBound = inSpeErrorBoundC*norm + commonBound*fabs(det);
	if (fabs(det) >= errorBound)
		return det;

	return inSphereExact(a, b, c, d, e);
}

template<class FT> FT Predicator<FT>::inSphereExact(const VectorBase<FT>& a, const VectorBase<FT>& b, const VectorBase<FT>& c, const VectorBase<FT>& d, const VectorBase<FT>& e) const{
	FT axby1, axby0, aybx1, aybx0, ab[4];
	arthemetricer.twoProduct(a.x, b.y, axby1, axby0);
	arthemetricer.twoProduct(a.y, b.x, aybx1, aybx0);
	arthemetricer.twoTwoDiff(axby1, axby0, aybx1, aybx0, ab[3], ab[2], ab[1], ab[0]);

	FT bxcy1, bxcy0, bycx1, bycx0, bc[4];
	arthemetricer.twoProduct(b.x, c.y, bxcy1, bxcy0);
	arthemetricer.twoProduct(b.y, c.x, bycx1, bycx0);
	arthemetricer.twoTwoDiff(bxcy1, bxcy0, bycx1, bycx0, bc[3], bc[2], bc[1], bc[0]);

	FT cxdy1, cxdy0, cydx1, cydx0, cd[4];
	arthemetricer.twoProduct(c.x, d.y, cxdy1, cxdy0);
	arthemetricer.twoProduct(c.y, d.x, cydx1, cydx0);
	arthemetricer.twoTwoDiff(cxdy1, cxdy0, cydx1, cydx0, cd[3], cd[2], cd[1], cd[0]);

	FT dxey1, dxey0, dyex1, dyex0, de[4];
	arthemetricer.twoProduct(d.x, e.y, dxey1, dxey0);
	arthemetricer.twoProduct(d.y, e.x, dyex1, dyex0);
	arthemetricer.twoTwoDiff(dxey1, dxey0, dyex1, dyex0, de[3], de[2], de[1], de[0]);

	FT exay1, exay0, eyax1, eyax0, ea[4];
	arthemetricer.twoProduct(e.x, a.y, exay1, exay0);
	arthemetricer.twoProduct(e.y, a.x, eyax1, eyax0);
	arthemetricer.twoTwoDiff(exay1, exay0, eyax1, eyax0, ea[3], ea[2], ea[1], ea[0]);

	FT axcy1, axcy0, aycx1, aycx0, ac[4];
	arthemetricer.twoProduct(a.x, c.y, axcy1, axcy0);
	arthemetricer.twoProduct(a.y, c.x, aycx1, aycx0);
	arthemetricer.twoTwoDiff(axcy1, axcy0, aycx1, aycx0, ac[3], ac[2], ac[1], ac[0]);

	FT bxdy1, bxdy0, bydx1, bydx0, bd[4];
	arthemetricer.twoProduct(b.x, d.y, bxdy1, bxdy0);
	arthemetricer.twoProduct(b.y, d.x, bydx1, bydx0);
	arthemetricer.twoTwoDiff(bxdy1, bxdy0, bydx1, bydx0, bd[3], bd[2], bd[1], bd[0]);

	FT cxey1, cxey0, cyex1, cyex0, ce[4];
	arthemetricer.twoProduct(c.x, e.y, cxey1, cxey0);
	arthemetricer.twoProduct(c.y, e.x, cyex1, cyex0);
	arthemetricer.twoTwoDiff(cxey1, cxey0, cyex1, cyex0, ce[3], ce[2], ce[1], ce[0]);

	FT dxay1, dxay0, dyax1, dyax0, da[4];
	arthemetricer.twoProduct(d.x, a.y, dxay1, dxay0);
	arthemetricer.twoProduct(d.y, a.x, dyax1, dyax0);
	arthemetricer.twoTwoDiff(dxay1, dxay0, dyax1, dyax0, da[3], da[2], da[1], da[0]);

	FT exby1, exby0, eybx1, eybx0, eb[4];
	arthemetricer.twoProduct(e.x, b.y, exby1, exby0);
	arthemetricer.twoProduct(e.y, b.x, eybx1, eybx0);
	arthemetricer.twoTwoDiff(exby1, exby0, eybx1, eybx0, eb[3], eb[2], eb[1], eb[0]);


	FT abc[24], bcd[24], cde[24], dea[24], eab[24], abd[24], bce[24], cda[24], deb[24], eac[24];
	int abcLen, bcdLen, cdeLen, deaLen, eabLen, abdLen, bceLen, cdaLen, debLen, eacLen;
	FT temp8a[8], temp8b[9], temp16[16];
	int temp8aLen, temp8bLen, temp16Len;

	//abc
	temp8aLen = arthemetricer.scaleExpansion(bc, a.z, temp8a, 4);
	temp8bLen = arthemetricer.scaleExpansion(ac, -b.z, temp8b, 4);
	temp16Len = arthemetricer.fastExpansionSum(temp8a, temp8b, temp16, temp8aLen, temp8bLen);
	temp8aLen = arthemetricer.scaleExpansion(ab, c.z, temp8a, 4);
	abcLen = arthemetricer.fastExpansionSum(temp16, temp8a, abc, temp16Len, temp8aLen);

	//bcd
	temp8aLen = arthemetricer.scaleExpansion(cd, b.z, temp8a, 4);
	temp8bLen = arthemetricer.scaleExpansion(bd, -c.z, temp8b, 4);
	temp16Len = arthemetricer.fastExpansionSum(temp8a, temp8b, temp16, temp8aLen, temp8bLen);
	temp8aLen = arthemetricer.scaleExpansion(bc, d.z, temp8a, 4);
	bcdLen = arthemetricer.fastExpansionSum(temp16, temp8a, bcd, temp16Len, temp8aLen);

	//cde
	temp8aLen = arthemetricer.scaleExpansion(de, c.z, temp8a, 4);
	temp8bLen = arthemetricer.scaleExpansion(ce, -d.z, temp8b, 4);
	temp16Len = arthemetricer.fastExpansionSum(temp8a, temp8b, temp16, temp8aLen, temp8bLen);
	temp8aLen = arthemetricer.scaleExpansion(cd, e.z, temp8a, 4);
	cdeLen = arthemetricer.fastExpansionSum(temp16, temp8a, cde, temp16Len, temp8aLen);

	//dea
	temp8aLen = arthemetricer.scaleExpansion(ea, d.z, temp8a, 4);
	temp8bLen = arthemetricer.scaleExpansion(da, -e.z, temp8b, 4);
	temp16Len = arthemetricer.fastExpansionSum(temp8a, temp8b, temp16, temp8aLen, temp8bLen);
	temp8aLen = arthemetricer.scaleExpansion(de, a.z, temp8a, 4);
	deaLen = arthemetricer.fastExpansionSum(temp16, temp8a, dea, temp16Len, temp8aLen);

	//eab
	temp8aLen = arthemetricer.scaleExpansion(ab, e.z, temp8a, 4);
	temp8bLen = arthemetricer.scaleExpansion(eb, -a.z, temp8b, 4);
	temp16Len = arthemetricer.fastExpansionSum(temp8a, temp8b, temp16, temp8aLen, temp8bLen);
	temp8aLen = arthemetricer.scaleExpansion(ea, b.z, temp8a, 4);
	eabLen = arthemetricer.fastExpansionSum(temp16, temp8a, eab, temp16Len, temp8aLen);

	//abd
	temp8aLen = arthemetricer.scaleExpansion(bd, a.z, temp8a, 4);
	temp8bLen = arthemetricer.scaleExpansion(da, b.z, temp8b, 4);
	temp16Len = arthemetricer.fastExpansionSum(temp8a, temp8b, temp16, temp8aLen, temp8bLen);
	temp8aLen = arthemetricer.scaleExpansion(ab, d.z, temp8a, 4);
	abdLen = arthemetricer.fastExpansionSum(temp16, temp8a, abd, temp16Len, temp8aLen);

	//bce
	temp8aLen = arthemetricer.scaleExpansion(ce, b.z, temp8a, 4);
	temp8bLen = arthemetricer.scaleExpansion(eb, c.z, temp8b, 4);
	temp16Len = arthemetricer.fastExpansionSum(temp8a, temp8b, temp16, temp8aLen, temp8bLen);
	temp8aLen = arthemetricer.scaleExpansion(bc, e.z, temp8a, 4);
	bceLen = arthemetricer.fastExpansionSum(temp16, temp8a, bce, temp16Len, temp8aLen);

	//cda
	temp8aLen = arthemetricer.scaleExpansion(da, c.z, temp8a, 4);
	temp8bLen = arthemetricer.scaleExpansion(ac, d.z, temp8b, 4);
	temp16Len = arthemetricer.fastExpansionSum(temp8a, temp8b, temp16, temp8aLen, temp8bLen);
	temp8aLen = arthemetricer.scaleExpansion(cd, a.z, temp8a, 4);
	cdaLen = arthemetricer.fastExpansionSum(temp16, temp8a, cda, temp16Len, temp8aLen);

	//deb
	temp8aLen = arthemetricer.scaleExpansion(eb, d.z, temp8a, 4);
	temp8bLen = arthemetricer.scaleExpansion(bd, e.z, temp8b, 4);
	temp16Len = arthemetricer.fastExpansionSum(temp8a, temp8b, temp16, temp8aLen, temp8bLen);
	temp8aLen = arthemetricer.scaleExpansion(de, b.z, temp8a, 4);
	debLen = arthemetricer.fastExpansionSum(temp16, temp8a, deb, temp16Len, temp8aLen);

	//eac
	temp8aLen = arthemetricer.scaleExpansion(ac, e.z, temp8a, 4);
	temp8bLen = arthemetricer.scaleExpansion(ce, a.z, temp8b, 4);
	temp16Len = arthemetricer.fastExpansionSum(temp8a, temp8b, temp16, temp8aLen, temp8bLen);
	temp8aLen = arthemetricer.scaleExpansion(ea, c.z, temp8a, 4);
	eacLen = arthemetricer.fastExpansionSum(temp16, temp8a, eac, temp16Len, temp8aLen);

	FT xyzw[96];
	int xyzwLen;
	FT temp48a[48], temp48b[48], temp192[192], temp384x[384], temp384y[384], temp384z[384], temp768xy[768];
	int temp48aLen, temp48bLen, xLen, yLen, zLen, xyLen;

	FT adet[1152];
	int adetLen;
	temp48aLen = arthemetricer.fastExpansionSum(cde, bce, temp48a, cdeLen, bceLen);
	temp48bLen = arthemetricer.fastExpansionSum(deb, bcd, temp48b, debLen, bcdLen);
	for (int i = 0; i < temp48bLen; i++)
		temp48b[i] = -temp48b[i];
	xyzwLen = arthemetricer.fastExpansionSum(temp48a, temp48b, xyzw, temp48aLen, temp48bLen);
	xLen = arthemetricer.scaleExpansion(xyzw, a.x, temp192, xyzwLen);
	xLen = arthemetricer.scaleExpansion(temp192, a.x, temp384x, xLen);
	yLen = arthemetricer.scaleExpansion(xyzw, a.y, temp192, xyzwLen);
	yLen = arthemetricer.scaleExpansion(temp192, a.y, temp384y, yLen);
	zLen = arthemetricer.scaleExpansion(xyzw, a.z, temp192, xyzwLen);
	zLen = arthemetricer.scaleExpansion(temp192, a.z, temp384z, zLen);
	xyLen = arthemetricer.fastExpansionSum(temp384x, temp384y, temp768xy, xLen, yLen);
	adetLen = arthemetricer.fastExpansionSum(temp768xy, temp384z, adet, xyLen, zLen);

	FT bdet[1152];
	int bdetLen;
	temp48aLen = arthemetricer.fastExpansionSum(dea, cda, temp48a, deaLen, cdaLen);
	temp48bLen = arthemetricer.fastExpansionSum(cde, eac, temp48b, cdeLen, eacLen);
	for (int i = 0; i < temp48bLen; i++)
		temp48b[i] = -temp48b[i];
	xyzwLen = arthemetricer.fastExpansionSum(temp48a, temp48b, xyzw, temp48aLen, temp48bLen);
	xLen = arthemetricer.scaleExpansion(xyzw, b.x, temp192, xyzwLen);
	xLen = arthemetricer.scaleExpansion(temp192, b.x, temp384x, xLen);
	yLen = arthemetricer.scaleExpansion(xyzw, b.y, temp192, xyzwLen);
	yLen = arthemetricer.scaleExpansion(temp192, b.y, temp384y, yLen);
	zLen = arthemetricer.scaleExpansion(xyzw, b.z, temp192, xyzwLen);
	zLen = arthemetricer.scaleExpansion(temp192, b.z, temp384z, zLen);
	xyLen = arthemetricer.fastExpansionSum(temp384x, temp384y, temp768xy, xLen, yLen);
	bdetLen = arthemetricer.fastExpansionSum(temp768xy, temp384z, bdet, xyLen, zLen);

	FT cdet[1152];
	int cdetLen;
	temp48aLen = arthemetricer.fastExpansionSum(deb, eab, temp48a, debLen, eabLen);
	temp48bLen = arthemetricer.fastExpansionSum(dea, abd, temp48b, deaLen, abdLen);
	for (int i = 0; i < temp48bLen; i++)
		temp48b[i] = -temp48b[i];
	xyzwLen = arthemetricer.fastExpansionSum(temp48a, temp48b, xyzw, temp48aLen, temp48bLen);
	xLen = arthemetricer.scaleExpansion(xyzw, c.x, temp192, xyzwLen);
	xLen = arthemetricer.scaleExpansion(temp192, c.x, temp384x, xLen);
	yLen = arthemetricer.scaleExpansion(xyzw, c.y, temp192, xyzwLen);
	yLen = arthemetricer.scaleExpansion(temp192, c.y, temp384y, yLen);
	zLen = arthemetricer.scaleExpansion(xyzw, c.z, temp192, xyzwLen);
	zLen = arthemetricer.scaleExpansion(temp192, c.z, temp384z, zLen);
	xyLen = arthemetricer.fastExpansionSum(temp384x, temp384y, temp768xy, xLen, yLen);
	cdetLen = arthemetricer.fastExpansionSum(temp768xy, temp384z, cdet, xyLen, zLen);

	FT ddet[1152];
	int ddetLen;
	temp48aLen = arthemetricer.fastExpansionSum(eac, abc, temp48a, eacLen, abcLen);
	temp48bLen = arthemetricer.fastExpansionSum(bce, eab, temp48b, bceLen, eabLen);
	for (int i = 0; i < temp48bLen; i++)
		temp48b[i] = -temp48b[i];
	xyzwLen = arthemetricer.fastExpansionSum(temp48a, temp48b, xyzw, temp48aLen, temp48bLen);
	xLen = arthemetricer.scaleExpansion(xyzw, d.x, temp192, xyzwLen);
	xLen = arthemetricer.scaleExpansion(temp192, d.x, temp384x, xLen);
	yLen = arthemetricer.scaleExpansion(xyzw, d.y, temp192, xyzwLen);
	yLen = arthemetricer.scaleExpansion(temp192, d.y, temp384y, yLen);
	zLen = arthemetricer.scaleExpansion(xyzw, d.z, temp192, xyzwLen);
	zLen = arthemetricer.scaleExpansion(temp192, d.z, temp384z, zLen);
	xyLen = arthemetricer.fastExpansionSum(temp384x, temp384y, temp768xy, xLen, yLen);
	ddetLen = arthemetricer.fastExpansionSum(temp768xy, temp384z, ddet, xyLen, zLen);

	FT edet[1152];
	int edetLen;
	temp48aLen = arthemetricer.fastExpansionSum(bcd, abd, temp48a, bcdLen, abdLen);
	temp48bLen = arthemetricer.fastExpansionSum(cda, abc, temp48b, cdaLen, abcLen);
	for (int i = 0; i < temp48bLen; i++)
		temp48b[i] = -temp48b[i];
	xyzwLen = arthemetricer.fastExpansionSum(temp48a, temp48b, xyzw, temp48aLen, temp48bLen);
	xLen = arthemetricer.scaleExpansion(xyzw, e.x, temp192, xyzwLen);
	xLen = arthemetricer.scaleExpansion(temp192, e.x, temp384x, xLen);
	yLen = arthemetricer.scaleExpansion(xyzw, e.y, temp192, xyzwLen);
	yLen = arthemetricer.scaleExpansion(temp192, e.y, temp384y, yLen);
	zLen = arthemetricer.scaleExpansion(xyzw, e.z, temp192, xyzwLen);
	zLen = arthemetricer.scaleExpansion(temp192, e.z, temp384z, zLen);
	xyLen = arthemetricer.fastExpansionSum(temp384x, temp384y, temp768xy, xLen, yLen);
	edetLen = arthemetricer.fastExpansionSum(temp768xy, temp384z, edet, xyLen, zLen);

	FT abdet[2304], cddet[2304], cdedet[3456], fin[5760];
	int abdetLen, cddetLen, cdedetLen, finLength;
	abdetLen = arthemetricer.fastExpansionSum(adet, bdet, abdet, adetLen, bdetLen);
	cddetLen = arthemetricer.fastExpansionSum(cdet, ddet, cddet, cdetLen, ddetLen);
	cdedetLen = arthemetricer.fastExpansionSum(cddet, edet, cdedet, cddetLen, edetLen);
	finLength = arthemetricer.fastExpansionSum(cdedet, abdet, fin, cdedetLen, abdetLen);

	return fin[finLength - 1];
}

template<class FT> FT Predicator<FT>::inOrthoSphere(const VectorBase<FT> &a, FT aWeight, const VectorBase<FT> &b, FT bWeight,
	const VectorBase<FT> &c, FT cWeight, const VectorBase<FT> &d, FT dWeight, const VectorBase<FT> &e, FT eWeight) const{
	VectorBase<FT> ea = a - e;
	VectorBase<FT> eb = b - e;
	VectorBase<FT> ec = c - e;
	VectorBase<FT> ed = d - e;

	FT eaxeby = ea.x*eb.y;
	FT eayebx = ea.y*eb.x;
	FT ab = eaxeby - eayebx;

	FT ebxecy = eb.x*ec.y;
	FT ebyecx = eb.y*ec.x;
	FT bc = ebxecy - ebyecx;

	FT ecxedy = ec.x*ed.y;
	FT ecyedx = ec.y*ed.x;
	FT cd = ecxedy - ecyedx;

	FT edxeay = ed.x*ea.y;
	FT edyeax = ed.y*ea.x;
	FT da = edxeay - edyeax;

	FT eaxecy = ea.x*ec.y;
	FT eayecx = ea.y*ec.x;
	FT ac = eaxecy - eayecx;

	FT ebxedy = eb.x*ed.y;
	FT ebyedx = eb.y*ed.x;
	FT bd = ebxedy - ebyedx;

	FT abc = ea.z*bc - eb.z*ac + ec.z*ab;
	FT bcd = eb.z*cd - ec.z*bd + ed.z*bc;
	FT cda = ec.z*da + ed.z*ac + ea.z*cd;
	FT dab = ed.z*ab + ea.z*bd + eb.z*da;

	FT aLift = ea.length2() + eWeight - aWeight;
	FT bLift = eb.length2() + eWeight - bWeight;
	FT cLift = ec.length2() + eWeight - cWeight;
	FT dLift = ed.length2() + eWeight - dWeight;

	FT det = -aLift*bcd + bLift*cda - cLift*dab + dLift*abc;

	FT eaxebyAbsolute = fabs(eaxeby);
	FT eayebxAbsolute = fabs(eayebx);
	FT ebxecyAbsolute = fabs(ebxecy);
	FT ebyecxAbsolute = fabs(ebyecx);
	FT ecxedyAbsolute = fabs(ecxedy);
	FT ecyedxAbsolute = fabs(ecyedx);
	FT edxeayAbsolute = fabs(edxeay);
	FT edyeaxAbsolute = fabs(edyeax);
	FT eaxecyAbsolute = fabs(eaxecy);
	FT eayecxAbsolute = fabs(eayecx);
	FT ebxedyAbsolute = fabs(ebxedy);
	FT ebyedxAbsolute = fabs(ebyedx);
	FT eazAbsolute = fabs(ea.z);
	FT ebzAbsolute = fabs(eb.z);
	FT eczAbsolute = fabs(ec.z);
	FT edzAbsolute = fabs(ed.z);

	abc = eazAbsolute*(ebxecyAbsolute + ebyecxAbsolute) + ebzAbsolute*(ecxedyAbsolute + eayecxAbsolute) + eczAbsolute*(eaxebyAbsolute + eayebxAbsolute);
	bcd = ebzAbsolute*(ecxedyAbsolute + ecyedxAbsolute) + ebzAbsolute*(ebxedyAbsolute + ebyedxAbsolute) + edzAbsolute*(ebxecyAbsolute + ebyecxAbsolute);
	cda = eczAbsolute*(edxeayAbsolute + edyeaxAbsolute) + edzAbsolute*(eaxecyAbsolute + eayecxAbsolute) + eazAbsolute*(ecxedyAbsolute + ecyedxAbsolute);
	dab = edzAbsolute*(eaxebyAbsolute + eayebxAbsolute) + eazAbsolute*(ebxedyAbsolute + ebyedxAbsolute) + ebzAbsolute*(edxeayAbsolute + edyeaxAbsolute);

	FT norm = fabs(aLift)*bcd + fabs(bLift)*cda + fabs(cLift)*dab + fabs(dLift)*abc;
	FT errorBound = inSpeErrorBoundA*norm;
	if (fabs(det) >= errorBound)
		return det;

	return inOrthoSphereAdaptive(a, aWeight, b, bWeight, c, cWeight, d, dWeight, e, eWeight, norm);
}

template<class FT> FT Predicator<FT>::inOrthoSphereAdaptive(const VectorBase<FT> &a, FT aWeight, const VectorBase<FT> &b, FT bWeight,
	const VectorBase<FT> &c, FT cWeight, const VectorBase<FT> &d, FT dWeight, const VectorBase<FT> &e, FT eWeight, FT norm) const{
	FT fin[1344];
	int finLength;
	FT det, errorBound;

	VectorBase<FT> ea = a - e;
	VectorBase<FT> eb = b - e;
	VectorBase<FT> ec = c - e;
	VectorBase<FT> ed = d - e;

	FT eaWeight = eWeight - aWeight;
	FT ebWeight = eWeight - bWeight;
	FT ecWeight = eWeight - cWeight;
	FT edWeight = eWeight - dWeight;

	FT eaxeby1, eaxeby0, eayebx1, eayebx0, ab[4];
	arthemetricer.twoProduct(ea.x, eb.y, eaxeby1, eaxeby0);
	arthemetricer.twoProduct(ea.y, eb.x, eayebx1, eayebx0);
	arthemetricer.twoTwoDiff(eaxeby1, eaxeby0, eayebx1, eayebx0, ab[3], ab[2], ab[1], ab[0]);

	FT ebxecy1, ebxecy0, ebyecx1, ebyecx0, bc[4];
	arthemetricer.twoProduct(eb.x, ec.y, ebxecy1, ebxecy0);
	arthemetricer.twoProduct(eb.y, ec.x, ebyecx1, ebyecx0);
	arthemetricer.twoTwoDiff(ebxecy1, ebxecy0, ebyecx1, ebyecx0, bc[3], bc[2], bc[1], bc[0]);

	FT ecxedy1, ecxedy0, ecyedx1, ecyedx0, cd[4];
	arthemetricer.twoProduct(ec.x, ed.y, ecxedy1, ecxedy0);
	arthemetricer.twoProduct(ec.y, ed.x, ecyedx1, ecyedx0);
	arthemetricer.twoTwoDiff(ecxedy1, ecxedy0, ecyedx1, ecyedx0, cd[3], cd[2], cd[1], cd[0]);

	FT edxeay1, edxeay0, edyeax1, edyeax0, da[4];
	arthemetricer.twoProduct(ed.x, ea.y, edxeay1, edxeay0);
	arthemetricer.twoProduct(ed.y, ea.x, edyeax1, edyeax0);
	arthemetricer.twoTwoDiff(edxeay1, edxeay0, edyeax1, edyeax0, da[3], da[2], da[1], da[0]);

	FT eaxecy1, eaxecy0, eayecx1, eayecx0, ac[4];
	arthemetricer.twoProduct(ea.x, ec.y, eaxecy1, eaxecy0);
	arthemetricer.twoProduct(ea.y, ec.x, eayecx1, eayecx0);
	arthemetricer.twoTwoDiff(eaxecy1, eaxecy0, eayecx1, eayecx0, ac[3], ac[2], ac[1], ac[0]);

	FT ebxedy1, ebxedy0, ebyedx1, ebyedx0, bd[4];
	arthemetricer.twoProduct(eb.x, ed.y, ebxedy1, ebxedy0);
	arthemetricer.twoProduct(eb.y, ed.x, ebyedx1, ebyedx0);
	arthemetricer.twoTwoDiff(ebxedy1, ebxedy0, ebyedx1, ebyedx0, bd[3], bd[2], bd[1], bd[0]);

	FT temp8a[8], temp8b[8], temp8c[8], temp16[16], temp24[24], temp48[48];
	int temp8aLen, temp8bLen, temp8cLen, temp16Len, temp24Len, temp48Len;
	FT xdet[96], ydet[96], zdet[96], wdet[48], xydet[192], xyzdet[288];
	int xdetLen, ydetLen, zdetLen, wdetLen, xydetLen, xyzdetLen;
	FT adet[336], bdet[336], cdet[336], ddet[336];
	int adetLen, bdetLen, cdetLen, ddetLen;

	//abc = ea.z*bc - eb.z*ac + ec.z*ab;
	temp8aLen = arthemetricer.scaleExpansion(bc, ea.z, temp8a, 4);
	temp8bLen = arthemetricer.scaleExpansion(ac, -eb.z, temp8b, 4);
	temp8cLen = arthemetricer.scaleExpansion(ab, ec.z, temp8c, 4);
	temp16Len = arthemetricer.fastExpansionSum(temp8a, temp8b, temp16, temp8aLen, temp8bLen);
	temp24Len = arthemetricer.fastExpansionSum(temp8c, temp16, temp24, temp8cLen, temp16Len);

	//(ed.leng2() + edWeight)*abc
	temp48Len = arthemetricer.scaleExpansion(temp24, ed.x, temp48, temp24Len);
	xdetLen = arthemetricer.scaleExpansion(temp48, ed.x, xdet, temp48Len);
	temp48Len = arthemetricer.scaleExpansion(temp24, ed.y, temp48, temp24Len);
	ydetLen = arthemetricer.scaleExpansion(temp48, ed.y, ydet, temp48Len);
	temp48Len = arthemetricer.scaleExpansion(temp24, ed.z, temp48, temp24Len);
	zdetLen = arthemetricer.scaleExpansion(temp48, ed.z, zdet, temp48Len);
	wdetLen = arthemetricer.scaleExpansion(temp24, edWeight, wdet, temp24Len);
	xydetLen = arthemetricer.fastExpansionSum(xdet, ydet, xydet, xdetLen, ydetLen);
	xyzdetLen = arthemetricer.fastExpansionSum(xydet, zdet, xyzdet, xydetLen, zdetLen);
	ddetLen = arthemetricer.fastExpansionSum(xyzdet, wdet, ddet, xyzdetLen, wdetLen);

	//bcd = eb.z*cd - ec.z*bd + ed.z*bc
	temp8aLen = arthemetricer.scaleExpansion(cd, eb.z, temp8a, 4);
	temp8bLen = arthemetricer.scaleExpansion(bd, -ec.z, temp8b, 4);
	temp8cLen = arthemetricer.scaleExpansion(bc, ed.z, temp8c, 4);
	temp16Len = arthemetricer.fastExpansionSum(temp8a, temp8b, temp16, temp8aLen, temp8bLen);
	temp24Len = arthemetricer.fastExpansionSum(temp8c, temp16, temp24, temp8cLen, temp16Len);

	//-(ea.length2() + eaWeight)*bcd
	temp48Len = arthemetricer.scaleExpansion(temp24, ea.x, temp48, temp24Len);
	xdetLen = arthemetricer.scaleExpansion(temp48, -ea.x, xdet, temp48Len);
	temp48Len = arthemetricer.scaleExpansion(temp24, ea.y, temp48, temp24Len);
	ydetLen = arthemetricer.scaleExpansion(temp48, -ea.y, ydet, temp48Len);
	temp48Len = arthemetricer.scaleExpansion(temp24, ea.z, temp48, temp24Len);
	zdetLen = arthemetricer.scaleExpansion(temp48, -ea.z, zdet, temp48Len);
	wdetLen = arthemetricer.scaleExpansion(temp24, -eaWeight, wdet, temp24Len);
	xydetLen = arthemetricer.fastExpansionSum(xdet, ydet, xydet, xdetLen, ydetLen);
	xyzdetLen = arthemetricer.fastExpansionSum(xydet, zdet, xyzdet, xydetLen, zdetLen);
	adetLen = arthemetricer.fastExpansionSum(xyzdet, wdet, adet, xyzdetLen, wdetLen);

	//cda = ec.z*da + ed.z*ac + ea.z*cd
	temp8aLen = arthemetricer.scaleExpansion(da, ec.z, temp8a, 4);
	temp8bLen = arthemetricer.scaleExpansion(ac, ed.z, temp8b, 4);
	temp8cLen = arthemetricer.scaleExpansion(cd, ea.z, temp8c, 4);
	temp16Len = arthemetricer.fastExpansionSum(temp8a, temp8b, temp16, temp8aLen, temp8bLen);
	temp24Len = arthemetricer.fastExpansionSum(temp8c, temp16, temp24, temp8cLen, temp16Len);

	//(eb.length2() + ebWeight)*cda
	temp48Len = arthemetricer.scaleExpansion(temp24, eb.x, temp48, temp24Len);
	xdetLen = arthemetricer.scaleExpansion(temp48, eb.x, xdet, temp48Len);
	temp48Len = arthemetricer.scaleExpansion(temp24, eb.y, temp48, temp24Len);
	ydetLen = arthemetricer.scaleExpansion(temp48, eb.y, ydet, temp48Len);
	temp48Len = arthemetricer.scaleExpansion(temp24, eb.z, temp48, temp24Len);
	zdetLen = arthemetricer.scaleExpansion(temp48, eb.z, zdet, temp48Len);
	wdetLen = arthemetricer.scaleExpansion(temp24, ebWeight, wdet, temp24Len);
	xydetLen = arthemetricer.fastExpansionSum(xdet, ydet, xydet, xdetLen, ydetLen);
	xyzdetLen = arthemetricer.fastExpansionSum(xydet, zdet, xyzdet, xydetLen, zdetLen);
	bdetLen = arthemetricer.fastExpansionSum(xyzdet, wdet, bdet, xyzdetLen, wdetLen);

	//FT dab = ed.z*ab + ea.z*bd + eb.z*da
	temp8aLen = arthemetricer.scaleExpansion(ab, ed.z, temp8a, 4);
	temp8bLen = arthemetricer.scaleExpansion(bd, ea.z, temp8b, 4);
	temp8cLen = arthemetricer.scaleExpansion(da, eb.z, temp8c, 4);
	temp16Len = arthemetricer.fastExpansionSum(temp8a, temp8b, temp16, temp8aLen, temp8bLen);
	temp24Len = arthemetricer.fastExpansionSum(temp8c, temp16, temp24, temp8cLen, temp16Len);

	//-(ec.length2() + ecWeight)*dab
	temp48Len = arthemetricer.scaleExpansion(temp24, ec.x, temp48, temp24Len);
	xdetLen = arthemetricer.scaleExpansion(temp48, -ec.x, xdet, temp48Len);
	temp48Len = arthemetricer.scaleExpansion(temp24, ec.y, temp48, temp24Len);
	ydetLen = arthemetricer.scaleExpansion(temp48, -ec.y, ydet, temp48Len);
	temp48Len = arthemetricer.scaleExpansion(temp24, ec.z, temp48, temp24Len);
	zdetLen = arthemetricer.scaleExpansion(temp48, -ec.z, zdet, temp48Len);
	wdetLen = arthemetricer.scaleExpansion(temp24, -ecWeight, wdet, temp24Len);
	xydetLen = arthemetricer.fastExpansionSum(xdet, ydet, xydet, xdetLen, ydetLen);
	xyzdetLen = arthemetricer.fastExpansionSum(xydet, zdet, xyzdet, xydetLen, zdetLen);
	cdetLen = arthemetricer.fastExpansionSum(xyzdet, wdet, cdet, xyzdetLen, wdetLen);

	FT abdet[672], cddet[672];
	int abdetLen, cddetLen;
	abdetLen = arthemetricer.fastExpansionSum(adet, bdet, abdet, adetLen, bdetLen);
	cddetLen = arthemetricer.fastExpansionSum(cdet, ddet, cddet, cdetLen, ddetLen);
	finLength = arthemetricer.fastExpansionSum(abdet, cddet, fin, abdetLen, cddetLen);

	det = arthemetricer.Estimate(fin, finLength);
	errorBound = inSpeErrorBoundB*norm;
	if (fabs(det) >= errorBound)
		return det;

	FT eaxError = arthemetricer.twoDiffError(a.x, e.x, ea.x);
	FT eayError = arthemetricer.twoDiffError(a.y, e.y, ea.y);
	FT eazError = arthemetricer.twoDiffError(a.z, e.z, ea.z);
	FT eawError = arthemetricer.twoDiffError(eWeight, aWeight, eaWeight);

	FT ebxError = arthemetricer.twoDiffError(b.x, e.x, eb.x);
	FT ebyError = arthemetricer.twoDiffError(b.y, e.y, eb.y);
	FT ebzError = arthemetricer.twoDiffError(b.z, e.z, eb.z);
	FT ebwError = arthemetricer.twoDiffError(eWeight, bWeight, ebWeight);

	FT ecxError = arthemetricer.twoDiffError(c.x, e.x, ec.x);
	FT ecyError = arthemetricer.twoDiffError(c.y, e.y, ec.y);
	FT eczError = arthemetricer.twoDiffError(c.z, e.z, ec.z);
	FT ecwError = arthemetricer.twoDiffError(eWeight, cWeight, ecWeight);

	FT edxError = arthemetricer.twoDiffError(d.x, e.x, ed.x);
	FT edyError = arthemetricer.twoDiffError(d.y, e.y, ed.y);
	FT edzError = arthemetricer.twoDiffError(d.z, e.z, ed.z);
	FT edwError = arthemetricer.twoDiffError(eWeight, dWeight, edWeight);

	if (eaxError == 0.0 && eayError == 0.0 && eazError == 0.0 && eawError == 0.0
		&& ebxError == 0.0 && ebyError == 0.0 && ebzError == 0.0 && ebwError == 0.0
		&& ecxError == 0.0 && ecyError == 0.0 && eczError == 0.0 && ecwError == 0.0
		&& edxError == 0.0 && edyError == 0.0 && edzError == 0.0 && edwError == 0.0)
		return det;

	FT abError = (ea.x*ebyError + eb.y*eaxError) - (ea.y*ebxError + eb.x*eayError);
	FT bcError = (eb.x*ecyError + ec.y*ebxError) - (eb.y*ecxError + ec.x*ebyError);
	FT cdError = (ec.x*edyError + ed.y*ecxError) - (ec.y*edxError + ed.x*ecyError);
	FT daError = (ed.x*eayError + ea.y*edxError) - (ed.y*eaxError + ea.x*edyError);
	FT acError = (ea.x*ecyError + ec.y*eaxError) - (ea.y*ecxError + ec.x*eayError);
	FT bdError = (eb.x*edyError + ed.y*ebxError) - (eb.y*edxError + ed.x*ebyError);

	FT abc = ea.z*bc[3] - eb.z*ac[3] + ec.z*ab[3];
	FT bcd = eb.z*cd[3] - ec.z*bd[3] + ed.z*bc[3];
	FT cda = ec.z*da[3] + ed.z*ac[3] + ea.z*cd[3];
	FT dab = ed.z*ab[3] + ea.z*bd[3] + eb.z*da[3];

	FT abcError = (ea.z*bcError - eb.z*acError + ec.z*abError) + (eazError*bc[3] - ebzError*ac[3] + eczError*ab[3]);
	FT bcdError = (eb.z*cdError - ec.z*bdError + ed.z*bcError) + (ebzError*cd[3] - eczError*bd[3] + edzError*bc[3]);
	FT cdaError = (ec.z*daError + ed.z*acError + ea.z*cdError) + (eczError*da[3] + edzError*ac[3] + eazError*cd[3]);
	FT dabError = (ed.z*abError + ea.z*bdError + eb.z*daError) + (edzError*ab[3] + eazError*bd[3] + ebzError*da[3]);

	FT eaLength2ErrorOne = ea.x*eaxError + ea.y*eayError + ea.z*eazError;
	FT ebLength2ErrorOne = eb.x*ebxError + eb.y*ebyError + eb.z*ebzError;
	FT ecLength2ErrorOne = ec.x*ecxError + ec.y*ecyError + ec.z*eczError;
	FT edLength2ErrorOne = ed.x*edxError + ed.y*edyError + ed.z*edzError;


	det += ((-(ea.length2() + eaWeight)*bcdError + (eb.length2() + ebWeight)*cdaError
		- (ec.length2() + ecWeight)*dabError + (ed.length2() + edWeight)*abcError)
		+ (-(FT(2)*eaLength2ErrorOne + eawError)*bcd + (FT(2)*ebLength2ErrorOne + ebwError)*cda
		- (FT(2)*ecLength2ErrorOne + ecwError)*dab + (FT(2)*edLength2ErrorOne + edwError)*abc));

	errorBound = inSpeErrorBoundC*norm + commonBound*fabs(det);
	if (fabs(det) >= errorBound)
		return det;

	return inOrthoSphereExact(a, aWeight, b, bWeight, c, cWeight, d, dWeight, e, eWeight);
}

template<class FT> FT Predicator<FT>::inOrthoSphereExact(const VectorBase<FT> &a, FT aWeight, const VectorBase<FT> &b, FT bWeight,
	const VectorBase<FT> &c, FT cWeight, const VectorBase<FT> &d, FT dWeight, const VectorBase<FT> &e, FT eWeight) const{
	FT axby1, axby0, aybx1, aybx0, ab[4];
	arthemetricer.twoProduct(a.x, b.y, axby1, axby0);
	arthemetricer.twoProduct(a.y, b.x, aybx1, aybx0);
	arthemetricer.twoTwoDiff(axby1, axby0, aybx1, aybx0, ab[3], ab[2], ab[1], ab[0]);

	FT bxcy1, bxcy0, bycx1, bycx0, bc[4];
	arthemetricer.twoProduct(b.x, c.y, bxcy1, bxcy0);
	arthemetricer.twoProduct(b.y, c.x, bycx1, bycx0);
	arthemetricer.twoTwoDiff(bxcy1, bxcy0, bycx1, bycx0, bc[3], bc[2], bc[1], bc[0]);

	FT cxdy1, cxdy0, cydx1, cydx0, cd[4];
	arthemetricer.twoProduct(c.x, d.y, cxdy1, cxdy0);
	arthemetricer.twoProduct(c.y, d.x, cydx1, cydx0);
	arthemetricer.twoTwoDiff(cxdy1, cxdy0, cydx1, cydx0, cd[3], cd[2], cd[1], cd[0]);

	FT dxey1, dxey0, dyex1, dyex0, de[4];
	arthemetricer.twoProduct(d.x, e.y, dxey1, dxey0);
	arthemetricer.twoProduct(d.y, e.x, dyex1, dyex0);
	arthemetricer.twoTwoDiff(dxey1, dxey0, dyex1, dyex0, de[3], de[2], de[1], de[0]);

	FT exay1, exay0, eyax1, eyax0, ea[4];
	arthemetricer.twoProduct(e.x, a.y, exay1, exay0);
	arthemetricer.twoProduct(e.y, a.x, eyax1, eyax0);
	arthemetricer.twoTwoDiff(exay1, exay0, eyax1, eyax0, ea[3], ea[2], ea[1], ea[0]);

	FT axcy1, axcy0, aycx1, aycx0, ac[4];
	arthemetricer.twoProduct(a.x, c.y, axcy1, axcy0);
	arthemetricer.twoProduct(a.y, c.x, aycx1, aycx0);
	arthemetricer.twoTwoDiff(axcy1, axcy0, aycx1, aycx0, ac[3], ac[2], ac[1], ac[0]);

	FT bxdy1, bxdy0, bydx1, bydx0, bd[4];
	arthemetricer.twoProduct(b.x, d.y, bxdy1, bxdy0);
	arthemetricer.twoProduct(b.y, d.x, bydx1, bydx0);
	arthemetricer.twoTwoDiff(bxdy1, bxdy0, bydx1, bydx0, bd[3], bd[2], bd[1], bd[0]);

	FT cxey1, cxey0, cyex1, cyex0, ce[4];
	arthemetricer.twoProduct(c.x, e.y, cxey1, cxey0);
	arthemetricer.twoProduct(c.y, e.x, cyex1, cyex0);
	arthemetricer.twoTwoDiff(cxey1, cxey0, cyex1, cyex0, ce[3], ce[2], ce[1], ce[0]);

	FT dxay1, dxay0, dyax1, dyax0, da[4];
	arthemetricer.twoProduct(d.x, a.y, dxay1, dxay0);
	arthemetricer.twoProduct(d.y, a.x, dyax1, dyax0);
	arthemetricer.twoTwoDiff(dxay1, dxay0, dyax1, dyax0, da[3], da[2], da[1], da[0]);

	FT exby1, exby0, eybx1, eybx0, eb[4];
	arthemetricer.twoProduct(e.x, b.y, exby1, exby0);
	arthemetricer.twoProduct(e.y, b.x, eybx1, eybx0);
	arthemetricer.twoTwoDiff(exby1, exby0, eybx1, eybx0, eb[3], eb[2], eb[1], eb[0]);


	FT abc[24], bcd[24], cde[24], dea[24], eab[24], abd[24], bce[24], cda[24], deb[24], eac[24];
	int abcLen, bcdLen, cdeLen, deaLen, eabLen, abdLen, bceLen, cdaLen, debLen, eacLen;
	FT temp8a[8], temp8b[9], temp16[16];
	int temp8aLen, temp8bLen, temp16Len;

	//abc
	temp8aLen = arthemetricer.scaleExpansion(bc, a.z, temp8a, 4);
	temp8bLen = arthemetricer.scaleExpansion(ac, -b.z, temp8b, 4);
	temp16Len = arthemetricer.fastExpansionSum(temp8a, temp8b, temp16, temp8aLen, temp8bLen);
	temp8aLen = arthemetricer.scaleExpansion(ab, c.z, temp8a, 4);
	abcLen = arthemetricer.fastExpansionSum(temp16, temp8a, abc, temp16Len, temp8aLen);

	//bcd
	temp8aLen = arthemetricer.scaleExpansion(cd, b.z, temp8a, 4);
	temp8bLen = arthemetricer.scaleExpansion(bd, -c.z, temp8b, 4);
	temp16Len = arthemetricer.fastExpansionSum(temp8a, temp8b, temp16, temp8aLen, temp8bLen);
	temp8aLen = arthemetricer.scaleExpansion(bc, d.z, temp8a, 4);
	bcdLen = arthemetricer.fastExpansionSum(temp16, temp8a, bcd, temp16Len, temp8aLen);

	//cde
	temp8aLen = arthemetricer.scaleExpansion(de, c.z, temp8a, 4);
	temp8bLen = arthemetricer.scaleExpansion(ce, -d.z, temp8b, 4);
	temp16Len = arthemetricer.fastExpansionSum(temp8a, temp8b, temp16, temp8aLen, temp8bLen);
	temp8aLen = arthemetricer.scaleExpansion(cd, e.z, temp8a, 4);
	cdeLen = arthemetricer.fastExpansionSum(temp16, temp8a, cde, temp16Len, temp8aLen);

	//dea
	temp8aLen = arthemetricer.scaleExpansion(ea, d.z, temp8a, 4);
	temp8bLen = arthemetricer.scaleExpansion(da, -e.z, temp8b, 4);
	temp16Len = arthemetricer.fastExpansionSum(temp8a, temp8b, temp16, temp8aLen, temp8bLen);
	temp8aLen = arthemetricer.scaleExpansion(de, a.z, temp8a, 4);
	deaLen = arthemetricer.fastExpansionSum(temp16, temp8a, dea, temp16Len, temp8aLen);

	//eab
	temp8aLen = arthemetricer.scaleExpansion(ab, e.z, temp8a, 4);
	temp8bLen = arthemetricer.scaleExpansion(eb, -a.z, temp8b, 4);
	temp16Len = arthemetricer.fastExpansionSum(temp8a, temp8b, temp16, temp8aLen, temp8bLen);
	temp8aLen = arthemetricer.scaleExpansion(ea, b.z, temp8a, 4);
	eabLen = arthemetricer.fastExpansionSum(temp16, temp8a, eab, temp16Len, temp8aLen);

	//abd
	temp8aLen = arthemetricer.scaleExpansion(bd, a.z, temp8a, 4);
	temp8bLen = arthemetricer.scaleExpansion(da, b.z, temp8b, 4);
	temp16Len = arthemetricer.fastExpansionSum(temp8a, temp8b, temp16, temp8aLen, temp8bLen);
	temp8aLen = arthemetricer.scaleExpansion(ab, d.z, temp8a, 4);
	abdLen = arthemetricer.fastExpansionSum(temp16, temp8a, abd, temp16Len, temp8aLen);

	//bce
	temp8aLen = arthemetricer.scaleExpansion(ce, b.z, temp8a, 4);
	temp8bLen = arthemetricer.scaleExpansion(eb, c.z, temp8b, 4);
	temp16Len = arthemetricer.fastExpansionSum(temp8a, temp8b, temp16, temp8aLen, temp8bLen);
	temp8aLen = arthemetricer.scaleExpansion(bc, e.z, temp8a, 4);
	bceLen = arthemetricer.fastExpansionSum(temp16, temp8a, bce, temp16Len, temp8aLen);

	//cda
	temp8aLen = arthemetricer.scaleExpansion(da, c.z, temp8a, 4);
	temp8bLen = arthemetricer.scaleExpansion(ac, d.z, temp8b, 4);
	temp16Len = arthemetricer.fastExpansionSum(temp8a, temp8b, temp16, temp8aLen, temp8bLen);
	temp8aLen = arthemetricer.scaleExpansion(cd, a.z, temp8a, 4);
	cdaLen = arthemetricer.fastExpansionSum(temp16, temp8a, cda, temp16Len, temp8aLen);

	//deb
	temp8aLen = arthemetricer.scaleExpansion(eb, d.z, temp8a, 4);
	temp8bLen = arthemetricer.scaleExpansion(bd, e.z, temp8b, 4);
	temp16Len = arthemetricer.fastExpansionSum(temp8a, temp8b, temp16, temp8aLen, temp8bLen);
	temp8aLen = arthemetricer.scaleExpansion(de, b.z, temp8a, 4);
	debLen = arthemetricer.fastExpansionSum(temp16, temp8a, deb, temp16Len, temp8aLen);

	//eac
	temp8aLen = arthemetricer.scaleExpansion(ac, e.z, temp8a, 4);
	temp8bLen = arthemetricer.scaleExpansion(ce, a.z, temp8b, 4);
	temp16Len = arthemetricer.fastExpansionSum(temp8a, temp8b, temp16, temp8aLen, temp8bLen);
	temp8aLen = arthemetricer.scaleExpansion(ea, c.z, temp8a, 4);
	eacLen = arthemetricer.fastExpansionSum(temp16, temp8a, eac, temp16Len, temp8aLen);

	FT xyzw[96];
	int xyzwLen;
	FT temp48a[48], temp48b[48], temp192[192], temp192w[192], temp384x[384], temp384y[384], temp384z[384], temp768xy[768], temp1152xyz[1152];
	int temp48aLen, temp48bLen, xLen, yLen, zLen, wLen, xyLen, xyzLen;

	FT adet[1344];
	int adetLen;
	temp48aLen = arthemetricer.fastExpansionSum(cde, bce, temp48a, cdeLen, bceLen);
	temp48bLen = arthemetricer.fastExpansionSum(deb, bcd, temp48b, debLen, bcdLen);
	for (int i = 0; i < temp48bLen; i++)
		temp48b[i] = -temp48b[i];
	xyzwLen = arthemetricer.fastExpansionSum(temp48a, temp48b, xyzw, temp48aLen, temp48bLen);
	xLen = arthemetricer.scaleExpansion(xyzw, a.x, temp192, xyzwLen);
	xLen = arthemetricer.scaleExpansion(temp192, a.x, temp384x, xLen);
	yLen = arthemetricer.scaleExpansion(xyzw, a.y, temp192, xyzwLen);
	yLen = arthemetricer.scaleExpansion(temp192, a.y, temp384y, yLen);
	zLen = arthemetricer.scaleExpansion(xyzw, a.z, temp192, xyzwLen);
	zLen = arthemetricer.scaleExpansion(temp192, a.z, temp384z, zLen);
	wLen = arthemetricer.scaleExpansion(xyzw, -aWeight, temp192w, xyzwLen);

	xyLen = arthemetricer.fastExpansionSum(temp384x, temp384y, temp768xy, xLen, yLen);
	xyzLen = arthemetricer.fastExpansionSum(temp768xy, temp384z, temp1152xyz, xyLen, zLen);
	adetLen = arthemetricer.fastExpansionSum(temp1152xyz, temp192w, adet, xyzLen, wLen);

	FT bdet[1344];
	int bdetLen;
	temp48aLen = arthemetricer.fastExpansionSum(dea, cda, temp48a, deaLen, cdaLen);
	temp48bLen = arthemetricer.fastExpansionSum(cde, eac, temp48b, cdeLen, eacLen);
	for (int i = 0; i < temp48bLen; i++)
		temp48b[i] = -temp48b[i];
	xyzwLen = arthemetricer.fastExpansionSum(temp48a, temp48b, xyzw, temp48aLen, temp48bLen);
	xLen = arthemetricer.scaleExpansion(xyzw, b.x, temp192, xyzwLen);
	xLen = arthemetricer.scaleExpansion(temp192, b.x, temp384x, xLen);
	yLen = arthemetricer.scaleExpansion(xyzw, b.y, temp192, xyzwLen);
	yLen = arthemetricer.scaleExpansion(temp192, b.y, temp384y, yLen);
	zLen = arthemetricer.scaleExpansion(xyzw, b.z, temp192, xyzwLen);
	zLen = arthemetricer.scaleExpansion(temp192, b.z, temp384z, zLen);
	wLen = arthemetricer.scaleExpansion(xyzw, -bWeight, temp192w, xyzwLen);

	xyLen = arthemetricer.fastExpansionSum(temp384x, temp384y, temp768xy, xLen, yLen);
	xyzLen = arthemetricer.fastExpansionSum(temp768xy, temp384z, temp1152xyz, xyLen, zLen);
	bdetLen = arthemetricer.fastExpansionSum(temp1152xyz, temp192w, bdet, xyzLen, wLen);

	FT cdet[1344];
	int cdetLen;
	temp48aLen = arthemetricer.fastExpansionSum(deb, eab, temp48a, debLen, eabLen);
	temp48bLen = arthemetricer.fastExpansionSum(dea, abd, temp48b, deaLen, abdLen);
	for (int i = 0; i < temp48bLen; i++)
		temp48b[i] = -temp48b[i];
	xyzwLen = arthemetricer.fastExpansionSum(temp48a, temp48b, xyzw, temp48aLen, temp48bLen);
	xLen = arthemetricer.scaleExpansion(xyzw, c.x, temp192, xyzwLen);
	xLen = arthemetricer.scaleExpansion(temp192, c.x, temp384x, xLen);
	yLen = arthemetricer.scaleExpansion(xyzw, c.y, temp192, xyzwLen);
	yLen = arthemetricer.scaleExpansion(temp192, c.y, temp384y, yLen);
	zLen = arthemetricer.scaleExpansion(xyzw, c.z, temp192, xyzwLen);
	zLen = arthemetricer.scaleExpansion(temp192, c.z, temp384z, zLen);
	wLen = arthemetricer.scaleExpansion(xyzw, -cWeight, temp192w, xyzwLen);

	xyLen = arthemetricer.fastExpansionSum(temp384x, temp384y, temp768xy, xLen, yLen);
	xyzLen = arthemetricer.fastExpansionSum(temp768xy, temp384z, temp1152xyz, xyLen, zLen);
	cdetLen = arthemetricer.fastExpansionSum(temp1152xyz, temp192w, cdet, xyzLen, wLen);

	FT ddet[1344];
	int ddetLen;
	temp48aLen = arthemetricer.fastExpansionSum(eac, abc, temp48a, eacLen, abcLen);
	temp48bLen = arthemetricer.fastExpansionSum(bce, eab, temp48b, bceLen, eabLen);
	for (int i = 0; i < temp48bLen; i++)
		temp48b[i] = -temp48b[i];
	xyzwLen = arthemetricer.fastExpansionSum(temp48a, temp48b, xyzw, temp48aLen, temp48bLen);
	xLen = arthemetricer.scaleExpansion(xyzw, d.x, temp192, xyzwLen);
	xLen = arthemetricer.scaleExpansion(temp192, d.x, temp384x, xLen);
	yLen = arthemetricer.scaleExpansion(xyzw, d.y, temp192, xyzwLen);
	yLen = arthemetricer.scaleExpansion(temp192, d.y, temp384y, yLen);
	zLen = arthemetricer.scaleExpansion(xyzw, d.z, temp192, xyzwLen);
	zLen = arthemetricer.scaleExpansion(temp192, d.z, temp384z, zLen);
	wLen = arthemetricer.scaleExpansion(xyzw, -dWeight, temp192w, xyzwLen);

	xyLen = arthemetricer.fastExpansionSum(temp384x, temp384y, temp768xy, xLen, yLen);
	xyzLen = arthemetricer.fastExpansionSum(temp768xy, temp384z, temp1152xyz, xyLen, zLen);
	ddetLen = arthemetricer.fastExpansionSum(temp1152xyz, temp192w, ddet, xyzLen, wLen);

	FT edet[1344];
	int edetLen;
	temp48aLen = arthemetricer.fastExpansionSum(bcd, abd, temp48a, bcdLen, abdLen);
	temp48bLen = arthemetricer.fastExpansionSum(cda, abc, temp48b, cdaLen, abcLen);
	for (int i = 0; i < temp48bLen; i++)
		temp48b[i] = -temp48b[i];
	xyzwLen = arthemetricer.fastExpansionSum(temp48a, temp48b, xyzw, temp48aLen, temp48bLen);
	xLen = arthemetricer.scaleExpansion(xyzw, e.x, temp192, xyzwLen);
	xLen = arthemetricer.scaleExpansion(temp192, e.x, temp384x, xLen);
	yLen = arthemetricer.scaleExpansion(xyzw, e.y, temp192, xyzwLen);
	yLen = arthemetricer.scaleExpansion(temp192, e.y, temp384y, yLen);
	zLen = arthemetricer.scaleExpansion(xyzw, e.z, temp192, xyzwLen);
	zLen = arthemetricer.scaleExpansion(temp192, e.z, temp384z, zLen);
	wLen = arthemetricer.scaleExpansion(xyzw, -eWeight, temp192w, xyzwLen);

	xyLen = arthemetricer.fastExpansionSum(temp384x, temp384y, temp768xy, xLen, yLen);
	xyzLen = arthemetricer.fastExpansionSum(temp768xy, temp384z, temp1152xyz, xyLen, zLen);
	edetLen = arthemetricer.fastExpansionSum(temp1152xyz, temp192w, edet, xyzLen, wLen);

	FT abdet[2688], cddet[2688], cdedet[4032], fin[6720];
	int abdetLen, cddetLen, cdedetLen, finLength;
	abdetLen = arthemetricer.fastExpansionSum(adet, bdet, abdet, adetLen, bdetLen);
	cddetLen = arthemetricer.fastExpansionSum(cdet, ddet, cddet, cdetLen, ddetLen);
	cdedetLen = arthemetricer.fastExpansionSum(cddet, edet, cdedet, cddetLen, edetLen);
	finLength = arthemetricer.fastExpansionSum(cdedet, abdet, fin, cdedetLen, abdetLen);

	return fin[finLength - 1];
}

template<class FT> bool Predicator<FT>::inSegmentRange(const VectorBase<FT>& u, const VectorBase<FT>& a, const VectorBase<FT>& b) const {
	VectorBase<FT> ab = b - a;
	VectorBase<FT> au = u - a;

	FT uabx = au.x * ab.x; FT uaby = au.y * ab.y; FT uabz = au.z * ab.z;
	FT uabLen = uabx + uaby + uabz;
	FT uabNorm = fabs(uabx) + fabs(uaby) + fabs(uabz);
	if (fabs(uabLen) >= inSegRangeErrorBoundA * uabNorm) {
		if (uabLen < FT(0))
			return false;
	}
	else {
		if (inSegmentRangeHalfAdaptive(u, a, b, ab, au, uabNorm) < FT(0))
			return false;
	}

	VectorBase<FT> ba = -ab;
	VectorBase<FT> bu = u - b;
	FT ubax = bu.x * ba.x; FT ubay = bu.y * ba.y; FT ubaz = bu.z * ba.z;
	FT ubaLen = ubax + ubay + ubaz;
	FT ubaNorm = fabs(ubax) + fabs(ubay) + fabs(ubaz);
	if (fabs(ubaLen) >= inSegRangeErrorBoundA * ubaNorm) {
		if (ubaLen < FT(0))
			return false;
	}
	else {
		if (inSegmentRangeHalfAdaptive(u, b, a, ba, bu, ubaNorm) < FT(0))
			return false;
	}

	return true;
}

template<class FT> FT Predicator<FT>::inSegmentRangeHalfAdaptive(const VectorBase<FT>& u, const VectorBase<FT>& a, const VectorBase<FT>& b,
	const VectorBase<FT>& ab, const VectorBase<FT>& au, FT norm) const {
	FT x, y, z, xErr, yErr, zErr;
	arthemetricer.twoProduct(ab.x, au.x, x, xErr);
	arthemetricer.twoProduct(ab.y, au.y, y, yErr);
	arthemetricer.twoProduct(ab.z, au.z, z, zErr);

	FT retEsti[6];
	arthemetricer.twoTwoTwoSum(x, xErr, y, yErr, z, zErr, retEsti[5], retEsti[4], retEsti[3], retEsti[2], retEsti[1], retEsti[0]);

	FT ret = arthemetricer.Estimate(retEsti, 6);

	if (fabs(ret) >= inSegRangeErrorBoundB * norm)
		return ret;

	FT auxErr = arthemetricer.twoDiffError(u.x, a.x, au.x);
	FT auyErr = arthemetricer.twoDiffError(u.y, a.y, au.y);
	FT auzErr = arthemetricer.twoDiffError(u.z, a.z, au.z);
	FT abxErr = arthemetricer.twoDiffError(b.x, a.x, ab.x);
	FT abyErr = arthemetricer.twoDiffError(b.y, a.y, ab.y);
	FT abzErr = arthemetricer.twoDiffError(b.z, a.z, ab.z);

	if (auxErr == FT(0) && auyErr == FT(0) && auzErr == FT(0) &&
		abxErr == FT(0) && abyErr == FT(0) && abzErr == FT(0))
		return ret;

	FT errorBound = commonBound * fabs(ret) + inSegRangeErrorBoundC * norm;
	ret += ((auxErr * ab.x + au.x * abxErr) + (auyErr * ab.y + au.y * abyErr) + (auzErr * ab.z + au.z * abzErr));

	if (fabs(ret) >= errorBound)
		return ret;
	
	//Exact
	FT s1, s0, t1, t0, u1, u0;
	FT v[4], w[6], c0[10], c1[14], c2[18], fin[24];
	arthemetricer.twoProduct(au.x, abxErr, s1, s0);
	arthemetricer.twoProduct(auxErr, ab.x, t1, t0);
	arthemetricer.twoTwoSum(s1, s0, t1, t0, v[3], v[2], v[1], v[0]);
	int c0Len = arthemetricer.fastExpansionSum(retEsti, v, c0, 6, 4);

	arthemetricer.twoProduct(au.y, abyErr, s1, s0);
	arthemetricer.twoProduct(auyErr, ab.y, t1, t0);
	arthemetricer.twoTwoSum(s1, s0, t1, t0, v[3], v[2], v[1], v[0]);
	int c1Len = arthemetricer.fastExpansionSum(c0, v, c1, c0Len, 4);

	arthemetricer.twoProduct(au.z, abzErr, s1, s0);
	arthemetricer.twoProduct(auzErr, ab.z, t1, t0);
	arthemetricer.twoTwoSum(s1, s0, t1, t0, v[3], v[2], v[1], v[0]);
	int c2Len = arthemetricer.fastExpansionSum(c1, v, c2, c1Len, 4);

	arthemetricer.twoProduct(auxErr, abxErr, s1, s0);
	arthemetricer.twoProduct(auyErr, abyErr, t1, t0);
	arthemetricer.twoProduct(auzErr, abzErr, u1, u0);
	arthemetricer.twoTwoTwoSum(s1, s0, t1, t0, u1, u0, w[5], w[4], w[3], w[2], w[1], w[0]);
	int finLen = arthemetricer.fastExpansionSum(c2, w, fin, c2Len, 6);

	return fin[finLen - 1];
}

template<class FT> FT Predicator<FT>::inOrthoCirclePerturbed(const VectorBase<FT> &a, FT aWeight, const VectorBase<FT> &b, FT bWeight,
	const VectorBase<FT> &c, FT cWeight, const VectorBase<FT> &d, FT dWeight, const VectorBase<FT> &above) const{
	FT test = inOrthoSphere(above, FT(0), a, aWeight, b, bWeight, c, cWeight, d, dWeight);

	if (test != FT(0))
		return test;

	int swaps = 0;
	int count = 0;
	VectorBase<FT> v[4];
	v[0] = a, v[1] = b, v[2] = c, v[3] = d;
	int num = 4;
	do {
		count = 0;
		num--;
		for (int i = 0; i < num; i++) {
			if (v[i] > v[i + 1]) {
				std::swap(v[i], v[i + 1]);
				count++;
			}
		}
		swaps += count;
	} while (count > 0);

	FT ori = orient3d(v[1], v[2], v[3], above);
	if ((swaps % 2) != 0) ori = -ori;
	return ori;
}

template<class FT> FT Predicator<FT>::inOrthoSpherePerturbed(const VectorBase<FT> &a, FT aWeight, const VectorBase<FT> &b, FT bWeight,
	const VectorBase<FT> &c, FT cWeight, const VectorBase<FT> &d, FT dWeight, const VectorBase<FT> &e, FT eWeight) const{
	FT test = inOrthoSphere(a, aWeight, b, bWeight, c, cWeight, d, dWeight, e, eWeight);

	if (test != 0.0)
		return test;

	int swaps = 0;
	int count = 0;
	VectorBase<FT> v[5];
	v[0] = a, v[1] = b, v[2] = c, v[3] = d, v[4] = e;
	int num = 5;
	do {
		count = 0;
		num--;
		for (int i = 0; i < num; i++) {
			if (v[i] > v[i + 1]) {
				std::swap(v[i], v[i + 1]);
				count++;
			}
		}
		swaps += count;
	} while (count > 0);

	FT oriA = orient3d(v[1], v[2], v[3], v[4]);
	if (oriA != 0.0) {
		if ((swaps % 2) != 0) oriA = -oriA;
		return oriA;
	}

	FT oriB = -orient3d(v[0], v[2], v[3], v[4]);
	if ((swaps % 2) != 0) oriB = -oriB;
	return oriB;
}

template<class FT> bool Predicator<FT>::inHalfSpace3D(const VectorBase<FT> &u, const VectorBase<FT> &a, const VectorBase<FT>& b, const VectorBase<FT> &c) const{
	FT orient = orient3d(u, a, b, c);
	if (orient > FT(0))
		return true;
	else if (orient == FT(0))
		return inSphere(a + Geometer::triangleNormal(a, b, c), a, b, c, u) > FT(0);
	else
		return false;
}

template<class FT> bool Predicator<FT>::inOrthoHalfSpace3D(const VectorBase<FT> &u, FT uWeight, const VectorBase<FT> &a, FT aWeight, const VectorBase<FT>& b, 
	        FT bWeight, const VectorBase<FT> &c, FT cWeight) const{
	FT orient = orient3d(u, a, b, c);
	if (orient > FT(0))
		return true;
	else if (orient == FT(0))
		return inOrthoCirclePerturbed(a, aWeight, b, bWeight, c, cWeight, u, uWeight, a + Geometer::triangleNormal(a, b, c)) > FT(0);
	else
		return false;
}

template<class FT> bool Predicator<FT>::inHalfSpace2D(const VectorBase<FT>& u, const VectorBase<FT>& b, const VectorBase<FT>& c, const VectorBase<FT>& above) const{
	FT orient = orient2d(u, b, c, above);
	if (orient > FT(0))
		return true;
	else if (orient == FT(0)){
		return inSegmentRange(u, b, c);
	}
	else
		return false;
}

template<class FT> bool Predicator<FT>::Intersection(const VectorBase<FT>& p, const VectorBase<FT>& q, const VectorBase<FT>& r,
	const VectorBase<FT>& a, const VectorBase<FT>& b, const VectorBase<FT>& c) const {

	FT dp = orient3d(p, a, b, c);
	FT dq = orient3d(q, a, b, c);
	FT dr = orient3d(r, a, b, c);

	if (dp != 0 || dq != 0 || dr != 0) {
		const VectorBase<FT> *sMin1 = NULL, *sMax1 = NULL, *tMin1 = NULL, *tMax1 = NULL;
		int condition = (dp >= 0) + ((dq >= 0) << 1) + ((dr >= 0) << 2);

		switch (condition) {
		case 0:
			return false;
		case 1: //positive p
			sMin1 = &p; tMin1 = &q; sMax1 = &r; tMax1 = &p;
			break;
		case 2: //positive q
			sMin1 = &q; tMin1 = &r; sMax1 = &p; tMax1 = &q;
			break;
		case 3: //negative r
			sMin1 = &q; tMin1 = &r; sMax1 = &r; tMax1 = &p;
			break;
		case 4: //positive r
			sMin1 = &r; tMin1 = &p; sMax1 = &q; tMax1 = &r;
			break;
		case 5: //negative q
			sMin1 = &p; tMin1 = &q; sMax1 = &q; tMax1 = &r;
			break;
		case 6: //negative p
			sMin1 = &r; tMin1 = &p; sMax1 = &p; tMax1 = &q;
			break;
		case 7:
		{
			VectorBase<FT> normal = Geometer::triangleNormal(a, b, c);
			FT nx = fabs(normal.x), ny = fabs(normal.y), nz = fabs(normal.z);
			Plane plane = Plane::Plane_Arbitary;
			if (nx != 0 || ny != 0 || nz != 0)
				plane = nz > ny ? (nz > nx ? Plane::Plane_XY : Plane::Plane_YZ) : (ny > nx ? Plane::Plane_XZ : Plane::Plane_YZ);

			int zeroCondition = (dp == 0) + ((dq == 0) << 1) + ((dr == 0) << 2);
			switch (zeroCondition) {
			case 0:
				return false;
			case 1:
				return orientCoplane(p, a, b, plane) >= 0 &&
					orientCoplane(p, b, c, plane) >= 0 &&
					orientCoplane(p, c, a, plane) >= 0;
			case 2:
				return orientCoplane(q, a, b, plane) >= 0 &&
					orientCoplane(q, b, c, plane) >= 0 &&
					orientCoplane(q, c, a, plane) >= 0;
			case 3:
				return intersectionTriSegCoplane(a, b, c, p, q, plane);
			case 4:
				return orientCoplane(r, a, b, plane) >= 0 &&
					orientCoplane(r, b, c, plane) >= 0 &&
					orientCoplane(r, c, a, plane) >= 0;
			case 5:
				return intersectionTriSegCoplane(a, b, c, r, p, plane);
			case 6:
				return intersectionTriSegCoplane(a, b, c, q, r, plane);
			default:
				Severe("Unexpected Case Triangle-Triangle Predicator::Intersection");
				break;
			}
			break;
		}
		default:
			Severe("Unexpected Case Triangle-Triangle Predicator::Intersection");
			break;
		}

		FT da = orient3d(a, p, q, r);
		FT db = orient3d(b, p, q, r);
		FT dc = orient3d(c, p, q, r);
		Assert(da != 0 || db != 0 || dc != 0);

		const VectorBase<FT> *sMin2 = NULL, *sMax2 = NULL, *tMin2 = NULL, *tMax2 = NULL;
		condition = (da >= 0) + ((db >= 0) << 1) + ((dc >= 0) << 2);

		switch (condition) {
		case 0:
			return false;
		case 1: //positive a
			sMin2 = &a; tMin2 = &b; sMax2 = &c; tMax2 = &a;
			break;
		case 2: //positive b
			sMin2 = &b; tMin2 = &c; sMax2 = &a; tMax2 = &b;
			break;
		case 3: //negative c
			sMin2 = &b; tMin2 = &c; sMax2 = &c; tMax2 = &a;
			break;
		case 4: //positive c
			sMin2 = &c; tMin2 = &a; sMax2 = &b; tMax2 = &c;
			break;
		case 5: //negative b
			sMin2 = &a; tMin2 = &b; sMax2 = &b; tMax2 = &c;
			break;
		case 6: //negative a
			sMin2 = &c; tMin2 = &a; sMax2 = &a; tMax2 = &b;
			break;
		case 7:
		{
			VectorBase<FT> normal = Geometer::triangleNormal(p, q, r);
			FT nx = fabs(normal.x), ny = fabs(normal.y), nz = fabs(normal.z);
			Plane plane = Plane::Plane_Arbitary;
			if (nx != 0 || ny != 0 || nz != 0)
				plane = nz > ny ? (nz > nx ? Plane::Plane_XY : Plane::Plane_YZ) : (ny > nx ? Plane::Plane_XZ : Plane::Plane_YZ);

			int zeroCondition = (da == 0) + ((db == 0) << 1) + ((dc == 0) << 2);
			switch (zeroCondition) {
			case 0:
				return false;
			case 1:
				return orientCoplane(a, p, q, plane) >= 0 &&
					orientCoplane(a, q, r, plane) >= 0 &&
					orientCoplane(a, r, p, plane) >= 0;
			case 2:
				return orientCoplane(b, p, q, plane) >= 0 &&
					orientCoplane(b, q, r, plane) >= 0 &&
					orientCoplane(b, r, p, plane) >= 0;
			case 3:
				return intersectionTriSegCoplane(p, q, r, a, b, plane);
			case 4:
				return orientCoplane(c, p, q, plane) >= 0 &&
					orientCoplane(c, q, r, plane) >= 0 &&
					orientCoplane(c, r, p, plane) >= 0;
			case 5:
				return intersectionTriSegCoplane(p, q, r, c, a, plane);
			case 6: 
				return intersectionTriSegCoplane(p, q, r, b, c, plane);
			default:
				Severe("Unexpected Case Triangle-Triangle Predicator::Intersection");
				break;
			}
			break;
		}
		default:
			Severe("Unexpected Case Triangle-Triangle Predicator::Intersection");
			break;
		}

		return orient3d(*sMin1, *tMin1, *tMin2, *sMin2) <= 0 &&
			orient3d(*sMax1, *tMax1, *sMax2, *tMax2) <= 0;
	}
	else {
		VectorBase<FT> normal = Geometer::triangleNormal(p, q, r);
		bool zeroTest = (normal.x == 0 && normal.y == 0 && normal.z == 0);
		if (zeroTest) {
			normal = Geometer::triangleNormal(a, b, c);
			zeroTest = (normal.x == 0 && normal.y == 0 && normal.z == 0);
		}

		Plane plane = Plane::Plane_Arbitary;
		if (!zeroTest) {
			FT nx = fabs(normal.x), ny = fabs(normal.y), nz = fabs(normal.z);
			plane = nz > ny ? (nz > nx ? Plane::Plane_XY : Plane::Plane_YZ) : (ny > nx ? Plane::Plane_XZ : Plane::Plane_YZ);
		}

		if (orientCoplane(p, b, a, plane) > 0) {
			if (orientCoplane(p, c, a, plane) > 0) {
				if (orientCoplane(p, b, c, plane) > 0)
					return intersectionTestEdge(p, q, r, a, b, plane);
				
				return intersectionTestVertex(p, q, r, b, c, a, plane);
			}
			return intersectionTestVertex(p, q, r, a, b, c, plane);
		}

		if (orientCoplane(p, c, b, plane) > 0) {
			if (orientCoplane(p, c, a, plane) > 0)
				return intersectionTestEdge(p, q, r, b, c, plane);

			return intersectionTestVertex(p, q, r, c, a, b, plane);
		}

		if (orientCoplane(p, a, c, plane) > 0)
			return intersectionTestEdge(p, q, r, c, a, plane);

		return true;
	}
}

template<class FT> bool Predicator<FT>::intersectionTestEdge(const VectorBase<FT>& p, const VectorBase<FT>& q, const VectorBase<FT>& r,
	const VectorBase<FT>& a, const VectorBase<FT>& b, Plane hint) const {

	if (orientCoplane(q, a, b, hint) >= 0) {
		if (orientCoplane(p, b, q, hint) >= 0)
			return orientCoplane(p, q, a, hint) >= 0;

		return orientCoplane(p, b, r, hint) >= 0 && orientCoplane(q, r, b, hint) >= 0;
	}

	if (orientCoplane(r, a, b, hint) >= 0)
		return orientCoplane(q, r, a, hint) >= 0 && orientCoplane(p, b, r, hint) >= 0;

	return false;
}

template<class FT> bool Predicator<FT>::intersectionTestVertex(const VectorBase<FT>& p, const VectorBase<FT>& q, const VectorBase<FT>& r,
	const VectorBase<FT>& a, const VectorBase<FT>& b, const VectorBase<FT>& c, Plane hint) const {

	if (orientCoplane(q, a, b, hint) >= 0) {
		if (orientCoplane(c, a, q, hint) >= 0) {
			if (orientCoplane(p, b, q, hint) >= 0)
				return orientCoplane(p, c, q, hint) <= 0;

			return orientCoplane(p, b, r, hint) >= 0 && orientCoplane(q, r, b, hint) >= 0;
		}

		return orientCoplane(p, c, q, hint) <= 0 && orientCoplane(c, a, r, hint) >= 0 &&
			orientCoplane(q, r, c, hint) >= 0;
	}

	if (orientCoplane(r, a, b, hint) >= 0) {
		if (orientCoplane(q, r, a, hint) >= 0)
			return orientCoplane(r, p, b, hint) >= 0;

		return orientCoplane(q, r, c, hint) >= 0 && orientCoplane(c, a, r, hint) >= 0;
	}

	return false;
}

template<class FT> bool Predicator<FT>::Intersection(const VectorBase<FT>& a, const VectorBase<FT>& b, const VectorBase<FT>& c,
	const VectorBase<FT>& p, const VectorBase<FT>& q) const {

	FT dp = orient3d(p, a, b, c);
	FT dq = orient3d(q, a, b, c);

	if (dp != 0 || dq != 0) {
		int condition = (dp >= 0) + ((dq >= 0) << 1);

		switch (condition) {
		case 0:
			return false;
		case 1:
			return orient3d(p, q, a, b) >= 0 && orient3d(p, q, b, c) >= 0 && orient3d(p, q, c, a) >= 0;
		case 2:
			return orient3d(q, p, a, b) >= 0 && orient3d(q, p, b, c) >= 0 && orient3d(q, p, c, a) >= 0;
		case 3:
		{
			if(dp == 0) return orient3d(q, p, a, b) >= 0 && orient3d(q, p, b, c) >= 0 && orient3d(q, p, c, a) >= 0;
			if(dq == 0) return orient3d(p, q, a, b) >= 0 && orient3d(p, q, b, c) >= 0 && orient3d(p, q, c, a) >= 0;
			return false;
		}
		default:
			Severe("Unexpected Case Triangle-Segment Predicator::Intersection");
			break;
		}
		return false;
	}
	else {
		VectorBase<FT> normal = Geometer::triangleNormal(a, b, c);
		FT nx = fabs(normal.x), ny = fabs(normal.y), nz = fabs(normal.z);
		Plane plane = Plane::Plane_Arbitary;
		if(nx != 0 || ny != 0 || nz != 0)
			plane = nz > ny ? (nz > nx ? Plane::Plane_XY : Plane::Plane_YZ) : (ny > nx ? Plane::Plane_XZ : Plane::Plane_YZ);

		return intersectionTriSegCoplane(a, b, c, p, q, plane);
	}
}

template<class FT> bool Predicator<FT>::intersectionTriSegCoplane(const VectorBase<FT>& a, const VectorBase<FT>& b, const VectorBase<FT>& c,
	const VectorBase<FT>& p, const VectorBase<FT>& q, Plane hint) const {
	FT pqa = orientCoplane(p, q, a, hint);
	FT pqb = orientCoplane(p, q, b, hint);
	FT pqc = orientCoplane(p, q, c, hint);

	int condition = (pqa >= 0) + ((pqb >= 0) << 1) + ((pqc >= 0) << 2);

	switch (condition) {
	case 0:
		return false;
	case 1:
		return orientCoplane(p, c, a, hint) >= 0 && orientCoplane(q, a, b, hint) >= 0;
	case 2:
		return orientCoplane(p, a, b, hint) >= 0 && orientCoplane(q, b, c, hint) >= 0;
	case 3:
		return orientCoplane(p, c, a, hint) >= 0 && orientCoplane(q, b, c, hint) >= 0;
	case 4:
		return orientCoplane(p, b, c, hint) >= 0 && orientCoplane(q, c, a, hint) >= 0;
	case 5:
		return orientCoplane(p, b, c, hint) >= 0 && orientCoplane(q, a, b, hint) >= 0;
	case 6:
		return orientCoplane(p, a, b, hint) >= 0 && orientCoplane(q, c, a, hint) >= 0;
	case 7:
	{
		int zeroCondition = (pqa == 0) + ((pqb == 0) << 1) + ((pqc == 0) << 2);
		switch (zeroCondition) {
		case 0:
			return false;
		case 1:
			return orientCoplane(p, a, b, hint) >= 0 && orientCoplane(q, c, a, hint) >= 0;
		case 2:
			return orientCoplane(p, b, c, hint) >= 0 && orientCoplane(q, a, b, hint) >= 0;
		case 3:
			return orientCoplane(p, b, c, hint) >= 0 && orientCoplane(q, c, a, hint) >= 0;
		case 4:
			return orientCoplane(p, c, a, hint) >= 0 && orientCoplane(q, b, c, hint) >= 0;
		case 5:
			return orientCoplane(p, a, b, hint) >= 0 && orientCoplane(q, b, c, hint) >= 0;
		case 6:
			return orientCoplane(p, c, a, hint) >= 0 && orientCoplane(q, a, b, hint) >= 0;
		default:
			Severe("Unexpected Case Predicator::intersectionTriSegCoplane");
			break;
		}
		break;
	}
	default:
		Severe("Unexpected Case Predicator::intersectionTriSegCoplane");
		break;
	}
	return false; //never here
}

}

#if defined(_MSC_VER)
#pragma float_control(pop)
#pragma fp_contract(on)
#endif

#endif