#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_CORE_NUMERMETHOD_H
#define ODER_CORE_NUMERMETHOD_H

#include <math.h>
#include "oder.h"
#include "latool.h"
#include "sparseMatrix.h"

namespace ODER{
	template<class FT> void gaussianElimination3x3(FT *A, FT *rhs, FT *result){
		//a.k.a LU decomposition
		int rowPermute[3] = { 0, 1, 2 };
		int columnPermute[3] = { 0, 1, 2 };

		constexpr int order = 3;
		constexpr int count = order - 1;
		//forward substitution
		for (int i = 0; i < count; i++){
			FT largest = 0.0;
			int rowIndex = i;
			int columnIndex = i;
			for (int j = i; j < order; j++){
				int subRow = rowPermute[j];
				for (int k = i; k < order; k++){
					int subColumn = columnPermute[k];
					FT entry = fabs(A[order*subRow + subColumn]);
					if (entry > largest){
						largest = entry;
						rowIndex = j;
						columnIndex = k;
					}
				}
			}
			Assert(largest > 0.0);
			if (rowPermute[i] != rowPermute[rowIndex])
				std::swap(rowPermute[i], rowPermute[rowIndex]);
			if (columnPermute[i] != columnPermute[columnIndex])
				std::swap(columnPermute[i], columnPermute[columnIndex]);

			int row = rowPermute[i];
			int column = columnPermute[i];
			FT scale = FT(1.0) / A[order*row + column];
			for (int j = i + 1; j < order; j++){
				int subRow = rowPermute[j];
				FT lower = scale*A[order*subRow + column];
				for (int k = i + 1; k < order; k++){
					int subColumn = columnPermute[k];
					A[order*subRow + subColumn] -= lower * A[order*row + subColumn];
				}
				rhs[subRow] -= lower * rhs[row];
			}
		}

		//backward substitution
		for (int i = count; i >= 0; i--){
			int row = rowPermute[i];
			int column = columnPermute[i];
			FT dot = 0.0;
			for (int j = i + 1; j < order; j++){
				int subColumn = columnPermute[j];
				dot += A[order*row + subColumn] * result[subColumn];
			}
			result[column] = (rhs[row] - dot) / A[order*row + column];
		}
	}

	template<class FUN, class FT> FT findRoot(FUN f, FT start, FT end, FT tolerance = 0.0, FT epsilon = FT(2e-8)){
		FT a = start, b = end, c = start, d = 0.0, e = 0.0;
		FT fa = f(a), fb = f(b), fc = fa;

		constexpr int maxIter = 1000;
		for (int i = 0; i < maxIter; i++){
			if (fb > 0.0 && fc > 0.0 || fb < 0.0 && fc < 0.0){
				c = a; fc = fa;
				d = e = b - a;
			}
			if (fabs(fc) < fabs(fb)){
				a = b; b = c; c = a;
				fa = fb; fb = fc; fc = fa;
			}

			FT tol = FT(2.0)*FT(epsilon)*fabs(b) + tolerance;
			FT m = FT(0.5)*(c - b);

			if (fabs(m) <= tol || fb == FT(0.0)) return b;
			if (fabs(e) >= tol && fa > fb){
				//inverse quadratic
				FT s = fa / fb;
				FT p = 0.0, q = 0.0;
				if (a == c){
					p = FT(2.0)*m*s;
					q = FT(1.0) - s;
				}
				else{
					FT t = fa / fc;
					FT r = fb / fc;
					p = s * (FT(2.0)*t*(r - t)*m - (FT(1.0) - r)*(b - a));
					q = (t - FT(1.0))*(r - FT(1.0))*(s - FT(1.0));
				}
				if (p > 0.0) q = -q;
				else p = -p;

				if (FT(2.0)*p < std::min(FT(3.0)*q*m - fabs(tol*q), fabs(e*q))){
					e = d;
					d = p / q;
				}
				else{
					//bisection
					d = m;
					e = d;
				}
			}
			else{
				//bisection
				d = m;
				e = d;
			}
			a = b;
			fa = fb;
			if (fabs(d) > tol)
				b += d;
			else
				b += (m > 0 ? tol : -tol);
			fb = f(b);
		}
		Severe("max iteration in ODER::findRoot");
		return 0.0;
	}

	template<int blockLength, int blockWidth>
	void SpMDV(const BlockedSymSparseMatrix<blockLength, blockWidth>& mat, 
		const double *src, double *dest){
		Assert(src != dest);

		const int blockColumnCount = mat.numBlockColumn;
		const int columnCount = mat.numColumns;
		const int remainedColumnCount = mat.numRemainedColumn;

		const double *values = mat.values;
		const int *blockPcol = mat.blockPcol;
		const int *blockRows = mat.blockRows;

		constexpr int regularSize = blockLength * blockWidth;
		constexpr int diagSize = (blockLength + 1) * blockLength / 2;

		double destBlock[blockLength];
		double srcBlock[blockLength];

		const int blockBoundary = blockLength * blockColumnCount;
		for (int startColumn = 0; startColumn < blockBoundary; startColumn += blockLength){
			int index = *blockPcol++;
			int end = *blockPcol;

			if (index == end) continue;

			int row = blockRows[index];

			if (row == startColumn){
				int diagIndex = 0;
				for (int j = 0; j < blockLength; j++){
					int actualRow = row + j;
					double x = src[actualRow];
					double ret = values[diagIndex++] * x;
					for (int k = j + 1; k < blockLength; k++){
						double val = values[diagIndex++];
						dest[row + k] += val * x;
						ret += val * src[row + k];
					}
					dest[actualRow] += ret;
				}
				index++;
				values += diagSize;
			}

			if (index < end){
				for (int j = 0; j < blockLength; j++){
					srcBlock[j] = src[startColumn + j];
					destBlock[j] = 0.0;
				}
			}
			else continue;

			while (index < end - 1){
				row = blockRows[index++];
				for (int j = 0, offset = 0; j < blockLength; j++, offset += blockWidth){
					double x = srcBlock[j];
					double ret = 0.0;
					for (int k = 0; k < blockWidth; k++){
						int actualRow = row + k;
						double val = values[offset + k];
						dest[actualRow] += val * x;
						ret += val * src[actualRow];
					}
					destBlock[j] += ret;
				}
				values += regularSize;
			}

			row = blockRows[index];
			int mayDegenWidth = columnCount - row;
			if (mayDegenWidth >= blockWidth){
				for (int j = 0, offset = 0; j < blockLength; j++, offset += blockWidth){
					double x = srcBlock[j];
					double ret = 0.0;
					for (int k = 0; k < blockWidth; k++){
						int actualRow = row + k;
						double val = values[offset + k];
						dest[actualRow] += val * x;
						ret += val * src[actualRow];
					}
					destBlock[j] += ret;
				}
				values += regularSize;
			}
			else{
				for (int j = 0, offset = 0; j < blockLength; j++, offset += mayDegenWidth){
					double x = srcBlock[j];
					double ret = 0.0;
					for (int k = 0; k < mayDegenWidth; k++){
						int actualRow = row + k;
						double val = values[offset + k];
						dest[actualRow] += val * x;
						ret += val * src[actualRow];
					}
					destBlock[j] += ret;
				}
				values += blockLength*mayDegenWidth;
			}

			for (int j = 0; j < blockLength; j++)
				dest[startColumn + j] += destBlock[j];
		}

		for (int i = blockBoundary; i < columnCount; i++){
			int index = *blockPcol++;
			int end = *blockPcol;

			if (index == end) continue;

			int row = blockRows[index];
			double y = 0.0, x = src[i];
			if (row == i){
				y += x * (*values++);
				index++;
			}

			while (index < end){
				row = blockRows[index++];
				double val = *values++;

				y += val * src[row];
				dest[row] += val * x;
			}

			dest[i] += y;

		}
	}

	template<int blockLength, int blockWidth>
	void SpMSV(const BlockedSymSparseMatrix<blockLength, blockWidth>& mat,
		const std::vector<std::vector<std::pair<int, int>>>& fullIndices,
		const SparseVector& src, double *dest) {

		const double *values = mat.values;

		const auto end = src.cend();
		for (auto vecIter = src.cbegin(); vecIter != end; ++vecIter) {
			const int column = *vecIter.indexIterator;
			const double entry = *vecIter.valueIterator;

			const std::vector<std::pair<int, int>>& rowIndexPairs = fullIndices[column];

			for (const auto& pair : rowIndexPairs)
				dest[pair.first] += values[pair.second] * entry;
		}
	}

	template<int blockLength, int blockWidth>
	void SpMSV(const BlockedSymSparseMatrix<blockLength, blockWidth>& mat,
		const std::vector<std::vector<std::pair<int, int>>>& fullIndices,
		const SparseVector& src, FastSparseVector& dest){

		const double *values = mat.values;

		const auto end = src.cend();
		for (auto vecIter = src.cbegin(); vecIter != end; ++vecIter) {
			const int column = *vecIter.indexIterator;
			const double entry = *vecIter.valueIterator;

			const std::vector<std::pair<int, int>>& rowIndexPairs = fullIndices[column];

			for (const auto& pair : rowIndexPairs) 
				dest.Add(pair.first, values[pair.second] * entry);
		}
	}

	template<class FT> void eigenSym3x3(const FT *triMat, FT *eigenvalues, FT *eigenvectors) {
		static_assert(std::is_same<FT, float>::value || std::is_same<FT, double>::value, "ODER::eigenSym3x3 support IEEE 754-1985 floating point only");
		constexpr FT pi = M_PI;
		const FT aveTrace = (triMat[0] + triMat[3] + triMat[5]) / FT(3.0);
		const FT deviatoricDiag0 = triMat[0] - aveTrace;
		const FT deviatoricDiag1 = triMat[3] - aveTrace;
		const FT deviatoricDiag2 = triMat[5] - aveTrace;
		const FT invariant = (deviatoricDiag0  * deviatoricDiag0 + deviatoricDiag1  * deviatoricDiag1 + deviatoricDiag2  * deviatoricDiag2) * FT(0.5) +
			triMat[1] * triMat[1] + triMat[2] * triMat[2] + triMat[4] * triMat[4];
		const FT det = deviatoricDiag0 * (deviatoricDiag1 *  deviatoricDiag2 - triMat[4] * triMat[4]) +
			triMat[1] * (triMat[4] * triMat[2] - triMat[1] * deviatoricDiag2) +
			triMat[2] * (triMat[1] * triMat[4] - deviatoricDiag1 * triMat[2]);

		constexpr FT squareRoot3 = FT(1.7320508075688772935274463415058723669428052538103806280558069794519330169088);
		const FT sqrtTerm = sqrt(invariant);
		const FT angleCos = (FT(3) * FT(0.5) * squareRoot3) * ((det / invariant) / sqrtTerm);

		if (isnan(angleCos) || isinf(angleCos)) {
			eigenvalues[0] = eigenvalues[1] = eigenvalues[2] = aveTrace;
			Initiation(eigenvectors, 9);
			eigenvectors[0] = eigenvectors[4] = eigenvectors[8] = FT(1);
			return;
		}

		FT angle = acos(Clamp(angleCos, FT(-1), FT(1))) / FT(3);
		FT distinctEigenVal = (FT(2) / squareRoot3) * sqrtTerm * (angle < pi / FT(6) ? cos(angle) : cos(angle + (FT(2) / FT(3)) * pi));
		eigenvalues[0] = distinctEigenVal + aveTrace;

		const VectorBase<FT> r0(deviatoricDiag0 - distinctEigenVal, triMat[1], triMat[2]);
		const VectorBase<FT> r1(triMat[1], deviatoricDiag1 - distinctEigenVal, triMat[4]);
		const VectorBase<FT> r2(triMat[2], triMat[4], deviatoricDiag2 - distinctEigenVal);

		const FT r0Len = r0.length();
		const FT r1Len = r1.length();
		const FT r2Len = r2.length();

		VectorBase<FT> s0;
		VectorBase<FT> t1, t2;

		if (r0Len > r1Len) {
			if (r0Len > r2Len) {
				s0 = r0 / r0Len;
				t1 = r1 - (s0 * r1) * s0;
				t2 = r2 - (s0 * r2) * s0;
			}
			else {
				s0 = r2 / r2Len;
				t1 = r0 - (s0 * r0) * s0;
				t2 = r1 - (s0 * r1) * s0;
			}
		}
		else {
			if (r1Len > r2Len) {
				s0 = r1 / r1Len;
				t1 = r0 - (s0 * r0) * s0;
				t2 = r2 - (s0 * r2) * s0;
			}
			else {
				s0 = r2 / r2Len;
				t1 = r0 - (s0 * r0) * s0;
				t2 = r1 - (s0 * r1) * s0;
			}
		}


		const FT t1Len2 = t1.length2();
		const FT t2Len2 = t2.length2();

		const VectorBase<FT> s1 = t1Len2 > t2Len2 ? t1 / sqrt(t1Len2) : t2 / sqrt(t2Len2);
		VectorBase<FT> eigenVector0 = s0 % s1;
		memcpy(eigenvectors, &eigenVector0[0], sizeof(FT) * 3);

		Tensor2<FT> mat(deviatoricDiag0, triMat[1], triMat[2],
			triMat[1], deviatoricDiag1, triMat[4],
			triMat[2], triMat[4], deviatoricDiag2);

		const VectorBase<FT> trans0 = mat * s0;
		const VectorBase<FT> trans1 = mat * s1;
		const FT a11 = s0 * trans0;
		const FT a12 = s0 * trans1;
		const FT a22 = s1 * trans1; 

		const FT discriminant = sqrt((a11 - a22) * (a11 - a22) + FT(4) * a12 * a12);
		const FT eigenVal1 = ((a11 + a22) - (a11 > a22 ? discriminant : -discriminant)) / FT(2);
		const FT eigenVal2 = a11 + a22 - eigenVal1;
		eigenvalues[1] = eigenVal1 + aveTrace;
		eigenvalues[2] = eigenVal2 + aveTrace;

		VectorBase<FT> eigenVector1;
		const VectorBase<FT> u0 = (trans0 - eigenVal1 * s0) % eigenVector0;
		const VectorBase<FT> u1 = (trans1 - eigenVal1 * s1) % eigenVector0;

		const FT u0Len2 = u0.length2();
		const FT u1Len2 = u1.length2();

		if (u0Len2 != FT(0) && u0Len2 != FT(0)) eigenVector1 = u0Len2 > u1Len2 ? u0 / sqrt(u0Len2) : u1 / sqrt(u1Len2);
		else eigenVector1 = s0;

		memcpy(eigenvectors + 3, &eigenVector1[0], sizeof(FT) * 3);

		VectorBase<FT> eigenVector2 = eigenVector0 % eigenVector1;
		memcpy(eigenvectors + 6, &eigenVector2[0], sizeof(FT) * 3);
	}

	template<class FT> void properQRDecompostion3x3(const FT *mat, FT* ortho, FT *upperTri) {
		auto givens = [](FT a, FT b, FT& c, FT& s, FT& r){
			if (b == FT(0)) {
				c = a < FT(0) ? FT(-1) : FT(1);
				s = FT(0);
				r = fabs(a);
			}
			else if (a == FT(0)) {
				c = FT(0);
				s = b < FT(0) ? FT(1) : FT(-1);
				r = fabs(b);
			}
			else if (fabs(b) > fabs(a)) {
				FT tau = -a / b;
				FT u = sqrt(FT(1) + tau * tau);
				if (b < FT(0)) u = -u;
				s = -FT(1) / u;
				c = s * tau;
				r = b * u;
			}
			else {
				FT tau = -b / a;
				FT u = sqrt(FT(1) + tau * tau);
				if (a < FT(0)) u = -u;
				c = FT(1) / u;
				s = c * tau;
				r = a * u;
			}
		};

		FT Q[9], R[9];
		memcpy(R, mat, sizeof(FT) * 9);
		Initiation(Q, 9);
		//three givens transformation
		FT c = FT(0), s = FT(0), r = FT(0);
		givens(R[3], R[6], c, s, r);
		Q[4] = Q[8] = c; Q[5] = s; Q[7] = -s;
		R[3] = r; R[6] = FT(0);
		for (int i = 0; i < 2; i++) {
			FT p = c * R[4 + i] - s * R[7 + i];
			FT q = s * R[4 + i] + c * R[7 + i];
			R[4 + i] = p;
			R[7 + i] = q;
		}

		givens(R[0], R[3], c, s, r);
		//update R
		R[0] = r; R[3] = FT(0);
		for (int i = 0; i < 2; i++) {
			FT p = c * R[1 + i] - s * R[4 + i];
			FT q = s * R[1 + i] + c * R[4 + i];
			R[1 + i] = p;
			R[4 + i] = q;
		}
		//update Q
		Q[0] = c; Q[1] = s;
		Q[3] = Q[4] * -s; Q[4] *= c;
		Q[6] = Q[7] * -s; Q[7] *= c;

		givens(R[4], R[7], c, s, r);
		//update R
		R[4] = r; R[7] = FT(0);
		FT pp = c * R[5] - s * R[8];
		FT qq = s * R[5] + c * R[8];
		R[5] = pp; R[8] = qq;
		//update Q
		for (int i = 0; i < 3; i++) {
			FT p = c * Q[i * 3 + 1] - s * Q[i * 3 + 2];
			FT q = s * Q[i * 3 + 1] + c * Q[i * 3 + 2];
			Q[i * 3 + 1] = p;
			Q[i * 3 + 2] = q;
		}

		memcpy(ortho, Q, sizeof(FT) * 9);
		memcpy(upperTri, R, sizeof(FT) * 9);
	}

	template<class FT> FT Dot(int width, const FT *x, const FT *y) {
		FT result = FT(0);
		for (int i = 0; i < width; i++)
			result += x[i] * y[i];
		return result;
	}
}
#endif