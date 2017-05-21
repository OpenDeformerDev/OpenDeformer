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
			FT largest = FT(0);
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
			Assert(largest > FT(0));
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
			FT dot = FT(0);
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
				b += (m > 0.0 ? tol : -tol);
			fb = f(b);
		}
		Severe("max iteration in ODER::findRoot");
		return Scalar(0);
	}

	template<int blockLength, int blockWidth>
	void SpMDV(const BlockedSymSparseMatrix<blockLength, blockWidth>& mat, 
		const Scalar *src, Scalar *dest){
		Assert(src != dest);

		const int blockColumnCount = mat.numBlockColumn;
		const int columnCount = mat.numColumns;
		const int remainedColumnCount = mat.numRemainedColumn;

		const Scalar *values = mat.values;
		const int *blockPcol = mat.blockPcol;
		const int *blockRows = mat.blockRows;

		constexpr int regularSize = blockLength * blockWidth;
		constexpr int diagSize = (blockLength + 1) * blockLength / 2;

		Scalar destBlock[blockLength];
		Scalar srcBlock[blockLength];

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
					Scalar x = src[actualRow];
					Scalar ret = values[diagIndex++] * x;
					for (int k = j + 1; k < blockLength; k++){
						Scalar val = values[diagIndex++];
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
					destBlock[j] = Scalar(0);
				}
			}
			else continue;

			while (index < end - 1){
				row = blockRows[index++];
				for (int j = 0, offset = 0; j < blockLength; j++, offset += blockWidth){
					Scalar x = srcBlock[j];
					Scalar ret = Scalar(0);
					for (int k = 0; k < blockWidth; k++){
						int actualRow = row + k;
						Scalar val = values[offset + k];
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
					Scalar x = srcBlock[j];
					Scalar ret = Scalar(0);
					for (int k = 0; k < blockWidth; k++){
						int actualRow = row + k;
						Scalar val = values[offset + k];
						dest[actualRow] += val * x;
						ret += val * src[actualRow];
					}
					destBlock[j] += ret;
				}
				values += regularSize;
			}
			else{
				for (int j = 0, offset = 0; j < blockLength; j++, offset += mayDegenWidth){
					Scalar x = srcBlock[j];
					Scalar ret = Scalar(0);
					for (int k = 0; k < mayDegenWidth; k++){
						int actualRow = row + k;
						Scalar val = values[offset + k];
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
			Scalar y = Scalar(0), x = src[i];
			if (row == i){
				y += x * (*values++);
				index++;
			}

			while (index < end){
				row = blockRows[index++];
				Scalar val = *values++;

				y += val * src[row];
				dest[row] += val * x;
			}

			dest[i] += y;

		}
	}

	template<int blockLength, int blockWidth>
	void SpMSV(const BlockedSymSparseMatrix<blockLength, blockWidth>& mat,
		const std::vector<std::vector<std::pair<int, int>>>& fullIndices,
		const SparseVector& src, Scalar *dest) {

		const Scalar *values = mat.values;

		const auto end = src.cend();
		for (auto vecIter = src.cbegin(); vecIter != end; ++vecIter) {
			const int column = *vecIter.indexIterator;
			const Scalar entry = *vecIter.valueIterator;

			const std::vector<std::pair<int, int>>& rowIndexPairs = fullIndices[column];

			for (const auto& pair : rowIndexPairs)
				dest[pair.first] += values[pair.second] * entry;
		}
	}

	template<int blockLength, int blockWidth>
	void SpMSV(const BlockedSymSparseMatrix<blockLength, blockWidth>& mat,
		const std::vector<std::vector<std::pair<int, int>>>& fullIndices,
		const SparseVector& src, FastSparseVector& dest){

		const Scalar *values = mat.values;

		const auto end = src.cend();
		for (auto vecIter = src.cbegin(); vecIter != end; ++vecIter) {
			const int column = *vecIter.indexIterator;
			const Scalar entry = *vecIter.valueIterator;

			const std::vector<std::pair<int, int>>& rowIndexPairs = fullIndices[column];

			for (const auto& pair : rowIndexPairs) 
				dest.Add(pair.first, values[pair.second] * entry);
		}
	}

	template<class FT> void eigenSym3x3(const FT *triMat, FT *eigenvalues, FT *eigenvectors) {
		static_assert(std::is_same<FT, float>::value || std::is_same<FT, double>::value, "ODER::eigenSym3x3 support IEEE 754-1985 floating point only");
		constexpr FT pi = FT(M_PI);
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
		static_assert(std::is_same<FT, float>::value || std::is_same<FT, double>::value, "ODER::properQRDecompostion3x3 support IEEE 754-1985 floating point only");
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

	namespace PolarDecompostion3x3Internal {
		template<class FT> void Orthonormalize(FT mat[8]) {
			constexpr FT tol = std::numeric_limits<FT>::min();
			FT beta0 = FT(0), beta1 = FT(0);
			FT v0[3], v1[2];

			//first householder
			v0[0] = mat[1]; v0[1] = mat[2]; v0[2] = mat[3];
			FT norm2 = v0[0] * v0[0] + v0[1] * v0[1] + v0[2] * v0[2];
			if (norm2 <= tol) {
				v0[0] = v0[1] = v0[2] = FT(0);
			}
			else {
				FT c0 = mat[0];
				FT sigma = sqrt(c0 * c0 + norm2);
				if (c0 >= FT(0)) sigma = -sigma;

				FT factor = c0 - sigma;
				v0[0] /= factor; v0[1] /= factor; v0[2] /= factor;
				beta0 = (sigma - c0) / sigma;
			}

			//apply first householder
			FT product = mat[4] + v0[0] * mat[5] + v0[1] * mat[6] + v0[2] * mat[7];
			FT c0 = mat[5] - product * beta0 * v0[0];
			v1[0] = mat[6] - product * beta0 * v0[1];
			v1[1] = mat[7] - product * beta0 * v0[2];

			//second householder
			norm2 = v1[0] * v1[0] + v1[1] * v1[1];
			if (norm2 <= tol) {
				v1[0] = v1[1] = FT(0);
			}
			else {
				FT sigma = sqrt(c0 * c0 + norm2);
				if (c0 >= FT(0)) sigma = -sigma;

				FT factor = c0 - sigma;
				v1[0] /= factor; v1[1] /= factor;
				beta1 = (sigma - c0) / sigma;
			}

			mat[0] = FT(1) - beta0; mat[1] = -beta0 * v0[0]; mat[2] = -beta0 * v0[1]; mat[3] = -beta0 * v0[2];
			//the second column
			FT subMat[9], lhs[3];
			subMat[0] = FT(1) - beta0 * v0[0] * v0[0];
			subMat[1] = -beta0 * v0[0] * v0[1]; subMat[3] = subMat[1];
			subMat[2] = -beta0 * v0[0] * v0[2]; subMat[6] = subMat[2];
			subMat[4] = FT(1) - beta0 * v0[1] * v0[1];
			subMat[5] = -beta0 * v0[1] * v0[2]; subMat[7] = subMat[5];
			subMat[8] = FT(1) - beta0 * v0[2] * v0[2];

			lhs[0] = FT(1) - beta1; lhs[1] = -beta1 * v1[0]; lhs[2] = -beta1 * v1[1];

			mat[4] = mat[1] * lhs[0] + mat[2] * lhs[1] + mat[3] * lhs[2];
			mat[5] = subMat[0] * lhs[0] + subMat[1] * lhs[1] + subMat[2] * lhs[2];
			mat[6] = subMat[3] * lhs[0] + subMat[4] * lhs[1] + subMat[5] * lhs[2];
			mat[7] = subMat[6] * lhs[0] + subMat[7] * lhs[1] + subMat[8] * lhs[2];
		}
	}

	template<class FT> void polarDecompostion3x3(FT *mat, FT *ortho, FT *sspd) {
		static_assert(std::is_same<FT, float>::value || std::is_same<FT, double>::value, "ODER::polarDecompostion3x3 support IEEE 754-1985 floating point only");
		using namespace PolarDecompostion3x3Internal;
		FT A[9], B[16];
		FT fnorm = FT(0);
		for (int i = 0; i < 9; i++) fnorm += mat[i] * mat[i];
		fnorm = sqrt(fnorm);
		if (fnorm == FT(0)) {
			Initiation(ortho, 9);
			Initiation(sspd, 9);
			ortho[0] = ortho[4] = ortho[8] = FT(1);
			return;
		}
		for (int i = 0; i < 9; i++) A[i] = mat[i] / fnorm;

		FT detB = FT(0), adjAEntry = FT(0);
		adjAEntry = A[4] * A[8] - A[5] * A[7]; detB += adjAEntry * adjAEntry;
		adjAEntry = A[3] * A[8] - A[5] * A[6]; detB += adjAEntry * adjAEntry;
		adjAEntry = A[3] * A[7] - A[4] * A[6]; detB += adjAEntry * adjAEntry;
		adjAEntry = A[1] * A[8] - A[2] * A[7]; detB += adjAEntry * adjAEntry;
		adjAEntry = A[0] * A[8] - A[2] * A[6]; detB += adjAEntry * adjAEntry;
		adjAEntry = A[0] * A[7] - A[1] * A[6]; detB += adjAEntry * adjAEntry;
		adjAEntry = A[1] * A[5] - A[2] * A[4]; detB += adjAEntry * adjAEntry;
		adjAEntry = A[0] * A[5] - A[2] * A[3]; detB += adjAEntry * adjAEntry;
		adjAEntry = A[0] * A[4] - A[1] * A[3]; detB += adjAEntry * adjAEntry;
		detB = FT(1) - FT(4) * detB;

		//column first matrix
		B[0] = -A[0] - A[4] - A[8]; B[1] = -A[5] + A[7]; B[2] = -A[6] + A[2]; B[3] = -A[1] + A[3];
		B[4] = B[1]; B[5] = -A[0] + A[4] + A[8]; B[6] = -A[1] - A[3]; B[7] = -A[2] - A[6];
		B[8] = B[2]; B[9] = B[6];  B[10] = -A[4] + A[0] + A[8]; B[11] = -A[5] - A[7];
		B[12] = B[3]; B[13] = B[7]; B[14] = B[11]; B[15] = -A[8] + A[0] + A[4];
			 
		FT detA = FT(1), omega = FT(0);
		constexpr FT tau = FT(1e-4);
		bool quick = detB <= FT(1) - tau;

		if (quick) {
			int rowPermute[2] = { 1, 2 };
			//LU with partial pivoting
			//first column
			{
				int rowIndex = 0;
				FT largest = fabs(A[0]);
				if (fabs(A[3]) > largest) {
					largest = fabs(A[3]);
					rowIndex = 1;
				}
				if (fabs(A[6]) > largest) rowIndex = 2;

				if (rowIndex != 0) {
					detA = -detA;
					rowPermute[rowIndex - 1] = 0;
				}

				int row = rowIndex;
				FT divisor = A[row * 3 + 0];
				detA *= divisor;

				if (detA != FT(0)) {
					for (int j = 0; j < 2; j++) {
						int subRow = rowPermute[j];
						FT lower = A[subRow * 3 + 0] / divisor;
						for (int k = 1; k < 3; k++)
							A[subRow * 3 + k] -= lower * A[row * 3 + k];
					}
				}
			}

			//second column
			if (detA != FT(0)) {
				int row = 0, subRow = 0;
				if (fabs(A[rowPermute[0] * 3 + 1]) > fabs(A[rowPermute[1] * 3 + 1])) {
					row = rowPermute[0];
					subRow = rowPermute[1];
				}
				else {
					detA = -detA;
					row = rowPermute[1];
					subRow = rowPermute[0];
				}
				FT divisor = A[row * 3 + 1];
				detA *= divisor;
				FT lower = A[subRow * 3 + 1] / divisor;
				detA *= (A[subRow * 3 + 2] - lower * A[row * 3 + 2]);
			}
		}
		else {
			int rowPermute[3] = { 0, 1, 2 };
			int columnPermute[3] = { 0, 1, 2 };
			//LU with complete pivoting
			for (int i = 0; i < 2; i++) {
				FT largest = FT(0);
				int rowIndex = i;
				int columnIndex = i;
				for (int j = i; j < 3; j++) {
					int subRow = rowPermute[j];
					for (int k = i; k < 3; k++) {
						int subColumn = columnPermute[k];
						FT entry = fabs(A[subRow * 3 + subColumn]);
						if (entry > largest) {
							largest = entry;
							rowIndex = j;
							columnIndex = k;
						}
					}
				}

				if (rowPermute[i] != rowPermute[rowIndex]) {
					std::swap(rowPermute[i], rowPermute[rowIndex]);
					detA = -detA;
				}
				if (columnPermute[i] != columnPermute[columnIndex]) {
					std::swap(columnPermute[i], columnPermute[columnIndex]);
					detA = -detA;
				}

				int row = rowPermute[i];
				int column = columnPermute[i];

				FT divisor = A[row * 3 + column];
				detA *= divisor;
				if (detA == FT(0)) break;

				for (int j = i + 1; j < 3; j++) {
					int subRow = rowPermute[j];
					FT lower = A[subRow * 3 + column] / divisor;
					for (int k = i + 1; k < 3; k++) {
						int subColumn = columnPermute[k];
						A[subRow * 3 + subColumn] -= lower * A[row * 3 + subColumn];
					}
				}
			}
			detA *= A[rowPermute[2] * 3 + columnPermute[2]];

			omega = fabs(A[rowPermute[1] * 3 + columnPermute[1]]);
		}

		if (detA < FT(0)) {
			detA = -detA;
			for (int i = 0; i < 16; i++) B[i] = -B[i];
		}
		
		constexpr FT inv3 = FT(0.333333333333333333333333333333333333333333333333333333333333333);
		FT domEigen = FT(0);
		if (detB >= tau - inv3) {
			//analytical solution
			FT c = FT(8) * detA;
			FT delta0 = FT(1) + FT(3) * detB;
			FT delta1 = FT(-1) + (FT(27) / FT(16)) * c * c + FT(9) * detB;
			FT sqrtDelta0 = sqrt(delta0);
			FT alpha = delta1 / (delta0 * sqrtDelta0);
			FT z = (FT(4) * inv3) * (FT(1) + sqrtDelta0 * cos(acos(Clamp(alpha, FT(-1), FT(1))) * inv3));
			FT s = sqrt(z) * FT(0.5);
			domEigen = s + FT(0.5) * sqrt(std::max(FT(0), FT(4) - z + c / s));
		}
		else {
			//Newton's method
			domEigen = FT(1.7320508075688772935274463415058723669428052538103806280558069794519330169088);
			FT domEigenOld = domEigen;
			constexpr FT eps = std::is_same<FT, double>::value ? FT(1e-12) : FT(1e-6);

			FT c = FT(8) * detA;
			do {
				FT p = domEigen * (domEigen * (domEigen * domEigen - FT(2)) - c) + detB;
				FT dp = domEigen * (FT(4) * domEigen * domEigen - FT(4)) - c;
				domEigenOld = domEigen;
				domEigen -= p / dp;
			} while (domEigenOld - domEigen > eps);
		}

		B[0] += domEigen; B[5] += domEigen; B[10] += domEigen; B[15] += domEigen;

		//LDLT with symmetric pivoting
		FT L[6];
		int indices[4] = { 0, 1, 2, 3 };
		//the first column
		{
			int diagIndex = 0, trueDiagIndex = 0;
			FT largest = B[0];
			for (int j = 1; j < 4; j++) {
				int trueIndex = indices[j];
				FT diag = B[trueIndex * 4 + trueIndex];
				if (diag > largest) {
					largest = diag;
					diagIndex = j;
				}
			}

			if (diagIndex != 0) std::swap(indices[0], indices[diagIndex]);

			trueDiagIndex = indices[0];
			L[0] = B[trueDiagIndex * 4 + indices[1]] / largest;
			L[1] = B[trueDiagIndex * 4 + indices[2]] / largest;
			L[2] = B[trueDiagIndex * 4 + indices[3]] / largest;

			int rowIndex = 0, colIndex = 0; FT offDiag = FT(0);

			rowIndex = indices[1];
			B[rowIndex * 4 + rowIndex] -= L[0] * B[trueDiagIndex * 4 + rowIndex];
			colIndex = indices[2];
			offDiag = B[rowIndex * 4 + colIndex] - L[0] * B[trueDiagIndex * 4 + colIndex];
			B[rowIndex * 4 + colIndex] = offDiag;
			B[colIndex * 4 + rowIndex] = offDiag;
			colIndex = indices[3];
			offDiag = B[rowIndex * 4 + colIndex] - L[0] * B[trueDiagIndex * 4 + colIndex];
			B[rowIndex * 4 + colIndex] = offDiag;
			B[colIndex * 4 + rowIndex] = offDiag;

			rowIndex = indices[2];
			B[rowIndex * 4 + rowIndex] -= L[1] * B[trueDiagIndex * 4 + rowIndex];
			colIndex = indices[3];
			offDiag = B[rowIndex * 4 + colIndex] - L[1] * B[trueDiagIndex * 4 + colIndex];
			B[rowIndex * 4 + colIndex] = offDiag;
			B[colIndex * 4 + rowIndex] = offDiag;

			rowIndex = indices[3];
			B[rowIndex * 4 + rowIndex] -= L[2] * B[trueDiagIndex * 4 + rowIndex];
		}

		//the second column
		{
			int diagIndex = 1, trueDiagIndex = indices[diagIndex];
			FT largest = B[trueDiagIndex * 4 + trueDiagIndex];
			for (int j = 2; j < 4; j++) {
				int trueIndex = indices[j];
				FT diag = B[trueIndex * 4 + trueIndex];
				if (diag > largest) {
					largest = diag;
					diagIndex = j;
				}
			}

			if (diagIndex != 1) {
				std::swap(indices[1], indices[diagIndex]);
				std::swap(L[0], L[diagIndex - 1]);
			}

			trueDiagIndex = indices[1];
			L[3] = B[trueDiagIndex * 4 + indices[2]] / largest;
			L[4] = B[trueDiagIndex * 4 + indices[3]] / largest;

			int rowIndex = 0, colIndex = 0; FT offDiag = FT(0);

			rowIndex = indices[2];
			B[rowIndex * 4 + rowIndex] -= L[3] * B[trueDiagIndex * 4 + rowIndex];
			colIndex = indices[3];
			offDiag = B[rowIndex * 4 + colIndex] - L[3] * B[trueDiagIndex * 4 + colIndex];
			B[rowIndex * 4 + colIndex] = offDiag;
			B[colIndex * 4 + rowIndex] = offDiag;

			rowIndex = indices[3];
			B[rowIndex * 4 + rowIndex] -= L[4] * B[trueDiagIndex * 4 + rowIndex];
		}

		FT v[4];
		if (quick) {
			//the third column
			if (B[indices[2] * 4 + indices[2]] < B[indices[3] * 4 + indices[3]]) {
				std::swap(indices[2], indices[3]);
				std::swap(L[1], L[2]); std::swap(L[3], L[4]);
			}
			L[5] = B[indices[2] * 4 + indices[3]] / B[indices[2] * 4 + indices[2]];

			v[indices[0]] = L[0] * L[4] + L[1] * L[5] - L[0] * L[5] * L[3] - L[2];
			v[indices[1]] = L[5] * L[3] - L[4];
			v[indices[2]] = -L[5];
			v[indices[3]] = FT(1);
		}
		else {
			FT dd = B[indices[2] * 4 + indices[2]] * B[indices[3] * 4 + indices[3]] - B[indices[2] * 4 + indices[3]] * B[indices[2] * 4 + indices[3]];
			if (dd == FT(0)) {
				if (B[indices[2] * 4 + indices[2]] == FT(0) && B[indices[3] * 4 + indices[3]] == FT(0)
					&& B[indices[2] * 4 + indices[3]] == FT(0)) {
					v[indices[0]] = L[0] * L[4] - L[2];
					v[indices[1]] = -L[4];
					v[indices[2]] = FT(0); v[indices[3]] = FT(1);
				}
				else {
					FT aa = -B[indices[2] * 4 + indices[3]], bb = B[indices[2] * 4 + indices[2]];
					FT norm = sqrt(aa * aa + bb * bb);
					aa /= norm; bb /= norm;
					
					v[indices[3]] = bb; v[indices[2]] = aa;
					v[indices[1]] = -L[4] * bb - L[3] * aa;
					v[indices[0]] = -L[2] * bb - L[1] * aa - L[0] * v[indices[1]];
				}
			}
			else {
				FT ID[5];
				ID[0] = FT(1) / B[indices[0] * 4 + indices[0]];
				ID[1] = FT(1) / B[indices[1] * 4 + indices[1]];
				ID[2] = B[indices[3] * 4 + indices[3]];
				ID[3] = -B[indices[2] * 4 + indices[3]];
				ID[4] = B[indices[2] * 4 + indices[2]];

				FT IL[5];
				IL[0] = -L[0];
				IL[1] = L[0] * L[3] - L[1];
				IL[2] = L[0] * L[4] - L[2];
				IL[3] = -L[3];
				IL[4] = -L[4];

				constexpr FT eps = std::is_same<FT, double>::value ? FT(6.607e-8) : FT(3.388e-4);
				if (omega > eps) {
					FT vv[4] = { IL[2], IL[4], FT(0), FT(1) };

					int nit = std::is_same<FT, double>::value ? (int)ceil(FT(15) / (FT(16.86) + FT(2) * log10(omega))) :
						(int)ceil(FT(7) / (FT(8.12) + FT(2) * log10(omega)));

					for (int i = 0; i < nit; i++) {
						//normalize vv
						FT norm = sqrt(vv[0] * vv[0] + vv[1] * vv[1] + vv[2] * vv[2] + vv[3] * vv[3]);
						for (int i = 0; i < 4; i++) vv[i] /= norm;

						//vv = L^{-1}vv
						FT vv1 = vv[1];
						vv[1] = IL[0] * vv[0] + vv1;
						vv[2] = IL[1] * vv[0] + IL[3] * vv1 + vv[2];
						vv[3] = IL[2] * vv[0] + IL[4] * vv1 + vv[3];

						//vv = D^{-1}vv
						vv[0] *= ID[0]; vv[1] *= ID[1];
						FT vv2 = vv[2];
						vv[2] = (ID[2] * vv2 + ID[3] * vv[3]) / dd;
						vv[3] = (ID[3] * vv2 + ID[4] * vv[3]) / dd;

						//vv = L^{-T}vv
						vv[0] = vv[0] + IL[0] * vv[1] + IL[1] * vv[2] + IL[2] * vv[3];
						vv[1] = vv[1] + IL[3] * vv[2] + IL[4] * vv[3];
					}

					v[indices[0]] = vv[0]; v[indices[1]] = vv[1]; v[indices[2]] = vv[2]; v[indices[3]] = vv[3];
				}
				else {
					FT vv[8];
					vv[0] = IL[1]; vv[1] = IL[3]; vv[2] = FT(1); vv[3] = FT(0);
					vv[4] = IL[2]; vv[5] = IL[4]; vv[6] = FT(0); vv[7] = FT(1);

					Orthonormalize(vv);
					for (int i = 0; i < 2; i++) {
						//vv = L^{-1}vv
						FT vv1 = vv[1];
						vv[1] = IL[0] * vv[0] + vv1;
						vv[2] = IL[1] * vv[0] + IL[3] * vv1 + vv[2];
						vv[3] = IL[2] * vv[0] + IL[4] * vv1 + vv[3];

						FT vv5 = vv[5];
						vv[5] = IL[0] * vv[4] + vv5;
						vv[6] = IL[1] * vv[4] + IL[3] * vv5 + vv[6];
						vv[7] = IL[2] * vv[4] + IL[4] * vv5 + vv[7];

						//vv = D^{-1}vv
						vv[0] *= ID[0]; vv[1] *= ID[1];
						FT vv2 = vv[2];
						vv[2] = (ID[2] * vv2 + ID[3] * vv[3]) / dd;
						vv[3] = (ID[3] * vv2 + ID[4] * vv[3]) / dd;

						vv[4] *= ID[0]; vv[5] *= ID[2];
						FT vv6 = vv[6];
						vv[6] = (ID[2] * vv6 + ID[3] * vv[7]) / dd;
						vv[7] = (ID[3] * vv6 + ID[4] * vv[7]) / dd;

						//vv = L^{-T}vv
						vv[0] = vv[0] + IL[0] * vv[1] + IL[1] * vv[2] + IL[2] * vv[3];
						vv[1] = vv[1] + IL[3] * vv[2] + IL[4] * vv[3];
						vv[4] = vv[4] + IL[0] * vv[5] + IL[1] * vv[6] + IL[2] * vv[7];
						vv[5] = vv[5] + IL[3] * vv[6] + IL[4] * vv[7];
					}
					Orthonormalize(vv);

					FT BB[3], temp0[8], temp1[8];
					temp0[0] = vv[0] + L[0] * vv[1] + L[1] * vv[2] + L[2] * vv[3];
					temp0[1] = vv[1] + L[3] * vv[2] + L[4] * vv[3];
					temp0[2] = vv[2]; temp0[3] = vv[3];
					temp0[4] = vv[4] + L[0] * vv[5] + L[1] * vv[6] + L[2] * vv[7];
					temp0[5] = vv[5] + L[3] * vv[6] + L[4] * vv[7];
					temp0[6] = vv[6]; temp0[7] = vv[7];

					FT d00 = B[indices[0] * 4 + indices[0]], d11 = B[indices[1] * 4 + indices[1]];
					FT d22 = B[indices[2] * 4 + indices[2]], d33 = B[indices[3] * 4 + indices[3]];
					FT d32 = B[indices[2] * 4 + indices[3]];

					temp1[0] = temp0[0] * d00; 
					temp1[1] = temp0[1] * d11;
					temp1[2] = temp0[2] * d22 + temp0[3] * d32; 
					temp1[3] = temp0[2] * d32 + temp0[3] * d33;
					temp1[4] = temp0[4] * d00;
					temp1[5] = temp0[5] * d11;
					temp1[6] = temp0[6] * d22 + temp0[7] * d32; 
					temp1[7] = temp0[6] * d32 + temp0[7] * d33;

					BB[0] = temp1[0] * temp0[0] + temp1[1] * temp0[1] + temp1[2] * temp0[2] + temp1[3] * temp0[3];
					BB[1] = temp1[0] * temp0[4] + temp1[1] * temp0[5] + temp1[2] * temp0[6] + temp1[3] * temp0[7];
					BB[2] = temp1[4] * temp0[4] + temp1[5] * temp0[5] + temp1[6] * temp0[6] + temp1[7] * temp0[7];
					
					constexpr FT eps2 = std::is_same<FT, double>::value ? FT(1e-15) : FT(1e-7);
					if (fabs(BB[1]) > eps2) {
						FT r = (BB[0] - BB[2]) * FT(0.5);
						FT w0 = r - sqrt(r * r + BB[1] * BB[1]);
						v[indices[0]] = vv[0] * w0 + vv[4] * BB[1];
						v[indices[1]] = vv[1] * w0 + vv[5] * BB[1];
						v[indices[2]] = vv[2] * w0 + vv[6] * BB[1];
						v[indices[3]] = vv[3] * w0 + vv[7] * BB[1];
					}
					else {
						if (BB[0] < BB[2]) {
							v[indices[0]] = vv[0]; v[indices[1]] = vv[1]; v[indices[2]] = vv[2]; v[indices[3]] = vv[3];
						}
						else {
							v[indices[0]] = vv[4]; v[indices[1]] = vv[5]; v[indices[2]] = vv[6]; v[indices[3]] = vv[7];
						}
					}
				}
			}
		}

		FT norm = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2] + v[3] * v[3]);
		for (int i = 0; i < 4; i++) v[i] /= norm;

		FT v11 = v[1] * v[1], v22 = v[2] * v[2], v33 = v[3] * v[3];
		FT v12 = v[1] * v[2], v03 = v[0] * v[3], v13 = v[1] * v[3];
		FT v02 = v[0] * v[2], v23 = v[2] * v[3], v01 = v[0] * v[1];
		FT Q[9];
		Q[0] = FT(1) - FT(2) * (v22 + v33);
		Q[1] = FT(2) * (v12 - v03);
		Q[2] = FT(2) * (v13 + v02);
		Q[3] = FT(2) * (v12 + v03);
		Q[4] = FT(1) - FT(2) * (v11 + v33);
		Q[5] = FT(2) * (v23 - v01);
		Q[6] = FT(2) * (v13 - v02);
		Q[7] = FT(2) * (v23 + v01);
		Q[8] = FT(1) - FT(2) * (v11 + v22);

		if (detA < FT(0)) for (int i = 0; i < 9; i++) Q[i] = -Q[i];
		ortho[0] = Q[0]; ortho[1] = Q[3]; ortho[2] = Q[6];
		ortho[3] = Q[1]; ortho[4] = Q[4]; ortho[5] = Q[7];
		ortho[6] = Q[2]; ortho[7] = Q[5]; ortho[8] = Q[8];

		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++)
				sspd[i * 3 + j] = Q[i * 3 + 0] * mat[0 * 3 + j] 
				+ Q[i * 3 + 1] * mat[1 * 3 + j] + Q[i * 3 + 2] * mat[2 * 3 + j];
		}
	}


	template<class FT> FT Dot(int width, const FT *x, const FT *y) {
		FT result = FT(0);
		for (int i = 0; i < width; i++)
			result += x[i] * y[i];
		return result;
	}
}
#endif