#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_CORE_NUMERMETHOD_H
#define ODER_CORE_NUMERMETHOD_H

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

				if (FT(2.0)*p < min(FT(3.0)*q*m - fabs(tol*q), fabs(e*q))){
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
	void SpMV(const BlockedSymSparseMatrix<blockLength, blockWidth>& mat,
		const DenseVector& src, DenseVector& dest){
		Assert(src.getWidth() == dest.getWidth() && src.getWidth() == mat.getNumColumns());

		const int blockColumnCount = mat.numBlockColumn;
		const int columnCount = mat.numColumns;
		const int remainedColumnCount = mat.numRemainedColumn;

		const double *values = mat.values;
		const int *blockPcol = mat.blockPcol;

		constexpr int regularSize = blockLength*blockWidth;
		constexpr int diagSize = (blockLength + 1) * blockLength / 2;
		constexpr int entryCount = blockLength > blockWidth ? blockLength : blockWidth;

		double destBlock[blockWidth];
		double srcBlock[blockWidth];

		const int blockBoundary = blockLength*blockColumnCount;
		for (int startColumn = 0; startColumn < blockBoundary; startColumn += blockLength){
			int index = *blockPcol++;
			int end = *blockPcol;

			if (end - index == 0) continue;

			int row = mat.blockRows[index];

			if (row == startColumn){
				int diagIndex = 0;
				int inf = row + blockLength;
				for (int j = row; j < inf; j++){
					double x = src[j];
					double ret = values[diagIndex++] * x;
					for (int k = j + 1; k < inf; k++){
						double val = values[diagIndex++] * x;
						dest[k] += val;
						ret += val;
					}
					dest[j] += ret;
				}
				index++;
				values += diagSize;
			}

			if (index <= end){
				for (int j = 0; j < blockWidth; j++){
					srcBlock[j] = src[startColumn + j];
					destBlock[j] = 0.0;
				}
			}
			else continue;

			while (index != end){
				row = mat.blockRows[index++];
				for (int j = 0, offset = 0; j < blockWidth; j++, offset += blockWidth){
					int actualRow = row + j;
					double x = src[actualRow];
					double ret = 0.0;
					for (int k = 0; k < blockWidth; k++){
						double val = values[offset + k];
						dest[actualRow] += val * x;
						ret += val * srcBlock[k];
					}
					destBlock[j] += ret;
				}
				values += regularSize;
			}

			row = mat.blockRows[index];
			int mayDegenWidth = row - columnCount;
			if (mayDegenWidth >= blockWidth){
				for (int j = 0, offset = 0; j < blockLength; j++, offset += blockWidth){
					int actualRow = row + j;
					double x = src[actualRow];
					double ret = 0.0;
					for (int k = 0; k < blockWidth; k++){
						double val = values[offset + k];
						dest[actualRow] += val * x;
						ret += val * srcBlock[k];
					}
					destBlock[j] += ret;
					values += regularSize;
				}
			}
			else{
				for (int j = 0, offset = 0; j < blockLength; j++, offset += mayDegenWidth){
					int actualRow = row + j;
					double x = src[actualRow];
					double ret = 0.0;
					for (int k = 0; k < mayDegenWidth; k++){
						double val = values[offset + k];
						dest[actualRow] += val * x;
						ret += val * srcBlock[k];
					}
					destBlock[j] += ret;
					values += blockLength*mayDegenWidth;
				}
			}

			for (int j = 0; j < blockWidth; j++)
				dest[startColumn + j] += destBlock[j];
		}

		for (int i = blockBoundary; i < columnCount; i++, i++){
			int index = *blockPcol++;
			int end = *blockPcol;

			int row = mat.blockRows[index];
			double y = 0.0, x = src[i];
			if (row == i){
				y += x * (*values++);
				index++;
			}

			while (index < end){
				row = mat.blockRows[index++];
				double val = *values++;

				y += val * src[row];
				dest[row] += val * x;
			}

			dest[i] += y;

		}
	}
}
#endif