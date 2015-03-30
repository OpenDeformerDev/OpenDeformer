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

	template<int blockLength, int blockWidth, class LhsVec, class RhsVec>
	void SpMDV(const BlockedSymSparseMatrix<blockLength, blockWidth>& mat, 
		const LhsVec& src, RhsVec& dest){
		static_assert(std::is_same<std::remove_const_t<LhsVec>, double *>::value || 
			std::is_same<std::decay_t<decltype(std::declval<LhsVec>().operator[](std::declval<int>()))>, double>::value,
			"only type with operator[] return double or double pointer supported for ODER::SpMV src yet");
		static_assert(std::is_same<RhsVec, double *>::value || 
			std::is_same<std::decay_t<decltype(std::declval<RhsVec>().operator[](std::declval<int>()))>, double>::value,
			"only type with operator[] return double or non-const double pointer supported for ODER::SpMV dest yet");

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
						double val = values[diagIndex++] * x;
						dest[row + k] += val;
						ret += val;
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
		const SparseVector& src, SparseVector& dest){
		constexpr int regularSize = blockLength * blockWidth;
		constexpr int diagSize = (blockLength + 1) * blockLength / 2;

		const int blockColumnCount = mat.numBlockColumn;
		const int columnCount = mat.numColumns;

		const double *values = mat.values;
		const int *blockPcol = mat.blockPcol;
		const int *blockRows = mat.blockRows;
		const int *blockColumnOris = mat.blockColumnOris;

		for (auto indexValPair : src){
			const int column = indexValPair.first;
			const double entry = indexValPair.second;

			int blockIndex = column / blockLength;
			int blockStartColumn = blockIndex * blockLength;
			for (int i = 0, startColumn = 0; i < blockIndex; i++, startColumn += blockLength){
				int start = blockPcol[i], end = blockPcol[i + 1];
				if (start != end){
					int blockRow = ((column - blockLength) / blockWidth) * blockWidth + blockLength;
					int offset = column - blockRow;
					const int *blockRowIndexStart = &blockRows[start], *blockRowIndexEnd = &blockRows[end];

					auto iter = std::lower_bound(blockRowIndexStart, blockRowIndexEnd, blockRow);
					if (*iter == blockRow){
						const double* block = values + blockColumnOris[i];
						block += diagSize + (iter - blockRowIndexStart - 1)*regularSize;
						if (*blockRowIndexStart != startColumn)
							block += regularSize - diagSize;

						for (int j = 0; j < blockLength; j++)
							dest.Add(startColumn + j, block[j * blockWidth + offset] * entry);
					}
				}
			}

			if (blockIndex < blockColumnCount){
				int start = blockPcol[blockIndex], end = blockPcol[blockIndex + 1];
				if (start != end){
					const int *blockRowIndex = &blockRows[start];
					int offset = column - blockStartColumn;
					const double* block = values + blockColumnOris[blockIndex];
					if (*blockRowIndex == blockStartColumn){
						for (int i = offset; i > 0; i--)
							dest.Add(column - i, block[mat.diagIndices[offset - i] + i] * entry);

						int entryCount = blockLength - offset;
						for (int i = 0; i < entryCount; i++)
							dest.Add(column + i, block[mat.diagIndices[offset] + i] * entry);
						block += diagSize;
						blockRowIndex++;
					}
					const int *blockRowEnd = &blockRows[end];
					if (blockRowIndex == blockRowEnd) continue;

					int adjust = offset * blockWidth;
					while (blockRowIndex != blockRowEnd - 1){
						int row = *blockRowIndex;
						for (int i = 0; i < blockWidth; i++)
							dest.Add(row + i, block[adjust + i] * entry);
						block += regularSize;
						blockRowIndex++;
					}

					int row = *blockRowIndex;
					int mayDegenWidth = columnCount - row;
					if (mayDegenWidth >= blockWidth){
						for (int i = 0; i < blockWidth; i++)
							dest.Add(row + i, block[adjust + i] * entry);
					}
					else{
						adjust = offset * mayDegenWidth;
						for (int i = 0; i < mayDegenWidth; i++)
							dest.Add(row + i, block[adjust + i] * entry);
					}
				}
			}
			else{
				const int index = blockIndex + column - blockStartColumn;
				const int start = blockPcol[index], end = blockPcol[index + 1];

				const double *vals = values + blockColumnOris[index];
				const int *rowIndex = &blockRows[start];

				for (int i = 0; i < end - start; i++)
					dest.Add(rowIndex[i], vals[i] * entry);
			}

		}
	}
}
#endif