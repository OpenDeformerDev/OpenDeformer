#include "stdafx.h"
#include "inCholeskyPreconditioner.h"
#include "sparseMatrix.h"
#include "latool.h"

namespace ODER{
	InCholeskyPreconditioner::InCholeskyPreconditioner(const BlockedSymSpMatrix& mat, double eps) : epsilon(eps){
		int columnCount = mat.getNumColumns();
		pcol.reserve(columnCount + 1);
		values.reserve(columnCount);
		rows.reserve(columnCount);

		incompleteCholeskyDecomposition(mat);
	}

	void InCholeskyPreconditioner::proccessSingleColumn(int columnIndex, const SparseVector& vec, double *diags){
		values.emplace_back(0.0);
		double diag = diags[columnIndex];
		auto iter = vec.cbegin(), vecEnd = vec.cend();
		int count = 0;
		while (++iter != vecEnd){
			int row = iter->first;
			double val = iter->second;
			if (val * val < epsilon * diag * diags[row]){
				val = fabs(val);
				double scale = sqrt(diag / diags[row]);
				diag += val * scale;
				diags[row] += val / scale;
			}
			else{
				values.emplace_back(val);
				rows.emplace_back(row);
				count++;
			}
		}
		int start = pcol[columnIndex];
		int end = start + count;
		pcol.emplace_back(end);
		values[start] = sqrt(diag);
		double invDiag = 1.0 / diag;
		for (int i = start; i < end; i++)
			values[i] *= invDiag;
	}

	void InCholeskyPreconditioner::resetPreconditionerSystem(const BlockedSymSpMatrix& mat){
		values.clear();
		rows.clear();
		pcol.clear();
		incompleteCholeskyDecomposition(mat);
	}

	void InCholeskyPreconditioner::incompleteCholeskyDecomposition(const BlockedSymSpMatrix& mat){
		int columnCount = mat.getNumColumns();
		pcol.emplace_back(0);

		double *diags = new double[columnCount];
		mat.getDiagonal(diags);

		SparseVector w;
		mat.getColumn(0, w);
		proccessSingleColumn(0, w, diags);

		for (int i = 1; i < columnCount; i++){
			w.Clear();
			mat.getColumn(i, w);
			for (int j = 0; j < i; j++){
				int* start = &rows[pcol[j]], *end = &rows[pcol[j + 1]];
				int* rowIter = std::lower_bound(start, end, i);
				if (rowIter != end && *rowIter == i){
					int offset = pcol[j] + rowIter - start;
					double upper = values[offset];
					while (rowIter < end)
						w.Add(*rowIter++, -upper * values[offset++]);
				}
				proccessSingleColumn(i, w, diags);
			}
		}
		delete[] diags;
	}

	void InCholeskyPreconditioner::solvePreconditionerSystem(const DenseVector& rhs, DenseVector& result) const{
		int width = rhs.getWidth();
		Assert(width != 0);

		int end = pcol[1];
		double sub = rhs[0] / values[0];
		result[0] = sub;
		for (int j = 1; j < end; j++){
			int row = rows[j];
			result[row] = rhs[row] - values[j] * sub;
		}

		for (int column = 1; column < width; column++){
			int start = pcol[column], end = pcol[column + 1];
			double sub = result[column] / values[start];
			result[column] = sub;
			for (int row = start + 1; row < end; row++){
				result[rows[row]] -= values[row] * sub;
			}
		}

		int lastDiag = width - 1;
		result[lastDiag] /= values[pcol[lastDiag]];

		for (int row = width - 1; row > 0; row--){
			double dot = 0.0;
			int start = pcol[row - 1], end = pcol[row];
			for (int column = start + 1; column < end; column++)
				dot += result[rows[column]] * values[column];
			int actualRow = rows[row];
			result[actualRow] = (result[actualRow] - dot) / values[start];
		}
	}
}