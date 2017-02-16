#include "stdafx.h"
#include "inCholeskyPreconditioner.h"
#include "sparseMatrix.h"
#include "latool.h"
#include "datastructure.h"

namespace ODER{
	InCholeskyPreconditioner::InCholeskyPreconditioner(const BlockedSymSpMatrix& mat, double eps) : epsilon(eps){
		int columnCount = mat.getNumColumns();
		pcol.reserve(columnCount + 1);
		values.reserve(columnCount);
		rows.reserve(columnCount);

		incompleteCholeskyDecomposition(mat);
	}

	void InCholeskyPreconditioner::proccessSingleColumn(int columnIndex, const FastSparseVector& vec, double *diags, RecycledList<std::pair<int, int>> *list){
		values.emplace_back(0.0);
		rows.emplace_back(columnIndex);
		auto indexIter = vec.indexBegin(), indexEnd = vec.indexEnd();
		double diag = diags[columnIndex];

		int count = 1;
		while (++indexIter != indexEnd){
			int row = *indexIter;
			double val = vec[row];
			if (val != 0.0){
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
		}

		int start = pcol[columnIndex];
		int end = start + count;
		pcol.emplace_back(end);

		if (count > 1) {
			int row = rows[start + 1];
			list[row - 1].emplace_back(std::make_pair(columnIndex, start + 1));
		}

		double diagRoot = sqrt(diag);
		values[start] = diagRoot;
		double invDiagRoot = 1.0 / diagRoot;
		for (int i = start + 1; i < end; i++)
			values[i] *= invDiagRoot;
	}

	void InCholeskyPreconditioner::resetPreconditionerSystem(const BlockedSymSpMatrix& mat){
		int columnCount = mat.getNumColumns();
		pcol.reserve(columnCount + 1);
		values.reserve(columnCount);
		rows.reserve(columnCount);
	}

	void InCholeskyPreconditioner::Preprocess(const BlockedSymSpMatrix& mat) {
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

		RecycledList<std::pair<int, int>> *list = new RecycledList<std::pair<int, int>>[columnCount - 1];

		FastSparseVector w; 
		mat.getColumn(0, w);
		proccessSingleColumn(0, w, diags, list);

		for (int i = 1; i < columnCount; i++){
			mat.getColumn(i, w);
			w.Set(i, diags[i]);

			auto iter = list[i - 1].cbegin();
			auto end = list[i - 1].cend();

			while (iter != end) {
				int j = iter->first;
				int trueIndex = iter->second;
				int nextIndex = trueIndex + 1;

				double upper = values[trueIndex];
				for (int k = nextIndex; k < pcol[j + 1]; k++)
					w.Add(rows[k], -upper * values[k]);

				if (nextIndex < pcol[j + 1])
					list[rows[nextIndex] - 1].emplace_back(std::make_pair(j, nextIndex));

				iter = list[i - 1].erase(iter);
			}

			proccessSingleColumn(i, w, diags, list);
		}

		delete[] diags;
		delete[] list;
	}

	void InCholeskyPreconditioner::solvePreconditionerSystem(int width, const double *rhs, double *result) const{
		Assert(width != 0);

		int end = pcol[1];
		double sub = rhs[0] / values[0];
		result[0] = sub;
		for (int j = 1; j < width; j++)
			result[j] = rhs[j];
		for (int j = 1; j < end; j++)
			result[rows[j]] -= values[j] * sub;

		for (int column = 1; column < width; column++){
			int start = pcol[column], end = pcol[column + 1];
			double sub = result[column] / values[start];
			result[column] = sub;
			for (int row = start + 1; row < end; row++){
				result[rows[row]] -= values[row] * sub;
			}
		}

		result[width - 1] /= values[pcol[width - 1]];

		for (int row = width - 2; row >= 0; row--){
			double dot = 0.0;
			int start = pcol[row], end = pcol[row + 1];
			for (int column = start + 1; column < end; column++)
				dot += result[rows[column]] * values[column];

			result[row] = (result[row] - dot) / values[start];
		}
	}
}