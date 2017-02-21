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

	void InCholeskyPreconditioner::proccessSingleColumn(int columnIndex, const FastSparseVector& vec, double *diags, 
		std::pair<int, int> *list, std::vector<std::pair<int, double>>& factorized){
		double diag = diags[columnIndex];
		factorized.clear();

		double epsilon2 = epsilon * epsilon;
		auto pairIter = vec.cbegin(), pairEnd = vec.cend();
		while (++pairIter != pairEnd){
			int row = *pairIter.indexIterator;
			double val = *pairIter.valueIterator;
			if (val != 0.0){
				if (val * val < epsilon2 * diag * diags[row]){
					val = fabs(val);
					double scale = sqrt(diag / diags[row]);
					diag += val * scale;
					diags[row] += val / scale;
				}
				else
					factorized.emplace_back(row, val);
			}
		}

		int start = pcol[columnIndex];
		int end = start + 1 + factorized.size();
		pcol.emplace_back(end);

		std::sort(factorized.begin(), factorized.end(), 
			[](const std::pair<int, double>& left, const std::pair<int, double>& right) { return left.first < right.first; });

		if (!factorized.empty()) {
			int row = factorized[0].first;
			list[columnIndex] = list[row];
			list[row] = std::make_pair(columnIndex, start + 1);
		}

		double diagRoot = sqrt(diag);
		values.emplace_back(diagRoot);
		rows.emplace_back(columnIndex);

		double invDiagRoot = 1.0 / diagRoot;
		for (auto pair : factorized) {
			int row = pair.first;
			double val = pair.second * invDiagRoot;

			rows.emplace_back(row);
			values.emplace_back(pair.second * invDiagRoot);
			diags[row] -= val * val;
		}
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

		std::pair<int, int> *list = new std::pair<int, int>[columnCount];
		Initiation(list, columnCount, std::make_pair(-1, -1));
		std::vector<std::pair<int, double>> factorized;
		factorized.reserve(32);

		FastSparseVector w(columnCount);
		mat.getColumn(0, w);
		proccessSingleColumn(0, w, diags, list, factorized);

		for (int i = 1; i < columnCount - 1; i++){
			mat.getColumn(i, w);

			auto pair = list[i];
			while (pair.first >= 0) {
				int j = pair.first;
				int trueIndex = pair.second;
				int nextIndex = trueIndex + 1;

				double upper = values[trueIndex];
				for (int k = nextIndex; k < pcol[j + 1]; k++)
					w.Add(rows[k], -upper * values[k]);

				pair = list[j];
				if (nextIndex < pcol[j + 1]) {
					int row = rows[nextIndex];
					list[j] = list[row];
					list[row] = std::make_pair(j, nextIndex);
				}
			}

			proccessSingleColumn(i, w, diags, list, factorized);
		}

		//proccess the last column
		values.emplace_back(sqrt(diags[columnCount - 1]));
		rows.emplace_back(columnCount - 1);
		pcol.emplace_back(pcol[columnCount - 1] + 1);

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