#include "stdafx.h"
#include "inLDLTPreconditioner.h"
#include "latool.h"
#include "sparseMatrix.h"
#include "numerMethod.h"

namespace ODER{
	InLDLTPreconditioner::InLDLTPreconditioner(const BlockedSymSpMatrix& mat, double sainvEps, double ldltEps)
	: sainvEpsilon(sainvEpsilon), ldltEpsilon(ldltEpsilon), matFullIndices(mat.getFullIndices()) {

		int columnCount = mat.getNumColumns();
		values.reserve(columnCount);
		rows.reserve(columnCount);
		pcol.reserve(columnCount);
		invDiagonal.resize(columnCount);
	}

	void InLDLTPreconditioner::resetPreconditionerSystem(const BlockedSymSpMatrix& mat){
		int columnCount = mat.getNumColumns();
		values.reserve(columnCount);
		rows.reserve(columnCount);
		pcol.reserve(columnCount);
		invDiagonal.resize(columnCount);
		matFullIndices = mat.getFullIndices();
	}

	void InLDLTPreconditioner::Preprocess(const BlockedSymSpMatrix& mat) {
		values.clear();
		rows.clear();
		pcol.clear();
		incompleteLDLTDecomposition(mat);
	}

	void InLDLTPreconditioner::incompleteLDLTDecomposition(const BlockedSymSpMatrix& mat){
		int columnCount = mat.getNumColumns();
		pcol.push_back(0);

		double *factors = &invDiagonal[0];
		SparseVector *vecs = new SparseVector[columnCount];

		//work space
		FastSparseVector temp(columnCount);

		for (int i = 0; i < columnCount; i++) {
			vecs[i].Reserve(32);
			vecs[i].emplaceBack(i, 1.0);
		}

		int count = 0;
		for (int i = 0; i < columnCount - 1; i++){
			SpMSV(mat, matFullIndices, vecs[i], temp);

			double diag = vecs[i] * temp;
			double inv = 1.0 / diag;
			invDiagonal[i] = inv;

			for (int j = i + 1; j < columnCount; j++){
				double entry = vecs[j] * temp;
				double lower = entry * inv;
				if (fabs(lower) > ldltEpsilon){
					values.push_back(lower);
					rows.push_back(j);
					++count;
				}
				factors[j] = entry;
			}
			pcol.push_back(count);

			//clean work space
			temp.Clear();

			for (int j = i + 1; j < columnCount; j++){
				double factor = factors[j];
				if (factor != 0){
					factor *= inv;
					//copy vecs[j];
					for (auto iter = vecs[j].cbegin(); iter != vecs[j].cend(); ++iter)
						temp.Set(*(iter.indexIterator), *(iter.valueIterator));

					for (auto iter = vecs[i].cbegin(); iter != vecs[i].cend(); ++iter)
						temp.Add(*(iter.indexIterator), -factor * (*(iter.valueIterator)));

					//drop vecs[j]
					vecs[j].Clear();
					for (auto colIter = temp.indexBegin(); colIter != temp.indexEnd(); ++colIter) {
						int col = *colIter;
						double val = temp[col];
						if (fabs(val) > sainvEpsilon) vecs[j].emplaceBack(col, val);
					}

					temp.Clear();
				}
			}
		}

		int lastColumn = columnCount - 1;
		SpMSV(mat, matFullIndices, vecs[lastColumn], temp);
		invDiagonal[lastColumn] = 1.0 / (vecs[lastColumn] * temp);

		delete[] vecs;
	}

	void InLDLTPreconditioner::solvePreconditionerSystem(int width, const double *rhs, double *result) const{
		Assert(width != 0);

		int end = pcol[1];
		double sub = rhs[0];
		result[0] = sub;
		for (int j = 1; j < width; j++)
			result[j] = rhs[j];

		for (int j = 0; j < end; j++){
			int row = rows[j];
			result[row] = rhs[row] - values[j] * sub;
		}

		for (int column = 1; column < width - 1; column++){
			int start = pcol[column], end = pcol[column + 1];
			double sub = result[column];
			for (int row = start; row < end; row++){
				result[rows[row]] -= values[row] * sub;
			}
		}

		for (int i = 0; i < width; i++)
			result[i] *= invDiagonal[i];

		for (int row = width - 2; row >= 0; row--){
			double dot = 0.0;
			int start = pcol[row], end = pcol[row + 1];
			for (int column = start; column < end; column++)
				dot += result[rows[column]] * values[column];

			result[row] -= dot;
		}
	}

}