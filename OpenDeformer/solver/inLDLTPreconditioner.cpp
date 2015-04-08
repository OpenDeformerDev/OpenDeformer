#include "stdafx.h"
#include "inLDLTPreconditioner.h"
#include "latool.h"
#include "sparseMatrix.h"
#include "numerMethod.h"

namespace ODER{
	InLDLTPreconditioner::InLDLTPreconditioner(const BlockedSymSpMatrix& mat, double sainvEps, double ldltEps){
		sainvEpsilon = sainvEps; ldltEpsilon = ldltEps;

		int columnCount = mat.getNumColumns();
		values.reserve(columnCount);
		rows.reserve(columnCount);
		pcol.reserve(columnCount);

		invDiagonal = new double[columnCount];

		incompleteLDLTDecomposition(mat);
	}

	void InLDLTPreconditioner::resetPreconditionerSystem(const BlockedSymSpMatrix& mat){
		values.clear();
		rows.clear();
		pcol.clear();

		incompleteLDLTDecomposition(mat);
	}

	void InLDLTPreconditioner::incompleteLDLTDecomposition(const BlockedSymSpMatrix& mat){
		int columnCount = mat.getNumColumns();
		pcol.push_back(0);

		double *factors = invDiagonal;
		SparseVector *vecs = new SparseVector[columnCount];
		SparseVector temp;

		for (int i = 0; i < columnCount; i++)
			vecs[i].emplaceBack(i, 1.0);

		int count = 0;
		for (int i = 0; i < columnCount - 1; i++){
			SpMSV(mat, vecs[i], temp);
			
			double diag = temp * vecs[i];
			double inv = 1.0 / diag;
			invDiagonal[i] = inv;

			for (int j = i + 1; j < columnCount; j++){
				double entry = temp * vecs[j];
				double lower = entry * inv;
				if (fabs(lower) >= ldltEpsilon){
					values.push_back(lower);
					rows.push_back(j);
					++count;
				}
				factors[j] = entry;
			}
			pcol.push_back(count);

			auto end = vecs[i].cend();
			for (int j = i + 1; j < columnCount; j++){
				double factor = factors[j];
				if (factor != 0.0){
					factor *= inv;
					auto jEnd = vecs[j].end();
					auto jIter = vecs[j].begin();
					for (auto iter = vecs[i].cbegin(); iter != end; ++iter){
						while (jIter != jEnd && jIter->first < iter->first) ++jIter;

						if (jIter != jEnd && jIter->first == iter->first){
							double entry = jIter->second - factor*iter->second;
							if (fabs(entry) > sainvEpsilon)
								jIter->second = entry;
							else
								jIter = vecs[j].Delete(jIter);
						}
						else{
							double entry = -factor*iter->second;
							if (fabs(entry) > sainvEpsilon){
								vecs[j].Set(jIter, iter->first, entry);
							}
						}

					}
				}
			}
			temp.Clear();
		}

		int lastColumn = columnCount - 1;
		SpMSV(mat, vecs[lastColumn], temp);
		invDiagonal[lastColumn] = 1.0 / (temp * vecs[lastColumn]);


		delete[] vecs;
	}



	void InLDLTPreconditioner::solvePreconditionerSystem(const DenseVector& rhs, DenseVector& result) const{
		int width = rhs.getWidth();
		Assert(width != 0);

		int end = pcol[1];
		double sub = rhs[0];
		result[0] = sub;
		for (int j = 1; j < width; j++)
			result[j] = rhs[j];
		for (int j = 1; j < end; j++){
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

		for (int row = width - 2; row > 0; row--){
			double dot = 0.0;
			int start = pcol[row - 1], end = pcol[row];
			for (int column = start; column < end; column++)
				dot += result[rows[column]] * values[column];
			result[rows[row]] -= dot;
		}
	}

}