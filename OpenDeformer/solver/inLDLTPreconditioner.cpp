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
		MemoryPool<InvMatRowListNode, std::alignment_of<InvMatRowListNode>::value> nodePool(columnCount);
		InvMatRowListNode *listPerRow = new InvMatRowListNode[columnCount];
		std::vector<InvMatRowListNode *> *invRowNodes = new std::vector<InvMatRowListNode *>[columnCount];

		//work space
		FastSparseVector vec(columnCount);
		FastSparseVector lowerColumn(columnCount);
		SparseVector invColumn;
		invColumn.Reserve(32);

		constexpr int preAlloced = 32;
		for (int i = 0; i < columnCount; i++) {
			InvMatRowListNode *node = nodePool.Alloc();
			node->row = i;  node->col = i; node->value = 1.0;
			listPerRow[i].prev = node;
			listPerRow[i].next = node;
			node->prev = &listPerRow[i];
			node->next = &listPerRow[i];

			invRowNodes[i].reserve(std::min(i + 1, preAlloced));
			invRowNodes[i].push_back(node);
		}

		int count = 0;
		pcol.push_back(0);
		for (int i = 0; i < columnCount - 1; i++){
			for (auto node : invRowNodes[i]) {
				invColumn.emplaceBack(node->row, node->value);
				node->prev->next = node->next;
			    node->next->prev = node->prev;
				nodePool.Dealloc(node);
			}
			invRowNodes[i].clear();
			invRowNodes[i].shrink_to_fit();

			SpMSV(mat, matFullIndices, invColumn, vec);
			double diag =  vec * invColumn;
			double inv = 1.0 / diag;
			invDiagonal[i] = inv;

			for (auto iter = vec.cbegin(); iter != vec.cend(); ++iter) {
				int index = *iter.indexIterator;
				double value = *iter.valueIterator;
				if (value != 0) {
					InvMatRowListNode *node = listPerRow[index].next;
					InvMatRowListNode *end = &listPerRow[index];
					while (node != end) {
						lowerColumn.Add(node->col, value * node->value);
						node = node->next;
					}
				}
			}

			for (auto iter = lowerColumn.cbegin(); iter != lowerColumn.cend(); ++iter) {
				int j = *iter.indexIterator;
				double lower = (*iter.valueIterator) * inv;
				if (fabs(lower) > ldltEpsilon){
					values.push_back(lower);
					rows.push_back(j);
					++count;
				}
			}
			pcol.push_back(count);

			vec.Clear();
			for (auto iter = lowerColumn.cbegin(); iter != lowerColumn.cend(); ++iter) {
				int j = *iter.indexIterator;
				double factor = *iter.valueIterator;
				if (factor != 0) {
					factor *= inv;
					//copy invRowNodes[j] and clean node
					for (auto node : invRowNodes[j]) {
						vec.emplaceBack(node->row, node->value);
						node->prev->next = node->next;
						node->next->prev = node->prev;
						nodePool.Dealloc(node);
					}
					invRowNodes[j].clear();

					for (auto iter = invColumn.cbegin(); iter != invColumn.cend(); ++iter)
						vec.Add(*(iter.indexIterator), -factor * (*(iter.valueIterator)));

					//drop invRowNodes[j]
					for (auto iter = vec.cbegin(); iter != vec.cend(); ++iter) {
						int index = *iter.indexIterator;
						double val = *iter.valueIterator;
						if (fabs(val) > sainvEpsilon) {
							//insert to linked list
							InvMatRowListNode *node = nodePool.Alloc();
							node->row = index; node->col = j; node->value = val;

							node->next = listPerRow[index].next;
							node->prev = &listPerRow[index];
							listPerRow[index].next->prev = node;
							listPerRow[index].next = node;
							invRowNodes[j].push_back(node);
						}
					}
					vec.Clear();
				}
			}
			invColumn.Clear();
			lowerColumn.Clear();
		}

		int lastColumn = columnCount - 1;
		for (auto node : invRowNodes[lastColumn])
			invColumn.emplaceBack(node->row, node->value);
		SpMSV(mat, matFullIndices, invColumn, vec);
		invDiagonal[lastColumn] = 1.0 / (vec * invColumn);

		delete[] listPerRow;
		delete[] invRowNodes;
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