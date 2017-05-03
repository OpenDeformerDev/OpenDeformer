#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_SOLVER_INLDLTPRECONDITIONER_H
#define ODER_SOLVER_INLDLTPRECONDITIONER_H

#include "preconditioner.h"
#include <vector>

namespace ODER{
	class InLDLTPreconditioner : public Preconditioner{
	public:
		InLDLTPreconditioner(Scalar sainvEpsilon, Scalar ldltEpsilon)
			:sainvEpsilon(sainvEpsilon), ldltEpsilon(ldltEpsilon) {}
		InLDLTPreconditioner(const BlockedSymSpMatrix& mat, Scalar sainvEpsilon, Scalar ldltEpsilon);
		InLDLTPreconditioner(const InLDLTPreconditioner&) = delete;
		InLDLTPreconditioner& operator=(const InLDLTPreconditioner&&) = delete;

		void resetPreconditionerSystem(const BlockedSymSpMatrix& mat);
		void solvePreconditionerSystem(int width, const Scalar *rhs, Scalar *result) const;
		void Preprocess(const BlockedSymSpMatrix& mat);

	private:
		void incompleteLDLTDecomposition(const BlockedSymSpMatrix& mat);

		struct InvMatRowListNode {
			InvMatRowListNode() : prev(NULL), next(NULL), row(-1), col(-1), value(0.0) {}
			InvMatRowListNode *prev;
			InvMatRowListNode *next;

			int row, col;
			Scalar value;
		};

		std::vector<Scalar> values;
		std::vector<int> rows;
		std::vector<int> pcol;
		std::vector<Scalar> invDiagonal;
		
		std::vector<std::vector<std::pair<int, int>>> matFullIndices;
		Scalar sainvEpsilon;
		Scalar ldltEpsilon;
	};
}


#endif