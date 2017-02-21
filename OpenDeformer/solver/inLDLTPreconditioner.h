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
		InLDLTPreconditioner(double sainvEpsilon, double ldltEpsilon)
			:sainvEpsilon(sainvEpsilon), ldltEpsilon(ldltEpsilon) {}
		InLDLTPreconditioner(const BlockedSymSpMatrix& mat, double sainvEpsilon, double ldltEpsilon);
		InLDLTPreconditioner(const InLDLTPreconditioner&) = delete;
		InLDLTPreconditioner& operator=(const InLDLTPreconditioner&&) = delete;

		void resetPreconditionerSystem(const BlockedSymSpMatrix& mat);
		void solvePreconditionerSystem(int width, const double *rhs, double *result) const;
		void Preprocess(const BlockedSymSpMatrix& mat);

	private:
		void incompleteLDLTDecomposition(const BlockedSymSpMatrix& mat);

		struct InvMatRowListNode {
			InvMatRowListNode() : prev(NULL), next(NULL), row(-1), col(-1), value(0.0) {}
			InvMatRowListNode *prev;
			InvMatRowListNode *next;

			int row, col;
			double value;
		};

		std::vector<double> values;
		std::vector<int> rows;
		std::vector<int> pcol;
		std::vector<double> invDiagonal;
		
		std::vector<std::vector<std::pair<int, int>>> matFullIndices;
		double sainvEpsilon;
		double ldltEpsilon;
	};
}


#endif