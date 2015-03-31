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
		InLDLTPreconditioner(const BlockedSymSpMatrix& mat, double sainvEpsilon, double ldltEpsilon);
		InLDLTPreconditioner(const InLDLTPreconditioner&) = delete;
		InLDLTPreconditioner& operator=(const InLDLTPreconditioner&&) = delete;

		void resetPreconditionerSystem(const BlockedSymSpMatrix& mat);
		void solvePreconditionerSystem(const DenseVector& rhs, DenseVector& result) const;
		~InLDLTPreconditioner(){ delete[] invDiagonal; }
	private:
		void incompleteLDLTDecomposition(const BlockedSymSpMatrix& mat);

		std::vector<double> values;
		std::vector<int> rows;
		std::vector<int> pcol;
		double *invDiagonal;

		double sainvEpsilon;
		double ldltEpsilon;
	};
}


#endif