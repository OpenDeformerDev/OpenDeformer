#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_SOLVER_STUBPRECONDITIONER_H
#define ODER_SOLVER_STUBPRECONDITIONER_H

#include "preconditioner.h"

namespace ODER {
	class StubPreconditioner : public Preconditioner {
	public:
		void solvePreconditionerSystem(int width, const double *rhs, double *result) const;
		void resetPreconditionerSystem(const BlockedSymSpMatrix& mat) {}
		void resetPreconditionerSystem(const SparseMatrix& mat) {}
		void Preprocess(const BlockedSymSpMatrix& mat) {}
		void Preprocess(const SparseMatrix& mat) {}
		~StubPreconditioner() = default;
	};
}

#endif
