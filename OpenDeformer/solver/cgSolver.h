#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_SOLVER_CGSOLVER_H
#define ODER_SOLVER_CGSOLVER_H

#include "oder.h"
#include "linearSolver.h"
#include "sparseMatrix.h"

namespace ODER{
	class CGSolver : public LinearSolver<BlockedSymSpMatrix>{
	public:
		CGSolver(Scalar tol, int maxIter, Preconditioner* precondition = NULL, const BlockedSymSpMatrix* m = NULL)
			:LinearSolver(m), tolerant(tol), maxIteration(maxIter), preconditioner(precondition){}
		CGSolver(const CGSolver&) = delete;
		CGSolver& operator=(const CGSolver&) = delete;
		CGSolver(CGSolver&&) = default;
		CGSolver& operator=(CGSolver&&) = default;
		void resetLinearSystem(const BlockedSymSpMatrix* m);
		void solveLinearSystem(const Scalar *rhs, Scalar *result) const;
		~CGSolver() = default;
	private:
		Scalar tolerant;
		int maxIteration;
		Preconditioner* preconditioner;
	};
}
#endif