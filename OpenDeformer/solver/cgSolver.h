#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_SOLVER_CGSOLVER_H
#define ODER_SOLVER_CGSOLVER_H

#include "oder.h"
#include "linearSolver.h"
#include "sparseMatrix.h"
#include "preconditioner.h"

namespace ODER{
	class CGSolver : public LinearSolver{
	public:
		CGSolver(BlockedSymSpMatrix&& m, double tolerant, const Preconditioner* preconditioner);
		CGSolver(const CGSolver&) = delete;
		CGSolver& operator=(const CGSolver&) = delete;
		CGSolver(CGSolver&&) = default;
		CGSolver& operator=(CGSolver&&) = default;
		void solveLinearSystem(const DenseVector& rhs, DenseVector& result) const;
		~CGSolver() = default;
	private:
		BlockedSymSpMatrix mat;
		double tol;
		const Preconditioner* preconditioner;
	};
}
#endif