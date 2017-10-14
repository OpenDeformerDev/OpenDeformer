#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_CORE_LINEARSLOVER_H
#define ODER_CORE_LINEARSLOVER_H

#include "oder.h"
#include "memory.h"

namespace ODER{
    template<class SpMatrix> class LinearSolver : public ReferenceCounted {
	public:
		LinearSolver(const SpMatrix* m = NULL) : mat(m){}
		virtual void resetLinearSystem(const SpMatrix* m){ mat = m; }
		virtual void solveLinearSystem(const Scalar *rhs, Scalar *result) = 0;
		virtual ~LinearSolver() = default;
	protected:
		const SpMatrix *mat;
	};
}

#endif