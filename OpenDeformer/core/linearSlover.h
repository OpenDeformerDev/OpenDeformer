#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_CORE_LINEARSLOVER
#define ODER_CORE_LINEARSLOVER

#include "oder.h"

namespace ODER{
    class LinearSolver{
	public:
		virtual void sloveLinearSystem(const SparseMatrix& mat, DenseVector& result){
			Severe("Unimplemented LinearSlover::sloveLinearSystem for SparseMatrix");
		}
		virtual void sloveLinearSystem(const BlockedSymSpMatrix &mat, DenseVector& result){
			Severe("Unimplemented LinearSlover::sloveLinearSystem for BlockedSymSpMatrix");
		}
		virtual ~LinearSolver() = default;
	private:
	};
}

#endif