#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_CORE_PRECONDITIONER_H
#define ODER_CORE_PRECONDITIONER_H

#include "oder.h"

namespace ODER{
	class Preconditioner{
	public:
		virtual void solvePreconditionerSystem(int width, const double *rhs, double *result) const = 0;
		virtual void resetPreconditionerSystem(const BlockedSymSpMatrix& mat){
			Severe("unimplemented Preconditioner::resetPreconditionerSystem called");
		}
		virtual void resetPreconditionerSystem(const SparseMatrix& mat){
			Severe("unimplemented Preconditioner::resetPreconditionerSystem called");
		}
		virtual void Preprocess(const BlockedSymSpMatrix& mat) {
			Severe("unimplemented Preconditioner::Preprocess called");
		}
		virtual void Preprocess(const SparseMatrix& mat) {
			Severe("unimplemented Preconditioner::Preprocess called");
		}
		virtual ~Preconditioner() = default;
	};
}

#endif