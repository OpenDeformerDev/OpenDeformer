#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_SOLVER_JACOBIPRECONDITIONER_H
#define ODER_SOLVER_JACOBIPRECONDITIONER_H

#include "preconditioner.h"
#include <vector>

namespace ODER {
	class JacobiPreconditioner : public Preconditioner {
	public:
		JacobiPreconditioner() = default;
		JacobiPreconditioner(const BlockedSymSpMatrix& mat);
		void solvePreconditionerSystem(int width, const Scalar *rhs, Scalar *result) const;
		void resetPreconditionerSystem(const BlockedSymSpMatrix& mat);
		void Preprocess(const BlockedSymSpMatrix& mat);
	private:
		void getInvJacobiMat(const BlockedSymSpMatrix& mat);

		std::vector<Scalar> invDiags;
	};
}
#endif
