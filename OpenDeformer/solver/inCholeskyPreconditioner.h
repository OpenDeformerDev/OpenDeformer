#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_SOLVER_INCHOLESKYPRECONDITIONER_H
#define ODER_SOLVER_INCHOLESKYPRECONDITIONER_H

#include "preconditioner.h"
#include <vector>

namespace ODER{
	class InCholeskyPreconditioner : public Preconditioner{
	public:
		InCholeskyPreconditioner(Scalar epsilon, Scalar relaxedScale) : epsilon(epsilon), relaxedScale(relaxedScale){}
		InCholeskyPreconditioner(const BlockedSymSpMatrix& mat, Scalar epsilon, Scalar relaxedScale);
		InCholeskyPreconditioner(const InCholeskyPreconditioner&&) = delete;
		InCholeskyPreconditioner& operator=(const InCholeskyPreconditioner&&) = delete;
		void resetPreconditionerSystem(const BlockedSymSpMatrix& mat);
		void solvePreconditionerSystem(int width, const Scalar *rhs, Scalar *result) const;
		void Preprocess(const BlockedSymSpMatrix& mat);
		~InCholeskyPreconditioner() = default;
	private:
		void incompleteCholeskyDecomposition(const BlockedSymSpMatrix& mat);
		void proccessSingleColumn(int columnIndex, const FastSparseVector& vec, Scalar *diags, 
			std::pair<int, int> *list, std::vector<std::pair<int, Scalar>>& factorized);

		std::vector<Scalar> values;
		std::vector<int> rows;
		std::vector<int> pcol;
		Scalar epsilon;
		Scalar relaxedScale;
	};
}

#endif