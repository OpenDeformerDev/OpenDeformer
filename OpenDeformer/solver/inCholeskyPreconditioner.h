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
		InCholeskyPreconditioner(const BlockedSymSpMatrix& mat, double epsilon);
		InCholeskyPreconditioner(const InCholeskyPreconditioner&&) = delete;
		InCholeskyPreconditioner& operator=(const InCholeskyPreconditioner&&) = delete;
		void resetPreconditionerSystem(const BlockedSymSpMatrix& mat);
		void solvePreconditionerSystem(const DenseVector& rhs, DenseVector& result) const;
		~InCholeskyPreconditioner() = default;
	private:
		void incompleteCholeskyDecomposition(const BlockedSymSpMatrix& mat);
		void proccessSingleColumn(int columnIndex, const SparseVector& vec, double *diags);


		std::vector<double> values;
		std::vector<int> rows;
		std::vector<int> pcol;
		double epsilon;
	};
}

#endif