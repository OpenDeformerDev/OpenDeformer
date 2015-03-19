#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_SOLVER_INCHOLESKYPRECONDITIONNER_H
#define ODER_SOLVER_INCHOLESKYPRECONDITIONNER_H

#include "preconditioner.h"
#include <vector>

namespace ODER{
	class InCholeskyPreconditioner : public Preconditioner{
	public:
		InCholeskyPreconditioner(const BlockedSymSpMatrix& mat, double epsilon);
		void solvePreconditionerSystem(const DenseVector& rhs, DenseVector& result) const;
	private:
		void proccessSingleColumn(int columnIndex, const SparseVector& vec, double *diags);

		std::vector<double> values;
		std::vector<int> rows;
		std::vector<int> pcol;
		double epsilon;
	};
}

#endif