#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_CORE_EIGENSOLVER_H
#define ODER_CORE_EIGENSOLVER_H

#include "oder.h"

namespace ODER{
	class EigenSolver{
	public:
		void getEigenValVectors(int dofs, const SparseMatrix &M, const SparseMatrix &K, Scalar *EigValues, Scalar *EigVectors) const;
	};
}

#endif
