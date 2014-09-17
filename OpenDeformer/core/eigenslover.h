#pragma once
#ifndef ODER_EIGENSLOVER
#define ODER_EIGENSLOVER

#include "oder.h"
#include "sparseMatrix.h"

namespace ODER{
	class EigenSlover{
	public:
		void getEigenValVectors(int dofs, const SparseMatrix &M, const SparseMatrix &K, double *EigValues, double *EigVectors) const;
	};
}

#endif
