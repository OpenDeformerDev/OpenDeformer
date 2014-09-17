#if defined(_MSC_VER)
#pragma once
#endif
#ifndef ODER_CORE_EIGENSLOVER_H
#define ODER_CORE_EIGENSLOVER_H

#include "oder.h"
#include "sparseMatrix.h"

namespace ODER{
	class EigenSlover{
	public:
		void getEigenValVectors(int dofs, const SparseMatrix &M, const SparseMatrix &K, double *EigValues, double *EigVectors) const;
	};
}

#endif
