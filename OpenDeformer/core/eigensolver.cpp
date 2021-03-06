#include "stdafx.h"
#include "eigenSolver.h"
#include "sparseMatrix.h"
#include "arlgsym.h"
#include "arlssym.h"

namespace ODER{
	void EigenSolver::getEigenValVectors(int dofs, const SparseMatrix &M, const SparseMatrix &K, Scalar *EigValues, Scalar *EigVectors) const{
		ARluSymMatrix<Scalar> mass(M.numColumns, M.pcol[M.numColumns], M.value, M.rows, M.pcol);
		ARluSymMatrix<Scalar> stiffness(K.numColumns, K.pcol[K.numColumns], K.value, K.rows, K.pcol);
		ARluSymGenEig<Scalar> slover(dofs, stiffness, mass);
		slover.SetShiftInvertMode(0.0);
		slover.EigenValVectors(EigVectors, EigValues);
	}
}