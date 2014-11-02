#include "stdafx.h"
#include "eigenslover.h"
#include "arlgsym.h"
#include "arlssym.h"

using namespace ODER;

void EigenSlover::getEigenValVectors(int dofs, const SparseMatrix &M, const SparseMatrix &K, double *EigValues, double *EigVectors) const{
	ARluSymMatrix<double> mass(M.numColumns, M.pcol[M.numColumns], M.value, M.rows, M.pcol);
	ARluSymMatrix<double> stiffness(K.numColumns, K.pcol[K.numColumns], K.value, K.rows, K.pcol);
	ARluSymGenEig<double> slover(dofs, stiffness, mass);
	slover.SetShiftInvertMode(0.0);
	slover.EigenValVectors(EigVectors, EigValues);
}