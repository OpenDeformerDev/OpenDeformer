#include "stdafx.h"
#include "jacobiPreconditioner.h"
#include "sparseMatrix.h"

namespace ODER {
	JacobiPreconditioner::JacobiPreconditioner(const BlockedSymSpMatrix& mat) :invDiags(mat.getNumColumns()) {}

	void JacobiPreconditioner::resetPreconditionerSystem(const BlockedSymSpMatrix& mat) {
		invDiags.resize(mat.getNumColumns());
	}

	void JacobiPreconditioner::Preprocess(const BlockedSymSpMatrix& mat) {
		getInvJacobiMat(mat);
	}

	void JacobiPreconditioner::getInvJacobiMat(const BlockedSymSpMatrix& mat) {
		mat.getDiagonal(&invDiags[0]);
		int columnCount = mat.getNumColumns();
		for (int i = 0; i < columnCount; i++)
			invDiags[i] = invDiags[i] > Scalar(0) ? Scalar(1.0) / invDiags[i] : Scalar(1.0);
	}

	void JacobiPreconditioner::solvePreconditionerSystem(int width, const Scalar *rhs, Scalar *result) const {
		for (int i = 0; i < width; i++)
			result[i] = invDiags[i] * rhs[i];
	}
}