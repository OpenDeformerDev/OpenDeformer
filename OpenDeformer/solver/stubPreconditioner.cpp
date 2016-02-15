#include "stdafx.h"
#include "stubPreconditioner.h"

namespace ODER {
	void StubPreconditioner::solvePreconditionerSystem(int width, const double *rhs, double *result) const {
		memcpy(result, rhs, width * sizeof(double));
	}
}