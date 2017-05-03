#include "stdafx.h"
#include "stubPreconditioner.h"

namespace ODER {
	void StubPreconditioner::solvePreconditionerSystem(int width, const Scalar *rhs, Scalar *result) const {
		memcpy(result, rhs, width * sizeof(Scalar));
	}
}