#include "stdafx.h"
#include "cgSolver.h"
#include "numerMethod.h"
#include "latool.h"

namespace ODER{
	void CGSolver::resetLinearSystem(const BlockedSymSpMatrix* m) {
		this->mat = m;
		if (preconditioner)
			preconditioner->resetPreconditionerSystem(*this->mat);
	}

	void CGSolver::solveLinearSystem(const Scalar *rhs, Scalar *result) {
		const int width = this->mat->getNumColumns();
		Scalar *memory = new Scalar[3 * width];
		Initiation(memory, 3 * width);
		Scalar *remainder = memory, *direction = memory + width, *temp = memory + 2 * width;
		// remainder = rhs - mat * result
		SpMDV(*(this->mat), result, remainder);
		for (int i = 0; i < width; i++)
			remainder[i] = rhs[i] - remainder[i];

		preconditioner->Preprocess(*this->mat);
		//M * d = r
		preconditioner->solvePreconditionerSystem(width, remainder, direction);
		Scalar delta = Dot(width, remainder, direction);

		Scalar epsilon = Dot(width, rhs, rhs) * tolerant * tolerant;
		Scalar remainderNorm2 = Dot(width, remainder, remainder);

		int iter = 0;
		while (iter++ < maxIteration && remainderNorm2 > epsilon) {
			// q = mat * direction
			Initiation(temp, width);
			SpMDV(*(this->mat), direction, temp);
			Scalar alpha = delta / Dot(width, direction, temp);
			//result = result + alpha * direction
			//remainder = remainder - alpha * q
			for (int i = 0; i < width; i++){
				result[i] += alpha * direction[i];
				remainder[i] -= alpha * temp[i];
			}

			// M * z = remainder
			preconditioner->solvePreconditionerSystem(width, remainder, temp);
			Scalar preDelta = delta;
			delta = Dot(width, remainder, temp);
			Scalar beta = delta / preDelta;
			//direction = beta * direction + z
			for (int i = 0; i < width; i++)
				direction[i] = beta * direction[i] + temp[i];

			remainderNorm2 = Dot(width, remainder, remainder);
		}

		delete[] memory;
	}
}