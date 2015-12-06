#include "stdafx.h"
#include "cgSolver.h"
#include "numerMethod.h"
#include "latool.h"
#include "preconditioner.h"

namespace ODER{
	void CGSolver::solveLinearSystem(const double *rhs, double *result) const{
		const int width = this->mat->getNumColumns();
		double *memory = new double[3 * width];
		Initiation(memory, 3 * width);
		double *remainder = memory, *direction = memory + width, *temp = memory + 2 * width;
		// remainder = rhs - mat * result
		SpMDV(*(this->mat), result, remainder);
		for (int i = 0; i < width; i++)
			remainder[i] = rhs[i] - remainder[i];
		//M * d = r
		preconditioner->solvePreconditionerSystem(width, remainder, direction);
		double delta = Dot(width, remainder, direction);

		double epsilon = delta * tolerant;

		while (delta > epsilon){
			// q = mat * direction
			Initiation(temp, width);
			SpMDV(*(this->mat), remainder, temp);
			double alpha = delta / Dot(width, direction, temp);
			//result = result + alpha * direction
			//remainder = remainder - alpha * q
			for (int i = 0; i < width; i++){
				result[i] += alpha * direction[i];
				remainder[i] -= alpha * temp[i];
			}

			// M * z = remainder
			preconditioner->solvePreconditionerSystem(width, remainder, temp);
			double preDelta = delta;
			delta = Dot(width, remainder, temp);
			double beta = delta / preDelta;
			//direction = beta * direction + z
			for (int i = 0; i < width; i++)
				direction[i] = beta * direction[i] + temp[i];
		}

		delete[] memory;
	}
}