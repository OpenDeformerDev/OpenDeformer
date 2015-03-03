#include "cgSolver.h"
#include "numerMethod.h"
#include "latool.h"

namespace ODER{
	void CGSolver::solveLinearSystem(const DenseVector& rhs, DenseVector& result) const{
		const int width = mat.getNumColumns();
		DenseVector remainder(width), direction(width), temp(width);
		// remainder = rhs - mat * result
		SpMV(mat, result, remainder);
		for (int i = 0; i < width; i++)
			remainder[i] = rhs[i] - remainder[i];
		//M * d = r
		preconditioner->solvePreconditionerSystem(remainder, direction);
		double delta = remainder * direction;
		double epsilon = delta * tol;

		for (int i = 0; i < width && delta > epsilon; i++){
			// q = mat * direction
			SpMV(mat, remainder, temp);
			double alpha = delta / (direction * temp);
			//result = result + alpha * direction
			for (int i = 0; i < width; i++)
				result[i] += alpha * direction[i];
			//remainder = remainder - alpha * q
			for (int i = 0; i < width; i++)
				remainder[i] -= alpha * temp[i];

			// M * z = remainder
			preconditioner->solvePreconditionerSystem(remainder, temp);
			double preDelta = delta;
			delta = remainder * temp;
			double beta = delta / preDelta;
			//direction = beta * direction + z
			for (int i = 0; i < width; i++)
				direction[i] = beta * direction[i] + temp[i];
		}
	}
}