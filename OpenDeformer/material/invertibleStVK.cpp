#include "stdafx.h"
#include "invertibleStVK.h"

namespace ODER {
	void InvertibleStVKMaterial::getEnergyGradient(const double *invariants, double *gradient) const {
		gradient[0] = 0.25 * lambda * (invariants[0] - 3.0) - 0.5 * mu;
		gradient[1] = 0.25 * mu;
		gradient[2] = 0.0;
	}

	void InvertibleStVKMaterial::getEnergyHassian(const double *invariants, double *hassian) const {
		Initiation(hassian, 6);
		hassian[0] = 0.25;
	}
}