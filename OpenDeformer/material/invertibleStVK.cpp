#include "stdafx.h"
#include "invertibleStVK.h"

namespace ODER {
	InvertibleStVKMaterial::InvertibleStVKMaterial(double rho, double youngs, double possion,
		double inversionTrashold) : InvertibleHyperelasticMaterial(rho, inversionTrashold) {
		lambda = possion * youngs / ((1.0 + possion) * (1.0 - 2.0 * possion));
		mu = youngs / (2.0 * (1.0 + possion));
	}

	void InvertibleStVKMaterial::getEnergyGradient(const double *invariants, double *gradient) const {
		gradient[0] = 0.25 * lambda * (invariants[0] - 3.0) - 0.5 * mu;
		gradient[1] = 0.25 * mu;
		gradient[2] = 0.0;
	}

	void InvertibleStVKMaterial::getEnergyHassian(const double *invariants, double *hassian) const {
		Initiation(hassian, 6);
		hassian[0] = 0.25 * lambda;
	}
}