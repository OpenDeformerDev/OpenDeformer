#include "stdafx.h"
#include "invertibleNeoHookean.h"

namespace ODER {
	InvertibleNeoHookeanMaterial::InvertibleNeoHookeanMaterial(double rho, double youngs, double possion,
		double inversionTrashold) : InvertibleHyperelasticMaterial(rho, inversionTrashold) {
		lambda = possion * youngs / ((1.0 + possion) * (1.0 - 2.0 * possion));
		mu = youngs / (2.0 * (1.0 + possion));
	}

	void InvertibleNeoHookeanMaterial::getEnergyGradient(const double *invariants, double *gradient) const {
		gradient[0] = 0.5 * mu;
		gradient[1] = 0.0;
		gradient[2] = (-0.5 * mu + 0.25 * lambda * log(invariants[2])) / invariants[2];
	}

	void InvertibleNeoHookeanMaterial::getEnergyHassian(const double *invariants, double *hassian) const {
		Initiation(hassian, 6);
		hassian[5] = (0.5 * mu + 0.25 * lambda * (1.0 - log(invariants[2]))) / (invariants[2] * invariants[2]);
	}
}