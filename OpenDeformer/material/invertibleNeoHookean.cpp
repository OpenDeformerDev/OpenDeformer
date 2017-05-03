#include "stdafx.h"
#include "invertibleNeoHookean.h"

namespace ODER {
	InvertibleNeoHookeanMaterial::InvertibleNeoHookeanMaterial(Scalar rho, Scalar youngs, Scalar possion,
		Scalar inversionTrashold) : InvertibleHyperelasticMaterial(rho, inversionTrashold) {
		lambda = possion * youngs / ((Scalar(1.0) + possion) * (Scalar(1.0) - Scalar(2.0) * possion));
		mu = youngs / (Scalar(2.0) * (Scalar(1.0) + possion));
	}

	void InvertibleNeoHookeanMaterial::getEnergyGradient(const Scalar *invariants, Scalar *gradient) const {
		gradient[0] = Scalar(0.5) * mu;
		gradient[1] = Scalar(0);
		gradient[2] = (Scalar(-0.5) * mu + Scalar(0.25) * lambda * log(invariants[2])) / invariants[2];
	}

	void InvertibleNeoHookeanMaterial::getEnergyHassian(const Scalar *invariants, Scalar *hassian) const {
		Initiation(hassian, 6);
		hassian[5] = (Scalar(0.5) * mu + Scalar(0.25) * lambda * (Scalar(1.0) - log(invariants[2]))) / (invariants[2] * invariants[2]);
	}
}