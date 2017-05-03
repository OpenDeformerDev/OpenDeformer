#include "stdafx.h"
#include "invertibleStVK.h"

namespace ODER {
	InvertibleStVKMaterial::InvertibleStVKMaterial(Scalar rho, Scalar youngs, Scalar possion,
		Scalar inversionTrashold) : InvertibleHyperelasticMaterial(rho, inversionTrashold) {
		lambda = possion * youngs / ((Scalar(1.0) + possion) * (Scalar(1.0) - Scalar(2.0) * possion));
		mu = youngs / (Scalar(2.0) * (Scalar(1.0) + possion));
	}

	void InvertibleStVKMaterial::getEnergyGradient(const Scalar *invariants, Scalar *gradient) const {
		gradient[0] = Scalar(0.25) * lambda * (invariants[0] - Scalar(3.0)) - Scalar(0.5) * mu;
		gradient[1] = Scalar(0.25) * mu;
		gradient[2] = Scalar(0);
	}

	void InvertibleStVKMaterial::getEnergyHassian(const Scalar *invariants, Scalar *hassian) const {
		Initiation(hassian, 6);
		hassian[0] = Scalar(0.25) * lambda;
	}
}