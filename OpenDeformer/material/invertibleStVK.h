#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_MATERIAL_INVERTIBLESTVK_H
#define ODER_MATERIAL_INVERTIBLESTVK_H

#include "invertibleHyperelastic.h"

namespace ODER {
	class InvertibleStVKMaterial : public InvertibleHyperelasticMaterial {
	public:
		InvertibleStVKMaterial(Scalar rho, Scalar yonugs, Scalar possion, Scalar inversionTrashold);
	private:
		void getEnergyGradient(const Scalar *invariants, Scalar *gradient) const;
		void getEnergyHassian(const Scalar *invariants, Scalar *hassian) const;

		Scalar lambda, mu;
	};
}

#endif
