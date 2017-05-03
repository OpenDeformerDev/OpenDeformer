#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_MATERIAL_INVERTIBLENEOHOOKE_H
#define ODER_MATERIAL_INVERTIBLENEOHOOKE_H

#include "invertibleHyperelastic.h"

namespace ODER {
	class InvertibleNeoHookeanMaterial : public InvertibleHyperelasticMaterial {
	public:
		InvertibleNeoHookeanMaterial(Scalar rho, Scalar yonugs, Scalar possion, Scalar inversionTrashold);
	private:
		void getEnergyGradient(const Scalar *invariants, Scalar *gradient) const;
		void getEnergyHassian(const Scalar *invariants, Scalar *hassian) const;

		Scalar lambda, mu;
	};
}


#endif
