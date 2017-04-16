#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_MATERIAL_INVERTIBLENEOHOOKE_H
#define ODER_MATERIAL_INVERTIBLENEOHOOKE_H

#include "invertibleHyperelastic.h"

namespace ODER {
	class InvertibleNeoHookeanMaterial : public InvertibleHyperelasticMaterial {
	public:
		InvertibleNeoHookeanMaterial(double rho, double yonugs, double possion, double inversionTrashold);
	private:
		void getEnergyGradient(const double *invariants, double *gradient) const;
		void getEnergyHassian(const double *invariants, double *hassian) const;

		double lambda, mu;
	};
}


#endif
