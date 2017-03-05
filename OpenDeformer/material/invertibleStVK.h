#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_MATERIAL_INVERTIBLESTVK_H
#define ODER_MATERIAL_INVERTIBLESTVK_H

#include "invertibleHyperelastic.h"

namespace ODER {
	class InvertibleStVKMaterial : public InvertibleHyperelasticMaterial {
	public:
		InvertibleStVKMaterial(double rho, double yonugs, double possion, double inversionTrashold);
	private:
		void getEnergyGradient(const double *invariants, double *gradient) const;
		void getEnergyHassian(const double *invariants, double *hassian) const;

		double lambda, mu;
	};
}

#endif
