#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_MATERIAL_INVERTIBLEHYPERELASTIC_H
#define ODER_MATERIAL_INVERTIBLEHYPERELASTIC_H

#include "fullOrderNonlinear.h"

namespace ODER {
	class InvertibleHyperelasticMaterial : public FullOrderNonlinearMaterial<BlockedSymSpMatrix>{
	public:
		InvertibleHyperelasticMaterial(double rho, double inversionTrashold) :
			FullOrderNonlinearMaterial(rho, MarterialType(Marterial_Isotropic | Marterial_Invertible)), trashold(inversionTrashold) {}
		void generateMatrixAndVirtualWorks(const Reference<Mesh> &mesh, const Reference<NodeIndexer> &indexer,
			const double *precomputes, const int *matrixIndices, BlockedSymSpMatrix& matrix, double *vws) const;
		double *getPrecomputes(const Reference<Mesh> &mesh) const;
		virtual ~InvertibleHyperelasticMaterial() = default;
	private:
		virtual void getEnergyGradient(const double *invariants, double *gradient) const = 0;
		virtual void getEnergyHassian(const double *invariants, double *hassian) const = 0;
		void getInvariants(const double *singularValues, double *invriants) const {
			double singularVal02 = singularValues[0] * singularValues[0];
			double singularVal12 = singularValues[1] * singularValues[1];
			double singularVal22 = singularValues[2] * singularValues[2];

			invriants[0] = singularVal02 + singularVal12 + singularVal22;
			invriants[1] = singularVal02 * singularVal02 + singularVal12 * singularVal12 + singularVal22 * singularVal22;
			invriants[2] = singularVal02 * singularVal12 * singularVal22;
		}
		void modifiedDeformGradient(const double *gradient, double *diags, double *leftOrthoMat, double *rightOrthoMat) const;
		void getPiolaKirchhoffStress(const double *diags, const double *leftOrthoMat, const double *rightOrthoMat, 
			const double *invariants, const double *eneryGradient, double *stress) const;

		double trashold;
	}; 
}

#endif


