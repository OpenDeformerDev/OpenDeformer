#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_MATERIAL_INVERTIBLEHYPERELASTIC_H
#define ODER_MATERIAL_INVERTIBLEHYPERELASTIC_H

#include "fullOrderNonlinear.h"

namespace ODER {
	class InvertibleHyperelasticMaterial : public FullOrderNonlinearMaterial<BlockedSymSpMatrix>{
	public:
		InvertibleHyperelasticMaterial(Scalar rho, Scalar inversionTrashold) :
			FullOrderNonlinearMaterial(rho, MarterialType(Marterial_Isotropic | Marterial_Invertible)), trashold(inversionTrashold) {}
		void generateMatrixAndVirtualWorks(const Reference<Mesh> &mesh, const Reference<NodeIndexer> &indexer,
			const Scalar *precomputes, const int *matrixIndices, BlockedSymSpMatrix& matrix, Scalar *vws) const;
		Scalar *getPrecomputes(const Reference<Mesh> &mesh) const;
		virtual ~InvertibleHyperelasticMaterial() = default;
	private:
		virtual void getEnergyGradient(const Scalar *invariants, Scalar *gradient) const = 0;
		virtual void getEnergyHassian(const Scalar *invariants, Scalar *hassian) const = 0;
		void getInvariants(const Scalar *singularValues, Scalar *invriants) const {
			Scalar singularVal02 = singularValues[0] * singularValues[0];
			Scalar singularVal12 = singularValues[1] * singularValues[1];
			Scalar singularVal22 = singularValues[2] * singularValues[2];

			invriants[0] = singularVal02 + singularVal12 + singularVal22;
			invriants[1] = singularVal02 * singularVal02 + singularVal12 * singularVal12 + singularVal22 * singularVal22;
			invriants[2] = singularVal02 * singularVal12 * singularVal22;
		}
		void modifiedDeformGradient(const Scalar *gradient, Scalar *diags, Scalar *leftOrthoMat, Scalar *rightOrthoMat) const;
		void getPiolaKirchhoffStress(const Scalar *diags, const Scalar *leftOrthoMat, const Scalar *rightOrthoMat, 
			const Scalar *invariants, const Scalar *eneryGradient, Scalar *stress) const;

		Scalar trashold;
	}; 
}

#endif


