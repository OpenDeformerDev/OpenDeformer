#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_MATERIAL_INVERTIBLEHYPERELASTIC_H
#define ODER_MATERIAL_INVERTIBLEHYPERELASTIC_H

#include "mechMaterial.h"
#include "element.h"

namespace ODER {
	class InvertibleHyperelasticMaterial : public MechMaterial {
		InvertibleHyperelasticMaterial(double rho, double inversionTrashold, const Reference<Mesh> &mesh);
		void generateTagentStiffnessMatrix(const Reference<Mesh> &mesh, const Reference<NodeIndexer> &indexer, 
			const double *u, BlockedSymSpMatrix& matrix) const;
		void generateInternalVirtualWorks(const Reference<Mesh> &mesh, const Reference<NodeIndexer> &indexer,
			const double *u, double *vws) const;
		void generateMatrixAndVirtualWorks(const Reference<Mesh> &mesh, const Reference<NodeIndexer> &indexer,
			const double *u, BlockedSymSpMatrix& matrix, double *vws);
		virtual ~InvertibleHyperelasticMaterial();
	private:
		virtual void getEnergyGradient(const double *invariants, double *gradient) const = 0;
		virtual void getEnergyHassian(const double *invariants, double *hassian) const = 0;
		void getInvariants(const double *singularValues, double invriants[3]) const {
			double singularVal02 = singularValues[0] * singularValues[0];
			double singularVal12 = singularValues[1] * singularValues[1];
			double singularVal22 = singularValues[2] * singularValues[2];

			invriants[0] = singularVal02 + singularVal12 + singularVal22;
			invriants[1] = singularVal02 * singularVal02 + singularVal12 * singularVal12 + singularVal22 * singularVal22;
			invriants[2] = singularVal02 * singularVal12 * singularVal22;
		}

		double *shapeFunctionDrivativesPrecomputed;
		double *deformationGradientPrecomputed;
		double *nodalVirtualWorkPrecomputed;

		InvertibleHyperelasticElement *element;
		double trashold;
	};
}

#endif


