#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_MATERIAL_INVERTIBLEHYPERELASTIC_H
#define ODER_MATERIAL_INVERTIBLEHYPERELASTIC_H

#include "mechMaterial.h"

namespace ODER {
	class InvertibleHyperelasticMaterial : public MechMaterial {
		InvertibleHyperelasticMaterial(double rho, double inversionTrashold) :
			MechMaterial(rho, MarterialType(Marterial_Isotropic | Marterial_NonLinear | Marterial_Invertible)), trashold(inversionTrashold){}
		void generateTagentStiffnessMatrix(const Reference<Mesh> &mesh, const Reference<NodeIndexer> &indexer, 
			const double *u, BlockedSymSpMatrix& matrix) const;
		void generateInternalVirtualWorks(const Reference<Mesh> &mesh, const Reference<NodeIndexer> &indexer,
			const double *u, double *vws) const;
		void generateMatrixAndVirtualWorks(const Reference<Mesh> &mesh, const Reference<NodeIndexer> &indexer,
			const double *u, BlockedSymSpMatrix& matrix, double *vws);
	private:
		double trashold;
	};
}

#endif


