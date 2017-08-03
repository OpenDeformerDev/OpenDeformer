#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_MATERIAL_COROTATIONALPLASTIC_H
#define ODER_MATERIAL_COROTATIONALPLASTIC_H

#include "fullOrderNonlinear.h"

namespace ODER {
	class CorotationalPlasticMaterial : public FullOrderNonlinearMaterial<BlockedSymSpMatrix> {
	public:
		CorotationalPlasticMaterial(Scalar rho, Scalar youngs, Scalar possion, 
			Scalar yieldCriterion, Scalar creepRatio, Scalar maxPlasticStrain,
			Scalar compressingTreshold);
		FullOrderNonlinearMaterialCache getPrecomputes(const Reference<Mesh> &mesh) const;
		void generateMatrixAndVirtualWorks(const Reference<Mesh> &mesh, const Reference<NodeIndexer> &indexer,
			const FullOrderNonlinearMaterialCache& cache, const SparseSymMatrixIndicesPerElementCache& matrixIndices, BlockedSymSpMatrix& matrix, Scalar *vws) const;
		virtual ~CorotationalPlasticMaterial() = default;
	private:
		void computeElasticStress(const Scalar *orthoMat, const Scalar *factoredDeformationGradient, Scalar* plasticStrain, Scalar *elasticStress) const;

		Scalar threshold;
		Scalar D[3];
		Scalar yieldCrierion;
		Scalar creepRatio;
		Scalar maxPlasticStrain;
	};
}

#endif
