#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_MATERIAL_COROTATIONALHYPERELASTIC_H
#define ODER_MATERIAL_COROTATIONALHYPERELASTIC_H

#include "fullOrderNonlinear.h"

namespace ODER {
	class CorotationalHyperelasticMaterial : public FullOrderNonlinearMaterial<BlockedSymSpMatrix> {
	public :
		CorotationalHyperelasticMaterial(Scalar rho, Scalar youngs, Scalar possion, Scalar compressingTreshold);
		FullOrderNonlinearMaterialCache getPrecomputes(const Reference<Mesh> &mesh) const;
		void generateMatrixAndVirtualWorks(const Reference<Mesh> &mesh, const Reference<NodeIndexer> &indexer,
			const FullOrderNonlinearMaterialCache& precomputes, const SparseSymMatrixIndicesPerElementCache& matrixIndices, BlockedSymSpMatrix& matrix, Scalar *vws) const;
		virtual ~CorotationalHyperelasticMaterial() = default;
	private:
		Scalar threshold;
		Scalar D[3];
	};
}

#endif
