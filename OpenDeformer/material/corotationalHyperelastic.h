#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_MATERIAL_COROTATIONALHYPERELASTIC_H
#define ODER_MATERIAL_COROTATIONALHYPERELASTIC_H

#include "fullOrderNonlinear.h"

namespace ODER {
	class CorotationalMaterial : public FullOrderNonlinearMaterial<BlockedSymSpMatrix> {
	public :
		CorotationalMaterial(Scalar rho, Scalar youngs, Scalar possion, Scalar compressingTreshold);
		Scalar *getPrecomputes(const Reference<Mesh> &mesh) const;
		void generateMatrixAndVirtualWorks(const Reference<Mesh> &mesh, const Reference<NodeIndexer> &indexer,
			const Scalar *precomputes, const int *matrixIndices, BlockedSymSpMatrix& matrix, Scalar *vws) const;
		virtual ~CorotationalMaterial() = default;
	private:
		Scalar threshold;
		Scalar D[3];
	};
}

#endif
