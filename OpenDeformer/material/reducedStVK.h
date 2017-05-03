#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_MATERIAL_STVK_H
#define ODER_MATERIAL_STVK_H

#include "reducedHyperelastic.h"

namespace ODER{
	class ReducedStVKMaterial :public ReducedHyperelasticMaterial{
	public:
		ReducedStVKMaterial(Scalar rho, Scalar lameFirst, Scalar lameSecond, int orders);
		void generateStiffnessMatrix(const Reference<Mesh> &mesh, const Reference<NodeIndexer> &indexer, SparseMatrixAssembler& matrix) const;
		void preprocessWithReduction(const Reference<Mesh> &mesh, const Reference<NodeIndexer> &indexer);
		void getNodeForces(const Reference<Mesh> &mesh, const Reference<NodeIndexer> &indexer, int order, int totalDofs, const Scalar *ds, Scalar *forces);
		~ReducedStVKMaterial();
	private:
		void getNodeDisplacements(const Scalar *ds, const int *nodeIndices, VectorBase<Scalar>& d) const;
		Scalar lambda, mu;
		Scalar D[3];
		Scalar *intergration[2];
		Scalar *stressNonlinear;
	};
}

#endif
