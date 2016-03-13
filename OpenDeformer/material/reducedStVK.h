#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_MATERIAL_STVK_H
#define ODER_MATERIAL_STVK_H

#include "reducedHyperelastic.h"

namespace ODER{
	class ReducedStVKMaterial :public ReducedHyperelasticMaterial{
	public:
		ReducedStVKMaterial(double rho, double lameFirst, double lameSecond, int orders);
		void generateStiffnessMatrix(const Reference<Mesh> &mesh, const Reference<NodeIndexer> &indexer, SparseMatrixAssembler& matrix) const;
		void preprocessWithReduction(const Reference<Mesh> &mesh, const Reference<NodeIndexer> &indexer);
		void getNodeForces(const Reference<Mesh> &mesh, const Reference<NodeIndexer> &indexer, int order, int totalDofs, const double *ds, double *forces);
		~ReducedStVKMaterial();
	private:
		void getNodeDisplacements(const double *ds, const int *nodeIndices, VectorBase<double>& d) const;
		double lambda, mu;
		double D[3];
		double *intergration[2];
		double *stressNonlinear;
	};
}

#endif
