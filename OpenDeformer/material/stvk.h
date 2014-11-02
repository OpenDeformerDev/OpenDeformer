#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_MATERIAL_STVK_H
#define ODER_MATERIAL_STVK_H

#include "mecMaterial.h"
#include "mesh.h"
#include "element.h"

namespace ODER{
	class StVKMaterial :public MecMaterial{
	public:
		StVKMaterial(double rho, double lameFirst, double lameSecond, int numElements, int orders);
		void generateStiffnessMatrix(const Reference<Mesh> &mesh, const Reference<NodeIndexer> &indexer, SparseMatrixAssembler& matrix) const;
		void preprocessWithReduction(const Reference<Mesh> &mesh, const Reference<NodeIndexer> &indexer);
		void getNodeForces(const Reference<Mesh> &mesh, const Reference<NodeIndexer> &indexer, int order, const double *ds, double *forces);
		int getNonlinearAsymptoticOrder(){ return orders; }
		~StVKMaterial();
	private:
		void getNodeDisplacements(const double *ds, const int *nodeIndices, VectorBase<double>& d) const;

		double lambda, mu;
		double D[3];
		int orders;
		double **intergration[2];
		double ***stressNonlinear;
	};
}

#endif
