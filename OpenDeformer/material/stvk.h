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
		StVKMaterial(double rho, double lameFirst, double lameSecond, const Reference<Mesh> &mesh, int orders);
		void generateStiffnessMatrix(const Reference<Mesh> &mesh, const NodeIndexer &indexer, SparseMatrixAssembler& matrix) const;
		void preprocessWithReduction(const Reference<Mesh> &mesh, const NodeIndexer &indexer, int dofs, const double *basises);
		void getNodeForces(double *ds, double *forces) const;
		~StVKMaterial();
	private:
		double lambda, mu;
		double D[3];
		int orders;
		double *intergration[2];
		double *loadParameter;
	};
}

#endif
