#pragma once
#include "mecMaterial.h"
#include "mesh.h"
#include "element.h"

namespace ODER{
	class HookeMaterial :public MecMaterial{
	public:
		HookeMaterial(double rho, double YOUNGS, double POSSION);
		void generateStiffnessMatrix(const Reference<Mesh> &mesh, const NodeIndexer &indexer, SparseMatrixAssembler& matrix) const;
	private:
		double youngs, possion;
		double D[3];
	};
}