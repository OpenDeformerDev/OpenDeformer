#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_MATERIAL_HOOKE_H
#define ODER_MATERIAL_HOOKE_H

#include "mechMaterial.h"
#include "mesh.h"
#include "element.h"

namespace ODER{
	class HookeMaterial :public MechMaterial{
	public:
		HookeMaterial(double rho, double YOUNGS, double POSSION);
		void generateStiffnessMatrix(const Reference<Mesh> &mesh, const Reference<NodeIndexer> &indexer, SparseMatrixAssembler& matrix) const;
	private:
		double youngs, possion;
		double D[3];
	};
}

#endif