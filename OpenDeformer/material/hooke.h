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
		HookeMaterial(Scalar rho, Scalar YOUNGS, Scalar POSSION);
		void generateStiffnessMatrix(const Reference<Mesh> &mesh, const Reference<NodeIndexer> &indexer, SparseMatrixAssembler& matrix) const;
	private:
		Scalar youngs, possion;
		Scalar D[3];
	};
}

#endif