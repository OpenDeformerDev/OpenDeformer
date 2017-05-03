#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_MATERIAL_REDUCEDHYPERELASTIC_H
#define ODER_MATERIAL_REDUCEDHYPERELASTIC_H

#include "mechMaterial.h"

namespace ODER{
	class ReducedHyperelasticMaterial :public MechMaterial{
	public:
		ReducedHyperelasticMaterial(Scalar rho, int numOrder) :MechMaterial(rho, MarterialType(Marterial_Isotropic | Marterial_NonLinear | Marterial_Reduced)), orderCount(std::max(numOrder, 1)){}
		virtual void generateStiffnessMatrix(const Reference<Mesh> &mesh, const Reference<NodeIndexer> &indexer, SparseMatrixAssembler& matrix) const = 0;
		virtual void getNodeForces(const Reference<Mesh> &mesh, const Reference<NodeIndexer> &indexer, int order, int totalDofs, const Scalar *ds, Scalar *forces) = 0;
		virtual void preprocessWithReduction(const Reference<Mesh> &mesh, const Reference<NodeIndexer> &indexer) = 0;
		int getNonlinearAsymptoticOrder(){ return orderCount; }
		virtual ~ReducedHyperelasticMaterial() = default;
	private:
		int orderCount;
	};
}


#endif