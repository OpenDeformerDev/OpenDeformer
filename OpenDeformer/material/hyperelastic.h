#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_MATERIAL_HYPERELASTIC_H
#define ODER_MATERIAL_HYPERELASTIC_H

#include "mechMaterial.h"

namespace ODER{
	class HyperelasticMaterial :public MechMaterial{
	public:
		HyperelasticMaterial(double rho, MarterialType t, int numOrder) :MechMaterial(rho, t), orderCount(std::max(numOrder, 1)){}
		virtual void getNodeForces(const Reference<Mesh> &mesh, const Reference<NodeIndexer> &indexer, int order, int totalDofs, const double *ds, double *forces) = 0;
		virtual void preprocessWithReduction(const Reference<Mesh> &mesh, const Reference<NodeIndexer> &indexer) = 0;
		int getNonlinearAsymptoticOrder(){ return orderCount; }

	private:
		int orderCount;
	};
}


#endif