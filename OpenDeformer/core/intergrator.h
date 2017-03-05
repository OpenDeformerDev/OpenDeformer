#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_CORE_INTERGRATOR_H
#define ODER_CORE_INTERGRATOR_H

#include "oder.h"
#include "mesh.h"
#include "nodeIndexer.h"

namespace ODER{
	class Intergrator{
	public:
		Intergrator(int DOFS, double massDamp, double stiffDamp, double ts)
			:dofs(DOFS), massDamping(massDamp), stiffnessDamping(stiffDamp), timeStep(ts){}

		virtual void setExternalVirtualWork(const Forcer& forcer) = 0;
		virtual void runOneTimeStep() = 0;
		virtual void updateMeshVerticesDisplacements(const Reference<NodeIndexer> &indexer, Reference<Mesh> &mesh) const = 0;
		virtual ~Intergrator() = default;

	protected:
		int dofs;
		double massDamping;
		double stiffnessDamping;
		double timeStep;
	};
}

#endif