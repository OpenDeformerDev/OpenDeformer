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
		Intergrator(int DOFS, Scalar massDamp, Scalar stiffDamp, Scalar ts)
			:dofs(DOFS), massDamping(massDamp), stiffnessDamping(stiffDamp), timeStep(ts){}

		virtual void setExternalVirtualWork(const Forcer& forcer) = 0;
		virtual void addExternalDampingForceMatrix(const DampingForcer& forcer) {
			Severe("Unimplemented features in Intergrator::setExternalDampingForceMatrix");
		}
		virtual void clearExternalDampingForceMatrix() {
			Severe("Unimplemented features in Intergrator::clearExternalDampingForceMatrix");
		}
		virtual void runOneTimeStep() = 0;
		virtual void updateMeshVerticesDisplacements(const Reference<NodeIndexer> &indexer, Reference<Mesh> &mesh) const = 0;
		virtual void getMeshVerticesVelocities(const Reference<NodeIndexer> &indexer, const Reference<Mesh> &mesh, Vector3 *velocities) const = 0;
		virtual ~Intergrator() = default;

	protected:
		int dofs;
		Scalar massDamping;
		Scalar stiffnessDamping;
		Scalar timeStep;
	};
}

#endif