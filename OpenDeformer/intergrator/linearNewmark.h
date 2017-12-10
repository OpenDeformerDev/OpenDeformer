#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_INTERGRATOR_NEWMARK_H
#define ODER_INTERGRATOR_NEWMARK_H

#include "intergrator.h"
#include "eigenSolver.h"
#include "hooke.h"

namespace ODER{
	class LinearNewmark : public Intergrator{
	public:
		LinearNewmark(Scalar beta, Scalar gamma, int DOFS, Scalar massDamp, Scalar stiffDamp, Scalar ts,
			const Reference<Mesh>& m, const Reference<NodeIndexer>& nodeIndexer, const HookeMaterial& mater);
		void setExternalVirtualWork(const Forcer& forcer);
		void runOneTimeStep();
		void updateMeshVerticesDisplacements(const Reference<NodeIndexer> &indexer, Reference<Mesh> &mesh) const;
		void getMeshVerticesVelocities(const Reference<NodeIndexer> &indexer, const Reference<Mesh> &mesh, Vector3 *velocities) const;
		~LinearNewmark();
	private:
		void getFullOrderProperties(const Scalar *reducedProperties, Scalar *fullOrderProperties) const;

		Scalar *d;
		Scalar *v;
		Scalar *a;
		Scalar *pre_d;
		Scalar *pre_v;
		Scalar *pre_a;
		Scalar *externalVirtualWork;

		int totalDofs;
		Scalar betaDeltaT2, gammaDeltaT, minusBetaDeltaT2, minusGammaDeltaT;
		Scalar *frequencies2;
		Scalar *basises;
	};
}

#endif