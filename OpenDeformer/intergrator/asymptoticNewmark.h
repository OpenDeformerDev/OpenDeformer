#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_INTERGRATOR_ASYMPTOTICNEWMARK_H
#define ODER_INTERGRATOR_ASYMPTOTICNEWMARK_H

#include "intergrator.h"
#include "eigenSolver.h"
#include "reducedHyperelastic.h"

namespace ODER{
	class AsymptoticNewmark : public Intergrator{
	public:
		AsymptoticNewmark(Scalar beta, Scalar gamma, int DOFS, Scalar massDamp, Scalar stiffDamp, Scalar ts,
			const Reference<Mesh>& m, const Reference<NodeIndexer>& nodeIndexer, ReducedHyperelasticMaterial* mater);
		void setExternalVirtualWork(const Forcer& forcer);
		void runOneTimeStep();
		void updateMeshVerticesDisplacements(const Reference<NodeIndexer> &indexer, Reference<Mesh> &mesh) const;
		void getMeshVerticesVelocities(const Reference<NodeIndexer> &indexer, const Reference<Mesh> &mesh, Vector3 *velocities) const {
			Severe("Unimplemented features in AsymptoticNewmark::getMeshVerticesVelocities");
		}
		~AsymptoticNewmark();

	private:
		void getRawDisplacements(Scalar *displacements) const;

		Scalar *d;
		Scalar *v;
		Scalar *a;
		Scalar *pre_d;
		Scalar *pre_v;
		Scalar *pre_a;
		Scalar *externalVirtualWork;

		Reference<Mesh> mesh;
		Reference<NodeIndexer> indexer;
		ReducedHyperelasticMaterial* material;
		int totalDofs;
		Scalar betaDeltaT2, gammaDeltaT, minusBetaDeltaT2, minusGammaDeltaT;
		Scalar *frequencies2;
		Scalar *basises;
		Scalar *loadFactors;
		Scalar *fullDisplacements;
		Scalar *vwBuffer;
	};
}

#endif