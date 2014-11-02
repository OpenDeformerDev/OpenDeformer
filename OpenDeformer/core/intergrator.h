#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_CORE_INTERGRATOR_H
#define ODER_CORE_INTERGRATOR_H

#include "oder.h"
#include "memory.h"
#include "mecMaterial.h"
#include "nodeIndexer.h"

namespace ODER{
	class Intergrator{
	public:
		Intergrator(int DOFS, double massDamp, double stiffDamp, double ts, const Reference<Mesh> m, const Reference<NodeIndexer>& nodeIndexer, const Reference<MecMaterial> &mater);

		virtual void setExternalVirtualWork(const Forcer& forcer) = 0;
		virtual void runOneTimeStep() = 0;
		virtual void getDisplacements(int displacementCount, double *displacements) const = 0;
		virtual void getVertexPositions(Vector *vertices, double *displacementBuffer = NULL) const = 0;
		virtual ~Intergrator();

	protected:
		double *d;
		double *v;
		double *a;
		double *pre_d;
		double *pre_v;
		double *pre_a;

		double* externalVirtualWork;

		int dofs;
		double massDamping;
		double stiffnessDamping;
		double timeStep;
		Reference<Mesh> mesh;
		Reference<NodeIndexer> indexer;
		Reference<MecMaterial> material;
	};
}

#endif