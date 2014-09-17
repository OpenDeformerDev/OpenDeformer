#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_CORE_INTERGRATOR_H
#define ODER_CORE_INTERGRATOR_H

#include "oder.h"
#include "memory.h"
#include "mecMaterial.h"
#include "forcer.h"

namespace ODER{
	class Intergrator{
	public:
		Intergrator(int DOFS, double massDamp, double stiffDamp, double ts);

		virtual void setExternalVirtualWork(const Forcer& forcer) = 0;
		virtual void runOneTimeStep() = 0;
		virtual void getDisplacements(const NodeIndexer &indexer, int displacementCount, double *displacements) const = 0;
		virtual void getVertexPositions(const Reference<Mesh> &mesh, const NodeIndexer &indexer, int verticesCount, Vector *vertices, double *displacementBuffer = NULL) const = 0;
		~Intergrator();

		int dofs;
		double massDamping;
		double stiffnessDamping;
		double timeStep;

	protected:
		double *d;
		double *v;
		double *a;
		double *pre_d;
		double *pre_v;
		double *pre_a;

		double* externalVirtualWork;
	};
}

#endif