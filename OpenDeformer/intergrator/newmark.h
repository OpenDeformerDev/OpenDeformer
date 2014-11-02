#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_INTERGRATOR_NEWMARK_H
#define ODER_INTERGRATOR_NEWMARK_H

#include "oder.h"
#include "intergrator.h"
#include "eigenslover.h"

namespace ODER{
	class NewmarkIntergrator : public Intergrator{
	public:
		NewmarkIntergrator(int totalDOFS, double beta, double gamma, int DOFS, double massDamp, double stiffDamp, double ts, 
			const Reference<Mesh> m, const Reference<NodeIndexer>& nodeIndexer, const Reference<MecMaterial> &mater);
		void setExternalVirtualWork(const Forcer& forcer);
		void runOneTimeStep();
		void getDisplacements(int displacementCount, double *displacements) const;
		void getVertexPositions(Vector *vertices, double *displacementBuffer = NULL) const;
		~NewmarkIntergrator();
	private:
		void getRawDisplacements(double *displacements) const;
		int totalDofs;
		double betaDeltaT2, gammaDeltaT, minusBetaDeltaT2, minusGammaDeltaT;
		double *frequencies2;
		double *basises;
		double *loadFactors;
	};
}

#endif