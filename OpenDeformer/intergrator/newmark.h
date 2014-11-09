#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_INTERGRATOR_NEWMARK_H
#define ODER_INTERGRATOR_NEWMARK_H

#include "intergrator.h"
#include "eigenslover.h"

namespace ODER{
	class LinearNewmark : public Intergrator{
	public:
		LinearNewmark(int totalDOFS, double beta, double gamma, int DOFS, double massDamp, double stiffDamp, double ts,
			const Reference<Mesh> m, const Reference<NodeIndexer>& nodeIndexer, const MecMaterial& mater);
		void setExternalVirtualWork(const Forcer& forcer);
		void runOneTimeStep();
		void getRawDisplacements(double *displacements) const;
		~LinearNewmark();
	private:
		int totalDofs;
		double betaDeltaT2, gammaDeltaT, minusBetaDeltaT2, minusGammaDeltaT;
		double *frequencies2;
		double *basises;
	};
}

#endif