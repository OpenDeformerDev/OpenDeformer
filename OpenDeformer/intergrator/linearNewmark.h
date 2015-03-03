#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_INTERGRATOR_NEWMARK_H
#define ODER_INTERGRATOR_NEWMARK_H

#include "intergrator.h"
#include "eigenSolver.h"

namespace ODER{
	class LinearNewmark : public Intergrator{
	public:
		LinearNewmark(double beta, double gamma, int DOFS, double massDamp, double stiffDamp, double ts,
			const Reference<Mesh> m, const Reference<NodeIndexer>& nodeIndexer, const MechMaterial& mater);
		void setExternalVirtualWork(const Forcer& forcer);
		void runOneTimeStep();
		void getRawDisplacements(double *displacements) const;
		~LinearNewmark();
	private:
		double *d;
		double *v;
		double *a;
		double *pre_d;
		double *pre_v;
		double *pre_a;
		double *externalVirtualWork;

		int totalDofs;
		double betaDeltaT2, gammaDeltaT, minusBetaDeltaT2, minusGammaDeltaT;
		double *frequencies2;
		double *basises;
	};
}

#endif