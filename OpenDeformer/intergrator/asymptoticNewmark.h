#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_INTERGRATOR_ASYMPTOTICNEWMARK_H
#define ODER_INTERGRATOR_ASYMPTOTICNEWMARK_H

#include "intergrator.h"
#include "eigenslover.h"
#include "hyperelastic.h"

namespace ODER{
	class AsymptoticNewmark : public Intergrator{
	public:
		AsymptoticNewmark(double beta, double gamma, int DOFS, double massDamp, double stiffDamp, double ts,
			const Reference<Mesh> m, const Reference<NodeIndexer>& nodeIndexer, HyperelasticMaterial* mater);
		void setExternalVirtualWork(const Forcer& forcer);
		void runOneTimeStep();
		void getRawDisplacements(double *displacements) const;
		~AsymptoticNewmark();

	private:
		void setFullDisplacement(int order, double loadFactor);
		void setReducedVirtualWorks(int order, const double *fullVirtualWork);

		double *d;
		double *v;
		double *a;
		double *pre_d;
		double *pre_v;
		double *pre_a;
		double *externalVirtualWork;

		Reference<Mesh> mesh;
		Reference<NodeIndexer> indexer;
		HyperelasticMaterial* material;
		int totalDofs;
		double betaDeltaT2, gammaDeltaT, minusBetaDeltaT2, minusGammaDeltaT;
		double *frequencies2;
		double *basises;
		double *loadFactors;
		double *fullDisplacements;
		double *vwBuffer;
	};
}

#endif