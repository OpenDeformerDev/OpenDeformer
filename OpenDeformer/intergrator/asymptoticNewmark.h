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
		void setFullDisplacement(int order);
		void setReducedVirtualWorks(int order);
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