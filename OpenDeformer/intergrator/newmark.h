#pragma once
#include "oder.h"
#include "intergrator.h"
#include "eigenslover.h"

namespace ODER{
	class NewmarkIntergrator : public Intergrator{
	public:
		NewmarkIntergrator(int totalDOFS, double beta, double gamma, const SparseMatrix &M, const SparseMatrix &K, int DOFS, double massDamp, double stiffDamp, double ts);
		void setExternalVirtualWork(const Forcer& forcer);
		void runOneTimeStep();
		void getDisplacements(const NodeIndexer &indexer, int displacementCount, double *displacements) const;
		void getVertexPositions(const Reference<Mesh> &mesh, const NodeIndexer &indexer, int verticesCount, Vector *vertices, double *displacementBuffer = NULL) const;
		~NewmarkIntergrator();
	private:
		void getRawDisplacements(double *displacements) const;
		int totalDofs;
		double betaDeltaT2, gammaDeltaT, minusBetaDeltaT2, minusGammaDeltaT;
		double *frequencies2;
		double *basises;
	};
}