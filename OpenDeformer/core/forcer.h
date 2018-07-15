#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_CORE_FORCER_H
#define ODER_CORE_FORCER_H

#include "oder.h"
#include "memory.h"
#include "mesh.h"
#include <map>

namespace ODER{
	class Forcer{
	public:
		void addBodyForce(const Reference<Mesh> &mesh, Scalar bodyForce[3], const Reference<NodeIndexer> &indexer);
		void addSurfaceForceByNode(const Reference<Mesh> &mesh, Scalar nodeForce[3], int nodeCounts, int *nodeIndex, const Reference<NodeIndexer> &indexer);
		void addSurfaceForceByElement(const Reference<Mesh> &mesh, Scalar elementForce[3], int surfaceCounts, const int *surfaceIndex, const Reference<NodeIndexer> &indexer);
		void getVirtualWorks(int dofs, int totalDofs, const Scalar *basises, Scalar *vws) const;
	private:
		void addConstraint(int constraintCounts);
		std::map<int, Scalar> virtualWorks;
	};

	class StaDampingForcer {
	public:
		virtual void addDampingMatrix(const SparseSymMatrixIndicesPerElementCache& matrixIndices, BlockedSymSpMatrix& mat) const = 0;
		virtual ~StaDampingForcer() = default;
	};

	template<class SpMatrix> class DynDampingForcer {
	public:
		virtual void addDampingMatrix(AggSpMatrix<SpMatrix>& mat) const = 0;
		virtual ~DynDampingForcer() = default;
	};
}

#endif