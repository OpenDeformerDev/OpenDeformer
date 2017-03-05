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
		void addBodyForce(const Reference<Mesh> &mesh, double bodyForce[3], const Reference<NodeIndexer> &indexer);
		void addSurfaceForceByNode(const Reference<Mesh> &mesh, double nodeForce[3], int nodeCounts, int *nodeIndex, const Reference<NodeIndexer> &indexer);
		void addSurfaceForceByElement(const Reference<Mesh> &mesh, double elementForce[3], int surfaceCounts, const int *surfaceIndex, const Reference<NodeIndexer> &indexer);
		void getVirtualWorks(int dofs, int totalDofs, const double *basises, double *vws) const;
	private:
		void addConstraint(int constraintCounts);
		std::map<int, double> virtualWorks;
	};
}

#endif