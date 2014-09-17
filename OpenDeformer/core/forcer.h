#pragma once
#include "oder.h"
#include "latool.h"
#include "memory.h"
#include "mesh.h"
#include "element.h"

namespace ODER{
	class Forcer{
	public:
		Forcer::Forcer(const Reference<Mesh> &m) : mesh(m){}
		void addBodyForce(double bodyForce[3], const NodeIndexer &indexer);
		void addSurfaceForceByNode(double nodeForce[3], int nodeCounts, int *nodeIndex, const NodeIndexer &indexer);
		void addSurfaceForceByElement(double elementForce[3], int surfaceCounts, int *surfaceIndex, const NodeIndexer &indexer);
		void getVirtualWorks(int dofs, int totalDofs, const double *basises, double *vws) const;
	private:
		void addConstraint(int constraintCounts);
		Reference<Mesh> mesh;
		map<int, double> virtualWorks;
	};
}