#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_CORE_NODEINDEXER_H
#define ODER_CORE_NODEINDEXER_H

#include "oder.h"
#include "memory.h"
#include <vector>

namespace ODER{
	class NodeIndexer : public ReferenceCounted {
	public:
		typedef std::vector<int>::const_iterator ConstrainIter;

		NodeIndexer(const Constrainer& constrainer);
		int getGlobalIndex(const Element& element, int localNodeIndex, int axis) const;
		int getGlobalIndex(const Facet& facet, int localVertIndex, int axis) const;
		int getGlobalIndex(int nodeIndex, int axis) const;
		void getElementNodesGlobalIndices(const Element& element, int nodeCounts, int *nodeIndices) const;
		void getElementNodesGlobalIndices(const Reference<Mesh> &mesh, int elementIndex, int nodeCounts, int *nodeIndices) const;
		int getMatrixOrder(const Reference<Mesh> &mesh) const;
		ConstrainIter getConstrainIterBegin() const{ return constrainIndices.begin(); }
		ConstrainIter getConstrainIterEnd() const{ return constrainIndices.end(); }

	private:
		std::vector<int> constrainIndices;
	};
}

#endif