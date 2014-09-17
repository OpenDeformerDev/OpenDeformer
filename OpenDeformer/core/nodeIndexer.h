#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_CORE_NODEINDEXER_H
#define ODER_CORE_NODEINDEXER_H

#include "oder.h"
#include "memory.h"

namespace ODER{
	class NodeIndexer{
	public:
		NodeIndexer(){}
		NodeIndexer(const Constrainer& constrainer);
		int getGlobalIndex(const Element& element, int localNodeIndex, int axis) const;
		int getGlobalIndex(const Facet& facet, int localVertIndex, int axis) const;
		int getGlobalIndex(int nodeIndex, int axis) const;
		int getConstrainIndices(const int** detector) const;
		int getMatrixOrder(const Reference<Mesh> &mesh) const;

	private:
		vector<int> constrainIndices;
	};
}

#endif