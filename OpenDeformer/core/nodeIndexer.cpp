#include "stdafx.h"
#include "nodeIndexer.h"
#include "constrainer.h"
#include "mesh.h"

namespace ODER{
	NodeIndexer::NodeIndexer(const Constrainer& constrainer){
		constrainer.outputConstrains(constrainIndices);
	}

	int NodeIndexer::getGlobalIndex(const Element& element, int localNodeIndex, int axis) const{
		int index = 3 * element.nodeIndexs[localNodeIndex] + axis;
		if (constrainIndices.empty())
			return index;
		else{
			int offset = std::lower_bound(constrainIndices.begin(), constrainIndices.end(), index) - constrainIndices.begin();
			if (offset < constrainIndices.size() && index == constrainIndices[offset])
				return -1;
			else if (offset == 0)
				return index;
			else
				return index - offset;
		}
	}
	int NodeIndexer::getGlobalIndex(const Facet& facet, int localVertIndex, int axis) const{
		int index = 3 * facet.vertIndexs[localVertIndex] + axis;
		if (constrainIndices.empty())
			return index;
		else{
			int offset = std::lower_bound(constrainIndices.begin(), constrainIndices.end(), index) - constrainIndices.begin();
			if (offset < constrainIndices.size() && index == constrainIndices[offset])
				return -1;
			else if (offset == 0)
				return index;
			else
				return index - offset;
		}
	}


	int NodeIndexer::getGlobalIndex(int nodeIndex, int axis) const{
		int index = 3 * nodeIndex + axis;
		if (constrainIndices.empty())
			return index;
		else{
			int offset = std::lower_bound(constrainIndices.begin(), constrainIndices.end(), index) - constrainIndices.begin();
			if (offset < constrainIndices.size() && index == constrainIndices[offset])
				return -1;
			else if (offset == 0)
				return index;
			else
				return index - offset;
		}
	}

	void NodeIndexer::getElementNodesGlobalIndices(const Element& element, int nodeCounts, int *nodeIndices) const{
		for (int nodeIndex = 0; nodeIndex < nodeCounts; nodeIndex++){
			for (int axis = 0; axis < 3; axis++)
				nodeIndices[nodeIndex * 3 + axis] = getGlobalIndex(element, nodeIndex, axis);
		}
	}

	int NodeIndexer::getConstrainIndices(const int **detector) const{
		*detector = &constrainIndices[0];
		return constrainIndices.size();
	}

	int NodeIndexer::getMatrixOrder(const Reference<Mesh> &mesh) const{
		return 3 * mesh->numNodes - constrainIndices.size();
	}
}