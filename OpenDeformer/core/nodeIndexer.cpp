#include "stdafx.h"
#include "nodeIndexer.h"
#include "constrainer.h"
#include "mesh.h"
#include "element.h"

namespace ODER{
	NodeIndexer::NodeIndexer(const Constrainer& constrainer){
		constrainer.outputConstrains(constrainIndices);
	}

	int NodeIndexer::getGlobalIndex(const Element& element, int localNodeIndex, int axis) const{
		int index = 3 * element.getNodeIndex(localNodeIndex) + axis;
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
		int index = 3 * facet.getVertIndex(localVertIndex)+ axis;
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

	void NodeIndexer::getElementNodesGlobalIndices(const Reference<Mesh> &mesh, int elementIndex, int nodeCounts, int *nodeIndices) const {
		const int *indices = mesh->getElementNodeReference(elementIndex);
		for (int nodeIndex = 0; nodeIndex < nodeCounts; nodeIndex++) {
			for (int axis = 0; axis < 3; axis++)
				nodeIndices[nodeIndex * 3 + axis] = getGlobalIndex(indices[nodeIndex], axis);
		}
	}

	int NodeIndexer::getMatrixOrder(const Reference<Mesh> &mesh) const{
		return 3 * mesh->getNodeCount() - constrainIndices.size();
	}
}