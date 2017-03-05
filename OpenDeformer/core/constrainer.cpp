#include "stdafx.h"
#include "constrainer.h"
#include "element.h"

namespace ODER{
	void Constrainer::setConstrain(int nodeIndex, ConstrainType constrainType){
		switch (constrainType){
		case Constrain_X:
			constrains.insert(nodeIndex * 3);
			break;
		case Constrain_Y:
			constrains.insert(nodeIndex * 3 + 1);
			break;
		case Constrain_Z:
			constrains.insert(nodeIndex * 3 + 2);
			break;
		case Constrain_XY:
			constrains.insert(nodeIndex * 3);
			constrains.insert(nodeIndex * 3 + 1);
			break;
		case Constrain_YZ:
			constrains.insert(nodeIndex * 3 + 1);
			constrains.insert(nodeIndex * 3 + 2);
			break;
		case Constrain_XZ:
			constrains.insert(nodeIndex * 3);
			constrains.insert(nodeIndex * 3 + 2);
			break;
		case Constrain_XYZ:
			constrains.insert(nodeIndex * 3);
			constrains.insert(nodeIndex * 3 + 1);
			constrains.insert(nodeIndex * 3 + 2);
			break;
		default:
			//Warning("Unexpected mode to constrain node, constain all coordinate"); 
			constrains.insert(nodeIndex * 3);
			constrains.insert(nodeIndex * 3 + 1);
			constrains.insert(nodeIndex * 3 + 2);
			break;
		}
	}

	void Constrainer::setConstrains(int nodeCounts, int *nodeIndices, ConstrainType *constrainTypes){
		if (!constrainTypes){
			for (int i = 0; i < nodeCounts; i++){
				constrains.insert(nodeIndices[i] * 3);
				constrains.insert(nodeIndices[i] * 3 + 1);
				constrains.insert(nodeIndices[i] * 3 + 2);
			}
		}
		else{
			for (int i = 0; i < nodeCounts; i++)
				setConstrain(nodeIndices[i], constrainTypes[i]);
		}
	}

	void Constrainer::setConstrains(const Vector3d& origin, const Vector3d& normal, ConstrainType constrainType){
		const double epsilon = 1e-4;
		Facet *facet = mesh->getFacet();

		for (int i = 0; i < mesh->getFacetCount(); i++){
			facet->setVertIndexs(i);
			for (int j = 0; j < mesh->getVertPerFacetCount(); j++){
				int index = facet->getVertIndex(j);
				Vector3d vert = mesh->getVertex(index);
				if (fabs((vert - origin)*normal) < epsilon){
					setConstrain(index, constrainType);
				}
			}
		}
		delete facet;
	}

	void Constrainer::outputConstrains(std::vector<int> &out) const{
		int size = constrains.size();
		out.reserve(size);
		for (auto constrain : constrains)
			out.push_back(constrain);
	}
}