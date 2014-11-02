#include "stdafx.h"
#include "forcer.h"
#include "nodeIndexer.h"

namespace ODER{
	void Forcer::addBodyForce(const Reference<Mesh> &mesh, double bodyForce[3], const Reference<NodeIndexer> &indexer){
		Element *element = mesh->getEmptyElement();
		const int numNodesPerElement = mesh->numNodesPerElement;
		double *result = new double[3 * numNodesPerElement];

		for (int i = 0; i < mesh->numElements; i++){
			element->setNodeIndexs(i);
			element->getBodyVirtualWorks(bodyForce, result);
			for (int j = 0; j < numNodesPerElement; j++){
				for (int k = 0; k < 3; k++){
					if (result[3 * j + k] != 0.0){
						int nodeIndexOffset = indexer->getGlobalIndex(*element, j, k);
						if (nodeIndexOffset >= 0){
							auto found = virtualWorks.find(nodeIndexOffset);
							if (found != virtualWorks.end())
								found->second += result[3 * j + k];
							else
								virtualWorks.insert(std::pair<int, double>(nodeIndexOffset, result[3 * j + k]));
						}
					}
				}
			}
		}
		delete element;
		delete[] result;
	}

	void Forcer::addSurfaceForceByNode(const Reference<Mesh> &mesh, double nodeForce[3], int nodeCounts, int *nodeIndex, const Reference<NodeIndexer> &indexer){
		for (int i = 0; i < nodeCounts; i++){
			for (int j = 0; j < 3; j++){
				if (nodeForce[j] != 0.0){
					int nodeIndexOffset = indexer->getGlobalIndex(nodeIndex[i], j);
					if (nodeIndexOffset >= 0){
						auto found = virtualWorks.find(nodeIndexOffset);
						if (found != virtualWorks.end())
							found->second += nodeForce[j];
						else
							virtualWorks.insert(std::pair<int, double>(nodeIndexOffset, nodeForce[j]));
					}
				}
			}
		}
	}

	void Forcer::addSurfaceForceByElement(const Reference<Mesh> &mesh, double surfaceForce[3], int surfaceCounts, int *surfaceIndex, const Reference<NodeIndexer> &indexer){
		Facet *facet = mesh->getEmptyFacet();
		const int numVertPerSur = mesh->numVertPerSur;
		double *result = new double[3 * numVertPerSur];

		for (int i = 0; i < surfaceCounts; i++){
			facet->setVertIndexs(surfaceIndex[i]);
			facet->getSurfVirtualWorks(surfaceForce, result);
			for (int j = 0; j < numVertPerSur; j++){
				for (int k = 0; k < 3; k++){
					if (result[3 * j + k] != 0.0){
						int nodeIndexOffset = indexer->getGlobalIndex(*facet, j, k);
						if (nodeIndexOffset >= 0){
							auto found = virtualWorks.find(nodeIndexOffset);
							if (found != virtualWorks.end())
								found->second += result[3 * j + k];
							else
								virtualWorks.insert(std::pair<int, double>(nodeIndexOffset, result[3 * j + k]));
						}
					}
				}
			}
		}
		delete facet;
		delete[] result;
	}

	void Forcer::getVirtualWorks(int dofs, int totalDofs, const double *basises, double *vws) const{
		memset(vws, 0, dofs*sizeof(double));
		if (virtualWorks.empty()) return;
		if (basises){
			const double *basis = basises;
			for (int i = 0; i < dofs; i++){
				for (auto entry : virtualWorks)
					vws[i] += entry.second*basis[entry.first];
				//point to start of next basis
				basis += totalDofs;
			}
		}
		else{
			for (auto entry : virtualWorks)
				vws[entry.first] = entry.second;
		}
	}

	void Forcer::addConstraint(int constraintIndex){
		auto iter = virtualWorks.lower_bound(constraintIndex);
		//check whether every entry is lower
		if (iter == virtualWorks.end()) return;
		//erase it(should never happen), just in case
		if (iter->first == constraintIndex)
			virtualWorks.erase(iter++);
		//change the indices which is higher than constraintIndex
		do{
			int newIndex = iter->first - 1;
			double vw = iter->second;
			virtualWorks.erase(iter++);
			virtualWorks.insert(iter, std::pair<int, double>(newIndex, vw));
		} while (iter != virtualWorks.end());
	}
}

