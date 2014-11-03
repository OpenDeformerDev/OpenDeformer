#include "stdafx.h"
#include "mesh.h"
#include "element.h"

namespace ODER{
	Mesh::Mesh(int numN, int numE, int numS, int numNPE, int numVPS)
		:numNodes(numN), numElements(numE), numSurfaces(numS), numNodesPerElement(numNPE), numVertPerSur(numVPS){
		vertices = new Vector[numNodes];
		elements = new int[numNodesPerElement * numElements];
		surfaces = new int[numVertPerSur * numSurfaces];
	}

	void Mesh::setElement(int elementIndex, int *nodeIndices){
		for (int i = 0; i < numNodesPerElement; i++){
			elements[numNodesPerElement * elementIndex + i] = nodeIndices[i];
		}
	}

	void Mesh::setFacet(int facetIndex, int *facetIndices){
		for (int i = 0; i < numVertPerSur; i++){
			surfaces[numVertPerSur * facetIndex + i] = facetIndices[i];
		}
	}

	int Mesh::getElementNodeIndex(int elementIndex, int nodeIndex){
		return elements[elementIndex*numNodesPerElement + nodeIndex];
	}

	const int* Mesh::getElementNodeReference(int elementIndex) const{
		return &elements[elementIndex*numNodesPerElement];
	}

	const int* Mesh::getFacetVertReference(int facetIndex) const{
		return &surfaces[facetIndex*numVertPerSur];
	}

	Mesh::~Mesh(){
		delete[] surfaces;
		delete[] elements;
		delete[] vertices;
	}
}
