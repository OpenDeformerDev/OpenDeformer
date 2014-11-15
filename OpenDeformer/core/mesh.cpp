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

	Mesh::~Mesh(){
		if (surfaces)
		    delete[] surfaces;
		if (elements)
		    delete[] elements;
		if (vertices)
		    delete[] vertices;
	}
}
