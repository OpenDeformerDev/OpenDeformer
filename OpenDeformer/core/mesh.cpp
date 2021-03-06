#include "stdafx.h"
#include "mesh.h"
#include "element.h"
#include "memory.h"

namespace ODER{
	Mesh::Mesh(int numN, int numE, int numS, int numNPE, int numVPS)
		:numNodes(numN), numElements(numE), numSurfaces(numS), numNodesPerElement(numNPE), numVertPerSur(numVPS){
		vertices = new Vector3[numNodes];
		elements = new int[numNodesPerElement * numElements];
		surfaces = new int[numVertPerSur * numSurfaces];
		displacements = new Vector3[numNodes];
		facetElementMap = NULL;
	}

	Mesh::Mesh(Mesh&& m) : numNodes(m.numNodes), numElements(m.numElements), numSurfaces(m.numSurfaces),
		numNodesPerElement(m.numNodesPerElement), numVertPerSur(m.numVertPerSur),
		vertices(m.vertices), elements(m.elements), surfaces(m.surfaces), 
		displacements(m.displacements), facetElementMap(m.facetElementMap){
		m.numNodes = 0; 
		m.numElements = 0; 
		m.numSurfaces = 0;
		m.vertices = NULL; 
		m.elements = NULL; 
		m.surfaces = NULL;
		m.displacements = NULL;
		m.facetElementMap = NULL;
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
		delete[] surfaces;
		delete[] elements;
		delete[] vertices;
		delete[] facetElementMap;
	}
}
