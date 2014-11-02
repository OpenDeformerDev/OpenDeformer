#include "stdafx.h"
#include "mesh.h"
#include "element.h"

namespace ODER{
	Mesh::Mesh(int numN, int numE, int numS, int numNPE, int numVPS)
		:numNodes(numN), numElements(numE), numSurfaces(numS), numNodesPerElement(numNPE), numVertPerSur(numVPS){
		vertices = new Vector[numNodes];
		elements = new int[numNodesPerElement * numElements];
		surfaces = new int[numVertPerSur * numS];
	}

	Mesh::~Mesh(){
		delete[] surfaces;
		delete[] elements;
		delete[] vertices;
	}
}
