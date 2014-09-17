#include "stdafx.h"
#include "mesh.h"
#include "element.h"

using namespace ODER;

Mesh::Mesh(int numN, int numE, int numS, int numNPE, int numVPS) 
    :numNodes(numN), numElements(numE), numSurfaces(numS), numNodesPerElement(numNPE), numVertPerSur(numVPS){
	vertices = new Vector[numNodes];
	elements = new int[numNodesPerElement * numElements];
	surfaces = new int[numVertPerSur * numS];
}

void Mesh::updateElementInfo(int elementIndex, Element *element) const{
	element->setNodeIndexs(elementIndex);
	element->setVolume();
	element->setBMatrixs();
}

void Mesh::updateFacetInfo(int surIndex, Facet *facet) const{
	facet->vertIndexs = &surfaces[numVertPerSur*surIndex];
}

Mesh::~Mesh(){
	delete[] surfaces;
	delete[] elements;
	delete[] vertices;
}
