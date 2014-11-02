#include "stdafx.h"
#include "mesh.h"
#include "element.h"

namespace ODER{
	Element::Element(Mesh *m, int index, MarterialType t){
		mesh = m;
		type = t;
		BMatrixs = NULL;
		nodeIndexs = &mesh->elements[mesh->numNodesPerElement*index];
	}

	Element::Element(Mesh *m, MarterialType t){
		mesh = m;
		type = t;
		BMatrixs = NULL;
	}

	void Element::setNodeIndexs(int elementIndex){
		nodeIndexs = &mesh->elements[mesh->numNodesPerElement*elementIndex];
	}

	Element::~Element(){
		if (BMatrixs)
			delete[] BMatrixs;
	}

	Facet::Facet(Mesh *m, int index) :mesh(m){
		vertIndexs = &mesh->surfaces[mesh->numVertPerSur*index];
	}

	void Facet::setVertIndexs(int vertIndex){
		vertIndexs = &mesh->surfaces[mesh->numVertPerSur*vertIndex];
	}
}