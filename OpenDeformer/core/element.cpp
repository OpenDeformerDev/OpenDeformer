#include "stdafx.h"
#include "mesh.h"
#include "element.h"

namespace ODER{
	Element::Element(Mesh *m, int index, MarterialType t){
		mesh = m;
		type = t;
		BMatrixs = NULL;
		nodeIndexs = mesh->getElementNodeReference(index);
	}

	Element::Element(Mesh *m, MarterialType t){
		mesh = m;
		type = t;
		BMatrixs = NULL;
	}

	void Element::setNodeIndexs(int elementIndex){
		nodeIndexs = mesh->getElementNodeReference(elementIndex);
	}

	int Element::getNodeIndex(int localIndex) const{
		return nodeIndexs[localIndex];
	}

	Element::~Element(){
		if (BMatrixs)
			delete[] BMatrixs;
	}

	Facet::Facet(Mesh *m, int index) :mesh(m){
		vertIndexs = mesh->getFacetVertReference(index);
	}

	int Facet::getVertIndex(int localIndex) const{
		return vertIndexs[localIndex];
	}
	void Facet::setVertIndexs(int vertIndex){
		vertIndexs = mesh->getFacetVertReference(vertIndex);
	}
}