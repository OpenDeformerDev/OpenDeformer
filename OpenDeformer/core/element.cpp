#include "stdafx.h"
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

	Element::~Element(){
		if (BMatrixs)
			delete[] BMatrixs;
	}

	Facet::Facet(Mesh *m, int index) :mesh(m){
		vertIndexs = mesh->getFacetVertReference(index);
	}
}