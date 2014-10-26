#include "stdafx.h"
#include "tetmesh.h"
#include "tetelement.h"

namespace ODER{
	Element* TetMesh::getEmptyElement() const{
		return new TetElement((TetMesh *)this);
	}

	Element* TetMesh::getEmptyMaterialElement(MarterialType type) const{
		return new TetElement((TetMesh *)this, type);
	}


	Facet* TetMesh::getEmptyFacet() const{
		return new TetFacet((TetMesh *)this);
	}

	Element* TetMesh::getElement(int elementIndex, MarterialType type) const{
		TetElement* ret = new TetElement((TetMesh *)this, elementIndex, type);
		ret->setVolume();
		return ret;
	}
}
