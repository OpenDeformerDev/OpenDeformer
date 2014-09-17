#include "stdafx.h"
#include "tetmesh.h"
#include "tetelement.h"
#include "nodeIndexer.h"

using namespace ODER;

void TetMesh::generateMassMatrix(const MecMaterial &material, const NodeIndexer &indexer, SparseMatrixAssembler &m) const{
	double density = material.getDensity();
	double inv = density / 20.0;

	double mab[10] = { 2.0, 
		               1.0, 2.0, 
		               1.0, 1.0, 2.0, 
		               1.0, 1.0, 1.0, 2.0 };

	double subMass[10];
	TetElement element((TetMesh*)this);

	for (int elementIndex = 0; elementIndex < numElements; elementIndex++){
		element.setNodeIndexs(elementIndex);
		element.setVolume();

		for (int i = 0; i < 10; i++){
			subMass[i] = mab[i] * inv * element.volume;
		}

		int k = 0;
		for (int aNodeIndex = 0; aNodeIndex < 4; aNodeIndex++){
			for (int bNodeIndex = 0; bNodeIndex <= aNodeIndex; bNodeIndex++){
				double entry = subMass[k++];
				for (int offset = 0; offset < 3; offset++){
					int i_index = indexer.getGlobalIndex(element, aNodeIndex, offset);
					int j_index = indexer.getGlobalIndex(element, bNodeIndex, offset);
					if (i_index >= 0 && j_index >= 0){
						int rowIndex = max(i_index, j_index);
						int columnIndex = min(i_index, j_index);
						m.addEntry(rowIndex, columnIndex, entry);
					}
				}
			}
		}
	}
}

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
