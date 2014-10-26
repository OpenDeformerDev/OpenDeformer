#include "stdafx.h"
#include "mecMaterial.h"
#include "element.h"
#include "mesh.h"
#include "nodeIndexer.h"

namespace ODER{
	void MecMaterial::generateMassMatrix(const Reference<Mesh> &mesh, const NodeIndexer &indexer, SparseMatrixAssembler& m) const{
		Element *element = mesh->getEmptyElement();
		int numNodesPerElement = mesh->numNodesPerElement;

		int entrys = ((1 + numNodesPerElement)*numNodesPerElement) >> 1;
		double *subMass = new double[entrys];
		for (int elementIndex = 0; elementIndex < mesh->numElements; elementIndex++){
			element->setNodeIndexs(elementIndex);
			element->setVolume();
			element->generateSubMassMatrix(subMass);

			int k = 0;
			for (int aNodeIndex = 0; aNodeIndex < numNodesPerElement; aNodeIndex++){
				for (int bNodeIndex = 0; bNodeIndex <= aNodeIndex; bNodeIndex++){
					double entry = density*subMass[k++];
					for (int offset = 0; offset < 3; offset++){
						int i_index = indexer.getGlobalIndex(*element, aNodeIndex, offset);
						int j_index = indexer.getGlobalIndex(*element, bNodeIndex, offset);
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
}