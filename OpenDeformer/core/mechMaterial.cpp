#include "stdafx.h"
#include "mechMaterial.h"
#include "element.h"
#include "mesh.h"
#include "nodeIndexer.h"
#include "sparseMatrix.h"

namespace ODER{
	void MechMaterial::generateMassMatrix(const Reference<Mesh> &mesh, const Reference<NodeIndexer> &indexer, SparseMatrixAssembler& m) const{
		GeometricElement *element = mesh->getGeometricElement();
		int numNodesPerElement = mesh->getNodePerElementCount();

		int entrys = ((1 + numNodesPerElement)*numNodesPerElement) >> 1;
		Scalar *subMass = new Scalar[entrys];
		for (int elementIndex = 0; elementIndex < mesh->getElementCount(); elementIndex++){
			element->setNodeIndices(elementIndex);
			element->generateSubMassMatrix(subMass);

			int k = 0;
			for (int aNodeIndex = 0; aNodeIndex < numNodesPerElement; aNodeIndex++){
				for (int bNodeIndex = 0; bNodeIndex <= aNodeIndex; bNodeIndex++){
					Scalar entry = density * subMass[k++];
					for (int offset = 0; offset < 3; offset++){
						int i_index = indexer->getGlobalIndex(*element, aNodeIndex, offset);
						int j_index = indexer->getGlobalIndex(*element, bNodeIndex, offset);
						if (i_index >= 0 && j_index >= 0){
							int rowIndex = std::max(i_index, j_index);
							int columnIndex = std::min(i_index, j_index);
							m.addEntry(rowIndex, columnIndex, entry);
						}
					}
				}
			}
		}
		delete[] subMass;
		delete element;
	}

	void MechMaterial::generateMassMatrix(const Reference<Mesh> &mesh, const Reference<NodeIndexer> &indexer, BlockedSymSpMatrixAssembler& m) const {
		GeometricElement *element = mesh->getGeometricElement();
		int numNodesPerElement = mesh->getNodePerElementCount();

		int entrys = ((1 + numNodesPerElement)*numNodesPerElement) >> 1;
		Scalar *subMass = new Scalar[entrys];
		for (int elementIndex = 0; elementIndex < mesh->getElementCount(); elementIndex++) {
			element->setNodeIndices(elementIndex);
			element->generateSubMassMatrix(subMass);

			int k = 0;
			for (int aNodeIndex = 0; aNodeIndex < numNodesPerElement; aNodeIndex++) {
				for (int bNodeIndex = 0; bNodeIndex <= aNodeIndex; bNodeIndex++) {
					Scalar entry = density*subMass[k++];
					for (int offset = 0; offset < 3; offset++) {
						int i_index = indexer->getGlobalIndex(*element, aNodeIndex, offset);
						int j_index = indexer->getGlobalIndex(*element, bNodeIndex, offset);
						if (i_index >= 0 && j_index >= 0) {
							int rowIndex = std::max(i_index, j_index);
							int columnIndex = std::min(i_index, j_index);
							m.addEntry(rowIndex, columnIndex, entry);
						}
					}
				}
			}
		}
		delete[] subMass;
		delete element;
	}

	void MechMaterial::generateMassMatrix(const Reference<Mesh> &mesh, const Reference<NodeIndexer> &indexer, 
		const SparseSymMatrixIndicesPerElementCache *matrixIndices, BlockedSymSpMatrix& matrix) const {
		GeometricElement *element = mesh->getGeometricElement();
		int numNodesPerElement = mesh->getNodePerElementCount();

		int entrys = ((1 + numNodesPerElement)*numNodesPerElement) >> 1;
		Scalar *subMass = new Scalar[entrys];
		for (int elementIndex = 0; elementIndex < mesh->getElementCount(); elementIndex++) {
			element->setNodeIndices(elementIndex);
			element->generateSubMassMatrix(subMass);

			const int *localIndices = matrixIndices->getElementMatIndices(elementIndex);
			int k = 0;
			for (int aNodeIndex = 0; aNodeIndex < numNodesPerElement; aNodeIndex++) {
				for (int bNodeIndex = 0; bNodeIndex <= aNodeIndex; bNodeIndex++) {
					Scalar entry = density * subMass[k++];
					for (int offset = 0; offset < 3; offset++) {
						int i_index = indexer->getGlobalIndex(*element, aNodeIndex, offset);
						int j_index = indexer->getGlobalIndex(*element, bNodeIndex, offset);
						if (i_index >= 0 && j_index >= 0) {
							int index = localIndices[element->getLocalMatrixIndex(aNodeIndex, offset, bNodeIndex, offset)];
							matrix.addEntry(index, entry);
						}
					}
				}
			}
		}
		delete[] subMass;
		delete element;
	}
}