#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_MATERIAL_FULLORDERNONLINEAR_H
#define ODER_MATERIAL_FULLORDERNONLINEAR_H

#include "mechMaterial.h"
#include "element.h"
#include "nodeIndexer.h"
#include "sparseMatrix.h"

namespace ODER {
	template<class SpMatrix> class FullOrderNonlinearMaterial : public MechMaterial {
	public:
		FullOrderNonlinearMaterial(double rho, MarterialType t) : MechMaterial(rho, MarterialType(t | Marterial_NonLinear)) {}
		virtual void generateMatrixAndVirtualWorks(const Reference<Mesh> &mesh, const Reference<NodeIndexer> &indexer,
			const double *u, const std::vector<std::unordered_map<int, int>> &matrixIndices,
			SpMatrix& matrix, double *vws) const = 0;
		void getMatrixStructure(const Reference<Mesh> &mesh, const Reference<NodeIndexer> &indexer, SparseMatrixAssembler& assmbler) const;
		void getMatrixStructure(const Reference<Mesh> &mesh, const Reference<NodeIndexer> &indexer, BlockedSymSpMatrixAssembler& assmbler) const;
	};

	template<class SparseMatrix> void FullOrderNonlinearMaterial<SparseMatrix>::getMatrixStructure(const Reference<Mesh> &mesh, const Reference<NodeIndexer> &indexer, SparseMatrixAssembler& assmbler) const {
		GeometricElement *element = mesh->getGeometricElement();
		int numNodesPerElement = mesh->getNodePerElementCount();

		int entrys = ((1 + numNodesPerElement) * numNodesPerElement) >> 1;
		for (int elementIndex = 0; elementIndex < mesh->getElementCount(); elementIndex++) {
			element->setNodeIndexs(elementIndex);

			int k = 0;
			for (int aNodeIndex = 0; aNodeIndex < numNodesPerElement; aNodeIndex++) {
				for (int bNodeIndex = 0; bNodeIndex <= aNodeIndex; bNodeIndex++) {
					for (int offset = 0; offset < 3; offset++) {
						int i_index = indexer->getGlobalIndex(*element, aNodeIndex, offset);
						int j_index = indexer->getGlobalIndex(*element, bNodeIndex, offset);
						if (i_index >= 0 && j_index >= 0) {
							int rowIndex = std::max(i_index, j_index);
							int columnIndex = std::min(i_index, j_index);
							assmbler.addEntry(rowIndex, columnIndex, 0.0);
						}
					}
				}
			}
		}
		delete element;
	}

	template<class SparseMatrix> void FullOrderNonlinearMaterial<SparseMatrix>::getMatrixStructure(const Reference<Mesh> &mesh, const Reference<NodeIndexer> &indexer, BlockedSymSpMatrixAssembler& assmbler) const {
		GeometricElement *element = mesh->getGeometricElement();
		int numNodesPerElement = mesh->getNodePerElementCount();

		int entrys = ((1 + numNodesPerElement) * numNodesPerElement) >> 1;
		for (int elementIndex = 0; elementIndex < mesh->getElementCount(); elementIndex++) {
			element->setNodeIndexs(elementIndex);

			int k = 0;
			for (int aNodeIndex = 0; aNodeIndex < numNodesPerElement; aNodeIndex++) {
				for (int bNodeIndex = 0; bNodeIndex <= aNodeIndex; bNodeIndex++) {
					for (int offset = 0; offset < 3; offset++) {
						int i_index = indexer->getGlobalIndex(*element, aNodeIndex, offset);
						int j_index = indexer->getGlobalIndex(*element, bNodeIndex, offset);
						if (i_index >= 0 && j_index >= 0) {
							int rowIndex = std::max(i_index, j_index);
							int columnIndex = std::min(i_index, j_index);
							assmbler.addEntry(rowIndex, columnIndex, 0.0);
						}
					}
				}
			}
		}
		delete element;
	}
}

#endif
