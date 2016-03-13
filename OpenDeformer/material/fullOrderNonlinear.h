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

		int *elementNodeIndices = (int *)alloca(3 * numNodesPerElement * sizeof(int));
		for (int elementIndex = 0; elementIndex < mesh->getElementCount(); elementIndex++) {
			element->setNodeIndexs(elementIndex);
			indexer->getElementNodesGlobalIndices(*element, numNodesPerElement, elementNodeIndices);

			for (int subRow = 0; subRow < numNodesPerElement * 3; subRow++) {
				int globalRow = elementNodeIndices[subRow];
				if (globalRow >= 0) {
					for (int subColumn = 0; subColumn <= subRow; subColumn++) {
						int globalColumn = elementNodeIndices[subColumn];
						if (globalColumn >= 0) {
							int rowIndex = std::max(globalRow, globalColumn);
							int columnIndex = std::min(globalRow, globalColumn);
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

		int *elementNodeIndices = (int *)alloca(3 * numNodesPerElement * sizeof(int));
		for (int elementIndex = 0; elementIndex < mesh->getElementCount(); elementIndex++) {
			element->setNodeIndexs(elementIndex);
			indexer->getElementNodesGlobalIndices(*element, numNodesPerElement, elementNodeIndices);

			for (int subRow = 0; subRow < numNodesPerElement * 3; subRow++) {
				int globalRow = elementNodeIndices[subRow];
				if (globalRow >= 0) {
					for (int subColumn = 0; subColumn <= subRow; subColumn++) {
						int globalColumn = elementNodeIndices[subColumn];
						if (globalColumn >= 0) {
							int rowIndex = std::max(globalRow, globalColumn);
							int columnIndex = std::min(globalRow, globalColumn);
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
