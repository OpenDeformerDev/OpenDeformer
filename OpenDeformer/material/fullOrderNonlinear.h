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
	struct FullOrderNonlinearMaterialCache {
		FullOrderNonlinearMaterialCache() : memory(NULL) {}
		FullOrderNonlinearMaterialCache(size_t byteCount) { memory = new uint8_t[byteCount]; }
		~FullOrderNonlinearMaterialCache() { delete[] memory; }
		const void *getMemoryBlock() const { return memory; }
		void *getMemoryBlock() { return memory; }

		FullOrderNonlinearMaterialCache(FullOrderNonlinearMaterialCache&& c) {
			memory = c.memory;
			c.memory = NULL;
		}
		FullOrderNonlinearMaterialCache &operator=(FullOrderNonlinearMaterialCache&& c) { 
			std::swap(memory, c.memory); 
			return *this;
		}

		FullOrderNonlinearMaterialCache(const FullOrderNonlinearMaterialCache&) = delete;
		FullOrderNonlinearMaterialCache &operator=(const FullOrderNonlinearMaterialCache&) = delete;
	private:
		uint8_t *memory;
	};

	template<class SpMatrix> class FullOrderNonlinearMaterial : public MechMaterial, public ReferenceCounted {
	public:
		FullOrderNonlinearMaterial(Scalar rho, MarterialType t) : MechMaterial(rho, MarterialType(t | Marterial_NonLinear)) {}
		virtual FullOrderNonlinearMaterialCache getPrecomputes(const Reference<Mesh> &mesh) const { return FullOrderNonlinearMaterialCache(); }
		virtual void generateMatrixAndVirtualWorks(const Reference<Mesh> &mesh, const Reference<NodeIndexer> &indexer,
			const FullOrderNonlinearMaterialCache& cache, const SparseSymMatrixIndicesPerElementCache& matrixIndices, SpMatrix& matrix, Scalar *vws) const = 0;
		void getMatrixStructure(const Reference<Mesh> &mesh, const Reference<NodeIndexer> &indexer, typename SpMatrix::Assembler& assmbler) const;
		SparseSymMatrixIndicesPerElementCache getMatrixIndicesPerElement(const Reference<Mesh> &mesh, const Reference<NodeIndexer> &indexer, const SpMatrix& mat) const;
	};

	template<class SpMatrix> void FullOrderNonlinearMaterial<SpMatrix>::getMatrixStructure(const Reference<Mesh> &mesh, const Reference<NodeIndexer> &indexer,
		typename SpMatrix::Assembler& assmbler) const {
		GeometricElement *element = mesh->getGeometricElement();
		int numNodesPerElement = mesh->getNodePerElementCount();

		int *elementNodeIndices = (int *)alloca(3 * numNodesPerElement * sizeof(int));
		for (int elementIndex = 0; elementIndex < mesh->getElementCount(); elementIndex++) {
			element->setNodeIndices(elementIndex);
			indexer->getElementNodesGlobalIndices(*element, numNodesPerElement, elementNodeIndices);

			for (int subCol = 0; subCol < numNodesPerElement * 3; subCol++) {
				int globalCol = elementNodeIndices[subCol];
				if (globalCol >= 0) {
					for (int subRow = subCol; subRow < numNodesPerElement * 3; subRow++) {
						int globalRow = elementNodeIndices[subRow];
						if (globalRow >= 0) {
							int rowIndex = std::max(globalRow, globalCol);
							int colIndex = std::min(globalRow, globalCol);
							assmbler.addEntry(rowIndex, colIndex, 0);
						}
					}
				}
			}
		}
		delete element;
	}

	template<class SpMatrix> SparseSymMatrixIndicesPerElementCache FullOrderNonlinearMaterial<SpMatrix>::getMatrixIndicesPerElement(const Reference<Mesh> &mesh, 
		const Reference<NodeIndexer> &indexer, const SpMatrix& mat) const {
		auto matrixIndices = mat.getIndices();

		GeometricElement *element = mesh->getGeometricElement();
		int numNodesPerElement = mesh->getNodePerElementCount();
		int numElement = mesh->getElementCount();
		SparseSymMatrixIndicesPerElementCache matIndices(numElement, numNodesPerElement);

		int *elementNodeIndices = (int *)alloca(3 * numNodesPerElement * sizeof(int));
		for (int elementIndex = 0; elementIndex < numElement; elementIndex++) {
			int *localIndices = matIndices.getElementMatIndices(elementIndex);
			element->setNodeIndices(elementIndex);
			indexer->getElementNodesGlobalIndices(*element, numNodesPerElement, elementNodeIndices);

			for (int subCol = 0; subCol < numNodesPerElement * 3; subCol++) {
				int globalCol = elementNodeIndices[subCol];
				if (globalCol >= 0) {
					for (int subRow = subCol; subRow < numNodesPerElement * 3; subRow++) {
						int globalRow = elementNodeIndices[subRow];
						if (globalRow >= 0) {
							int rowIndex = std::max(globalRow, globalCol);
							int colIndex = std::min(globalRow, globalCol);
							localIndices[element->getLocalMatrixIndex(subRow, subCol)] = matrixIndices[colIndex][rowIndex];
						}
						else
							localIndices[element->getLocalMatrixIndex(subRow, subCol)] = -1;
					}
				}
				else {
					for (int subRow = subCol; subRow < numNodesPerElement * 3; subRow++)
						localIndices[element->getLocalMatrixIndex(subRow, subCol)] = -1;
				}
			}
		}
		delete element;
		return matIndices;
	}
}

#endif
