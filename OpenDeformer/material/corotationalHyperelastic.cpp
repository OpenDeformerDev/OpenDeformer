#include "stdafx.h"
#include "corotationalHyperelastic.h"
#include "nodeIndexer.h"
#include "mesh.h"
#include "element.h"

namespace ODER {
	CorotationalHyperelasticMaterial::CorotationalHyperelasticMaterial(Scalar rho, Scalar youngs, Scalar possion, Scalar compressingTreshold)
		: FullOrderNonlinearMaterial(rho, MarterialType(Marterial_Isotropic | Marterial_Corotational)) {
		threshold = compressingTreshold;
		//Lame prameters
		Scalar lambda = possion * youngs / ((Scalar(1.0) + possion) * (Scalar(1.0) - Scalar(2.0) * possion));
		Scalar mu = youngs / (Scalar(2.0) * (Scalar(1.0) + possion));

		D[0] = lambda + Scalar(2.0) * mu;
		D[1] = lambda;
		D[2] = mu;
	}

	FullOrderNonlinearMaterialCache CorotationalHyperelasticMaterial::getPrecomputes(const Reference<Mesh> &mesh) const {
		CorotationalHyperelasticElement *element = dynamic_cast<CorotationalHyperelasticElement *>(mesh->getMaterialElement(type));
		const int elementCount = mesh->getElementCount();
		const int subMatEntry = element->getInitSubStiffMatEntryCount();
		const int deformGradientEntry = element->getDeformGradientsPreEntryCount();

		const int byteCount = sizeof(Scalar) * elementCount * (subMatEntry + deformGradientEntry);
		FullOrderNonlinearMaterialCache cache(byteCount);

		Scalar *precomputes = (Scalar *)cache.getMemoryBlock();
		Scalar *initSubStiffMat = precomputes;
		Scalar *deformationGradientPrecomputed = precomputes + elementCount * subMatEntry;

		for (int i = 0; i < elementCount; i++) {
			element->setNodeIndices(i);
			element->getPrecomputes(D, initSubStiffMat + i * subMatEntry,
				deformationGradientPrecomputed + i * deformGradientEntry);
		}

		delete element;
		return cache;
	}

	void CorotationalHyperelasticMaterial::generateMatrixAndVirtualWorks(const Reference<Mesh> &mesh, const Reference<NodeIndexer> &indexer,
		const FullOrderNonlinearMaterialCache& cache, const SparseSymMatrixIndicesPerElementCache& matrixIndices, BlockedSymSpMatrix& matrix, Scalar *vws) const {
		CorotationalHyperelasticElement *element = dynamic_cast<CorotationalHyperelasticElement *>(mesh->getMaterialElement(type));
		const int elementCount = mesh->getElementCount();
		const int nodePerElementCount = mesh->getNodePerElementCount();
		const int initSubMatEntry = element->getInitSubStiffMatEntryCount();
		const int deformGradientEntry = element->getDeformGradientsPreEntryCount();
		const int quadratureCount = element->getQuadraturePointCount();

		const int subMatrixEntryCount = ((3 * nodePerElementCount + 1) * 3 * nodePerElementCount) / 2;

		const int workingEntryCount = subMatrixEntryCount + 3 * nodePerElementCount + 9 * quadratureCount;
		Scalar *memory = new Scalar[workingEntryCount];
		Initiation(memory, workingEntryCount);

		Scalar *subMat = memory;
		Scalar *subVirtualWorks = memory + subMatrixEntryCount;
		Scalar *orthoMats = subVirtualWorks + 3 * nodePerElementCount;

		const Scalar *precomputes = (const Scalar *)cache.getMemoryBlock();
		const Scalar *initSubStiffMat = precomputes;
		const Scalar *deformationGradientPrecomputed = precomputes + elementCount * initSubMatEntry;

		int *elementNodeIndices = (int *)alloca(3 * nodePerElementCount * sizeof(int));

		for (int elementIndex = 0; elementIndex < elementCount; elementIndex++) {
			element->setNodeIndices(elementIndex);
			indexer->getElementNodesGlobalIndices(*element, nodePerElementCount, elementNodeIndices);
			element->generateProperOrthoMats(deformationGradientPrecomputed + elementIndex * deformGradientEntry, threshold, orthoMats);

			//generate subMat and subVirtualWorks
			element->generateSubStiffnessMatrixNodalVirtualWorks(orthoMats, initSubStiffMat + elementIndex * initSubMatEntry, subMat, subVirtualWorks);

			//assmeble to the matrix
			int entryIndex = 0;
			const int *localIndices = matrixIndices.getElementMatIndices(elementIndex);
			for (int subRow = 0; subRow < nodePerElementCount * 3; subRow++) {
				if (elementNodeIndices[subRow] >= 0) {
					for (int subColumn = subRow; subColumn < nodePerElementCount * 3; subColumn++) {
						if (elementNodeIndices[subColumn] >= 0)
							matrix.addEntry(localIndices[entryIndex], subMat[entryIndex]);

						entryIndex += 1;
					}
				}
				else
					entryIndex += nodePerElementCount * 3 - subRow;
			}

			//assmble to the vws
			for (int localIndex = 0; localIndex < nodePerElementCount * 3; localIndex++) {
				int globalIndex = elementNodeIndices[localIndex];
				if (globalIndex >= 0)
					vws[globalIndex] += subVirtualWorks[localIndex];
			}
		}


		delete element;
		delete[] memory;
	}
}