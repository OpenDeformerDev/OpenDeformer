#include "stdafx.h"
#include "corotationalPlastic.h"
#include "nodeIndexer.h"
#include "mesh.h"
#include "element.h"

namespace ODER {
	CorotationalPlasticMaterial::CorotationalPlasticMaterial(Scalar rho, Scalar youngs, Scalar possion,
		Scalar yieldCri, Scalar creep, Scalar maxStrain, Scalar compressingTreshold) 
		: FullOrderNonlinearMaterial(rho, MarterialType(Marterial_Isotropic | Marterial_Corotational | Marterial_Plastic)) {
		yieldCrierion = yieldCri;
		creepRatio = creep;
		maxPlasticStrain = maxStrain;
		threshold = compressingTreshold;
		//Lame prameters
		Scalar lambda = possion * youngs / ((Scalar(1.0) + possion) * (Scalar(1.0) - Scalar(2.0) * possion));
		Scalar mu = youngs / (Scalar(2.0) * (Scalar(1.0) + possion));

		D[0] = lambda + Scalar(2.0) * mu;
		D[1] = lambda;
		D[2] = mu;
	}

	FullOrderNonlinearMaterialCache CorotationalPlasticMaterial::getPrecomputes(const Reference<Mesh> &mesh) const {
		CorotationalPlasticElement *element = dynamic_cast<CorotationalPlasticElement *>(mesh->getMaterialElement(type));
		const int elementCount = mesh->getElementCount();
		const int subMatEntry = element->getInitSubStiffMatEntryCount();
		const int deformGradientEntry = element->getDeformGradientsPreEntryCount();
		const int drivativeEntry = element->getDirvateEntryCount();
		const int plasticStrainEntry = element->getQuadraturePointCount() * 9;

		const int byteCount = sizeof(Scalar) * elementCount * (subMatEntry + deformGradientEntry + drivativeEntry + plasticStrainEntry);
		FullOrderNonlinearMaterialCache cache(byteCount);

		Scalar *precomputes = (Scalar *)cache.getMemoryBlock();
		Scalar *initSubStiffMat = precomputes;
		Scalar *deformationGradientPrecomputed = precomputes + elementCount * subMatEntry;
		Scalar *drivativePrecomputed = deformationGradientPrecomputed + elementCount * deformGradientEntry;

		for (int i = 0; i < elementCount; i++) {
			element->setNodeIndices(i);
			element->getPrecomputes(D, initSubStiffMat + i * subMatEntry,
				deformationGradientPrecomputed + i * deformGradientEntry, drivativePrecomputed + i * drivativeEntry);
		}

		Scalar *plasticStrainsCache = drivativePrecomputed + elementCount * drivativeEntry;
		Initiation(plasticStrainsCache, elementCount * plasticStrainEntry);

		delete element;
		return cache;
	}

	void CorotationalPlasticMaterial::generateMatrixAndVirtualWorks(const Reference<Mesh> &mesh, const Reference<NodeIndexer> &indexer,
		const FullOrderNonlinearMaterialCache& cache, const SparseSymMatrixIndicesPerElementCache& matrixIndices, BlockedSymSpMatrix& matrix, Scalar *vws) const {
		CorotationalPlasticElement *element = dynamic_cast<CorotationalPlasticElement *>(mesh->getMaterialElement(type));
		const int elementCount = mesh->getElementCount();
		const int nodePerElementCount = mesh->getNodePerElementCount();
		const int initSubMatEntry = element->getInitSubStiffMatEntryCount();
		const int deformGradientEntry = element->getDeformGradientsPreEntryCount();
		const int drivativeEntry = element->getDirvateEntryCount();
		const int quadratureCount = element->getQuadraturePointCount();

		const int subMatrixEntryCount = ((3 * nodePerElementCount + 1) * 3 * nodePerElementCount) / 2;

		const int workingEntryCount = subMatrixEntryCount + 3 * nodePerElementCount + (9 + 9 + 9) * quadratureCount;
		Scalar *memory = new Scalar[workingEntryCount];
		Initiation(memory, workingEntryCount);

		Scalar *subMat = memory;
		Scalar *subVirtualWorks = memory + subMatrixEntryCount;
		Scalar *orthoMats = subVirtualWorks + 3 * nodePerElementCount;
		Scalar *factoredParts = orthoMats + 9 * quadratureCount;
		Scalar *elasticStresses = factoredParts + 9 * quadratureCount;

		Scalar *precomputes = (Scalar *)cache.getMemoryBlock();
		const Scalar *initSubStiffMats = precomputes;
		const Scalar *deformationGradientPrecomputed = precomputes + elementCount * initSubMatEntry;
		const Scalar *drivativePrecomputed = deformationGradientPrecomputed + elementCount * deformGradientEntry;
		Scalar *plasticStrainsCache = (Scalar *)drivativePrecomputed + elementCount * drivativeEntry;

		int *elementNodeIndices = (int *)alloca(3 * nodePerElementCount * sizeof(int));

		for (int elementIndex = 0; elementIndex < elementCount; elementIndex++) {
			element->setNodeIndices(elementIndex);
			indexer->getElementNodesGlobalIndices(*element, nodePerElementCount, elementNodeIndices);
			const Scalar *gradientPre = deformationGradientPrecomputed + elementIndex * deformGradientEntry;
			const Scalar *subMatPre = initSubStiffMats + elementIndex * initSubMatEntry;
			const Scalar *drivatePre = drivativePrecomputed + elementIndex * drivativeEntry;
			Scalar *plasticStrains = plasticStrainsCache + elementIndex * (9 * quadratureCount);

			element->generateDecomposedDeformationGradient(gradientPre, threshold, orthoMats, factoredParts);

			for (int quadrature = 0; quadrature < quadratureCount; quadrature++)
				computeElasticStress(orthoMats + 9 * quadrature, factoredParts + 9 * quadrature,
					plasticStrains + 9 * quadrature, elasticStresses + 9 * quadrature);


			//generate submat
			element->generateSubStiffnessMatrix(orthoMats, subMatPre, subMat);

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

			//generate virtual works
			element->generateNodalVirtualWorks(drivatePre, elasticStresses, subVirtualWorks);
			for (int localIndex = 0; localIndex < nodePerElementCount * 3; localIndex++) {
				int globalIndex = elementNodeIndices[localIndex];
				if (globalIndex >= 0)
					vws[globalIndex] += subVirtualWorks[localIndex];
			}
		}

		delete element;
		delete[] memory;
	}

	void CorotationalPlasticMaterial::computeElasticStress(const Scalar *orthoMat, const Scalar *factoredDeformationGradient,
		Scalar* plasticStrain, Scalar *elasticStress) const {
		Scalar lambda = D[1];
		Scalar mu2 = Scalar(2) * D[2];

		Scalar elasticStrain[9];
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				elasticStrain[i * 3 + j] = Scalar(0.5) * (factoredDeformationGradient[i * 3 + j] + factoredDeformationGradient[j * 3 + i]) - plasticStrain[i * 3 + j];
		elasticStrain[0] -= Scalar(1);
		elasticStrain[4] -= Scalar(1);
		elasticStrain[8] -= Scalar(1);

		Scalar stress[9];
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				stress[i * 3 + j] = mu2 * elasticStrain[i * 3 + j];
		Scalar diag = lambda * (elasticStrain[0] + elasticStrain[4] + elasticStrain[8]);
		stress[0] += diag;
		stress[4] += diag;
		stress[8] += diag;

		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				elasticStress[i * 3 + j] = orthoMat[i * 3 + 0] * stress[0 * 3 + j] + orthoMat[i * 3 + 1] * stress[1 * 3 + j] 
				+ orthoMat[i * 3 + 2] * stress[2 * 3 + j];

		Scalar elasticNorm = 0;
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				elasticNorm += elasticStrain[i] * elasticStrain[i];
		elasticNorm = sqrt(elasticNorm);

		//yield
		if (elasticNorm  > yieldCrierion) {
			Scalar ratio = creepRatio * (elasticNorm - yieldCrierion) / elasticNorm;
			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 3; j++)
					plasticStrain[i * 3 + j] += ratio * elasticStrain[i * 3 + j];

			Scalar plasticNorm = 0;
			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 3; j++)
					plasticNorm += plasticStrain[i] * plasticStrain[i];
			plasticNorm = sqrt(plasticNorm);

			//maximum plastic strain
			if (plasticNorm > maxPlasticStrain) {
				Scalar factor = maxPlasticStrain / plasticNorm;
				for (int i = 0; i < 3; i++)
					for (int j = 0; j < 3; j++)
						plasticStrain[i * 3 + j] *= factor;
			}
		}
	}
}