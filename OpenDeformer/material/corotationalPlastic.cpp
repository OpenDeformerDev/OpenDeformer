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
		CorotationalHyperelasticElement *element = dynamic_cast<CorotationalHyperelasticElement *>(mesh->getMaterialElement(type));
		const int elementCount = mesh->getElementCount();
		const int subMatEntry = element->getInitSubStiffMatEntryCount();
		const int deformGradientEntry = element->getDeformGradientsPreEntryCount();
		const int plasticStrainEntry = 9;

		const int byteCount = sizeof(Scalar) * elementCount * (subMatEntry + deformGradientEntry + plasticStrainEntry);
		FullOrderNonlinearMaterialCache cache(byteCount);

		Scalar *precomputes = (Scalar *)cache.getMemoryBlock();
		Scalar *initSubStiffMat = precomputes;
		Scalar *deformationGradientPrecomputed = precomputes + elementCount * subMatEntry;

		for (int i = 0; i < elementCount; i++) {
			element->setNodeIndices(i);
			element->getPrecomputes(D, initSubStiffMat + i * subMatEntry,
				deformationGradientPrecomputed + i * deformGradientEntry);
		}

		Scalar *plasticStrains = deformationGradientPrecomputed + elementCount * deformGradientEntry;
		Initiation(plasticStrains, plasticStrainEntry);

		delete[] element;
		return cache;
	}

	void CorotationalPlasticMaterial::computeElasticStress(const Scalar *orthoMat, const Scalar *factoredDeformationGradient,
		Scalar* plasticStrain, Scalar *elasticStress) {
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
				stress[i * 3 + j] += mu2 * elasticStrain[i * 3 + j];
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