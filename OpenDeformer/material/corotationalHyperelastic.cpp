#include "stdafx.h"
#include "corotationalHyperelastic.h"
#include "nodeIndexer.h"
#include "mesh.h"
#include "element.h"

namespace ODER {
	CorotationalMaterial::CorotationalMaterial(Scalar rho, Scalar youngs, Scalar possion, Scalar compressingTreshold) 
		: FullOrderNonlinearMaterial(rho, MarterialType(Marterial_Isotropic | Marterial_Corotational)) {
		threshold = compressingTreshold;
		//Lame prameters
		Scalar lambda = possion * youngs / ((Scalar(1.0) + possion) * (Scalar(1.0) - Scalar(2.0) * possion));
		Scalar mu = youngs / (Scalar(2.0) * (Scalar(1.0) + possion));

		D[0] = lambda + Scalar(2.0) * mu;
		D[1] = lambda;
		D[2] = mu;
	}

	Scalar *CorotationalMaterial::getPrecomputes(const Reference<Mesh> &mesh) const {
		CorotationalHyperelasticElement *element = dynamic_cast<CorotationalHyperelasticElement *>(mesh->getMaterialElement(type));
		const int elementCount = mesh->getElementCount();
		const int subMatEntry = element->getInitSubStiffMatEntryCount();
		const int deformGradientEntry = element->getDeformGradientsPreEntryCount();

		Scalar *precomputes = new Scalar[elementCount * (subMatEntry + deformGradientEntry)];
		Scalar *initSubStiffMat = precomputes;
		Scalar *deformationGradientPrecomputed = precomputes + elementCount * subMatEntry;

		for (int i = 0; i < elementCount; i++) {
			element->setNodeIndices(i);
			element->getPrecomputes(D, initSubStiffMat + i * subMatEntry,
				deformationGradientPrecomputed + i * deformGradientEntry);
		}

		delete element;
		return precomputes;
	}
}