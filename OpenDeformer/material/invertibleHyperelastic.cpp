#include "stdafx.h"
#include "invertibleHyperelastic.h"
#include "mesh.h"
#include "element.h"

namespace ODER{
	InvertibleHyperelasticMaterial::InvertibleHyperelasticMaterial(double rho, double inversionTrashold, const Reference<Mesh> &mesh) :
		MechMaterial(rho, MarterialType(Marterial_Isotropic | Marterial_NonLinear | Marterial_Invertible)), trashold(inversionTrashold) {

		element = dynamic_cast<InvertibleHyperelasticElement *>(mesh->getMaterialElement(type));
		const int drivativeEntry = element->getDirvateEntryCount();
		const int deformGradientEntry = element->getDeformGradientsPreEntryCount();
		const int nodalVirtualWorkEntry = element->getNodalVirtualWorksPreEntryCount();
		const int elementCount = mesh->getElementCount();

		shapeFunctionDrivativesPrecomputed = new double[elementCount * drivativeEntry];
		deformationGradientPrecomputed = new double[elementCount * deformGradientEntry];
		nodalVirtualWorkPrecomputed = new double[elementCount * nodalVirtualWorkEntry];

		for (int i = 0; i < elementCount; i++) {
			element->setNodeIndexs(i);
			element->getPrecomputes(shapeFunctionDrivativesPrecomputed + i * drivativeEntry,
				deformationGradientPrecomputed + i * deformGradientEntry, nodalVirtualWorkPrecomputed + i * nodalVirtualWorkEntry);
		}
	}

	InvertibleHyperelasticMaterial::~InvertibleHyperelasticMaterial() {
		delete element;
		delete[] shapeFunctionDrivativesPrecomputed;
		delete[] deformationGradientPrecomputed;
		delete[] nodalVirtualWorkPrecomputed;
	}
}

