#include "stdafx.h"
#include "invertibleHyperelastic.h"
#include "mesh.h"
#include "element.h"
#include "nodeIndexer.h"

namespace ODER{
	InvertibleHyperelasticMaterial::InvertibleHyperelasticMaterial(double rho, double inversionTrashold, const Reference<Mesh> &mesh) :
		MechMaterial(rho, MarterialType(Marterial_Isotropic | Marterial_NonLinear | Marterial_Invertible)), trashold(inversionTrashold) {

		element = dynamic_cast<InvertibleHyperelasticElement *>(mesh->getMaterialElement(type));
		const int drivativeEntry = element->getDirvateEntryCount();
		const int deformGradientEntry = element->getDeformGradientsPreEntryCount();
		const int nodalVirtualWorkEntry = element->getNodalVirtualWorksPreEntryCount();
		const int elementCount = mesh->getElementCount();
		const int nodePerElementCount = mesh->getNodePerElementCount();
		const int quadratureCount = element->getQuadraturePointCount();

		shapeFunctionDrivativesPrecomputed = new double[elementCount * drivativeEntry];
		deformationGradientPrecomputed = new double[elementCount * deformGradientEntry];
		nodalVirtualWorkPrecomputed = new double[elementCount * nodalVirtualWorkEntry];

		const int workingEntryCount = ((3 * nodePerElementCount + 1) * 3 * nodePerElementCount) / 2 +
			6 * nodePerElementCount + (9 + 3 + 9 + 9 + 3 + 6 + 9) * quadratureCount;
		memory = new double[workingEntryCount];

		for (int i = 0; i < elementCount; i++) {
			element->setNodeIndexs(i);
			element->getPrecomputes(shapeFunctionDrivativesPrecomputed + i * drivativeEntry,
				deformationGradientPrecomputed + i * deformGradientEntry, nodalVirtualWorkPrecomputed + i * nodalVirtualWorkEntry);
		}
	}

	void InvertibleHyperelasticMaterial::generateMatrixAndVirtualWorks(const Reference<Mesh> &mesh, const Reference<NodeIndexer> &indexer,
		const double *u, BlockedSymSpMatrix& matrix, double *vws) const {

		const int elementCount = mesh->getElementCount();
		const int nodePerElementCount = mesh->getNodePerElementCount();
		const int drivativeEntry = element->getDirvateEntryCount();
		const int deformGradientEntry = element->getDeformGradientsPreEntryCount();
		const int nodalVirtualWorkEntry = element->getNodalVirtualWorksPreEntryCount();
		const int quadratureCount = element->getQuadraturePointCount();
		const int subMatrixEntryCount = ((3 * nodePerElementCount + 1) * 3 * nodePerElementCount) / 2;

		Initiation(vws, 3 * elementCount);
		double *subMatrix = memory;
		double *subVirtualWorks = subMatrix + subMatrixEntryCount;
		double *nodeDisplacements = subVirtualWorks + 3 * nodePerElementCount;
		double *gradients = nodeDisplacements + 3 * nodePerElementCount;
		double *diags = gradients + 9 * quadratureCount;
		double *leftOrthoMats = diags + 3 * quadratureCount;
		double *rightOrthoMats = leftOrthoMats + 9 * quadratureCount;
		double *energyGradient = rightOrthoMats + 9 * quadratureCount;
		double *energyHassian = energyGradient + 3 * quadratureCount;
		double *stresses = energyHassian + 6 * quadratureCount;

		int *elementNodeIndices = (int *)alloca(3 * nodePerElementCount * sizeof(int));

		for (int elementIndex = 0; elementIndex < mesh->getElementCount(); elementIndex++) {
			element->setNodeIndexs(elementIndex);
			indexer->getElementNodesGlobalIndices(*element, nodePerElementCount, elementNodeIndices);
			const double *drivatePre = shapeFunctionDrivativesPrecomputed + elementIndex * drivativeEntry;
			const double *gradientPre = deformationGradientPrecomputed + elementIndex * deformGradientEntry;
			const double *virtualWorkPre = nodalVirtualWorkPrecomputed + elementIndex * nodalVirtualWorkEntry;

			//copy node displacements
			Initiation(nodeDisplacements, 3 * nodePerElementCount);
			for (int nodeIndex = 0; nodeIndex < nodePerElementCount; nodeIndex++) {
				if(elementNodeIndices[nodeIndex] >= 0)
				  memcpy(nodeDisplacements + 3 * nodeIndex, u + 3 * elementNodeIndices[nodeIndex], sizeof(double) * 3);
			}
			
			element->generateDeformationGradient(gradientPre, nodeDisplacements, gradients);
			for (int quadrature = 0; quadrature < quadratureCount; quadrature++) {
				modifiedDeformGradient(gradients + quadrature * 9, diags + quadrature * 3,
					leftOrthoMats + quadrature * 3, rightOrthoMats + quadrature * 3);

				double invarients[3];
				getInvariants(diags + quadrature * 3, invarients);
				getEnergyGradient(invarients, energyGradient + 3 * quadrature);
				getEnergyHassian(invarients, energyHassian + 6 * quadrature);
				getPiolaKirchhoffStress(diags + quadrature * 3, leftOrthoMats + quadrature * 3, rightOrthoMats + quadrature * 3,
					invarients, energyGradient + 3 * quadrature, stresses + 9 * quadrature);
			}

			//generate stifness matrix
			Initiation(subMatrix, subMatrixEntryCount);
			element->generateSubStiffnessMatrix(drivatePre, diags, leftOrthoMats, rightOrthoMats,
				energyGradient, energyHassian, subMatrix);
			//add matrix entries
			int entryIndex = 0;
			for (int subRow = 0; subRow < nodePerElementCount * 3; subRow++) {
				int globalRow = elementNodeIndices[subRow];
				if (globalRow >= 0) {
					for (int subColumn = 0; subColumn <= subRow; subColumn++) {
						int globalColumn = elementNodeIndices[subColumn];
						double matEntry = subMatrix[entryIndex++];
						if (globalColumn >= 0 && matEntry != 0.0)
							matrix.addEntry(globalRow, globalColumn, matEntry);					
					}
				}
				else
					entryIndex += subRow + 1;
			}

			//generate virtual works
			element->generateNodalVirtualWorks(virtualWorkPre, stresses, subVirtualWorks);
			for (int localIndex = 0; localIndex < nodePerElementCount * 3; localIndex++) {
				int globalIndex = elementNodeIndices[localIndex];
				if (globalIndex >= 0)
					vws[globalIndex] += subVirtualWorks[localIndex];
			}

		}
	}

	void InvertibleHyperelasticMaterial::modifiedDeformGradient(const double *gradient, double *diags, double *leftOrthoMat, double *rightOrthoMat) const {

	}

	void InvertibleHyperelasticMaterial::getPiolaKirchhoffStress(const double *diags, const double *leftOrthoMat, const double *rightOrthoMat,
		const double *invariants, const double *eneryGradient, double *stress) const {
		Initiation(stress, 9);
		double diagEntries[3];
		for (int i = 0; i < 3; i++)
			diagEntries[i] = 2.0 * eneryGradient[0] * diags[i] + 4.0 * eneryGradient[1] * diags[i] * diags[i] * diags[i] +
			    2.0 * eneryGradient[2] * invariants[2] / diags[i];

		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				stress[i * 3 + j] = leftOrthoMat[i * 3 + j] * diagEntries[j] * rightOrthoMat[i * 3 + j];
	}

	InvertibleHyperelasticMaterial::~InvertibleHyperelasticMaterial() {
		delete element;
		delete[] shapeFunctionDrivativesPrecomputed;
		delete[] deformationGradientPrecomputed;
		delete[] nodalVirtualWorkPrecomputed;
		delete[] memory;
	}
}

