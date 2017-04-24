#include "stdafx.h"
#include "invertibleHyperelastic.h"
#include "mesh.h"
#include "element.h"
#include "nodeIndexer.h"
#include "numerMethod.h"
#include "tetelement.h"
#include <iostream>

namespace ODER{
	double *InvertibleHyperelasticMaterial::getPrecomputes(const Reference<Mesh> &mesh) const {
		InvertibleHyperelasticElement *element = dynamic_cast<InvertibleHyperelasticElement *>(mesh->getMaterialElement(type));
		const int elementCount = mesh->getElementCount();
		const int drivativeEntry = element->getDirvateEntryCount();
		const int deformGradientEntry = element->getDeformGradientsPreEntryCount();

		double *precomputes = new double[elementCount * (drivativeEntry + deformGradientEntry)];
		double *shapeFunctionDrivativesPrecomputed = precomputes;
		double *deformationGradientPrecomputed = precomputes + elementCount * drivativeEntry;

		for (int i = 0; i < elementCount; i++) {
			element->setNodeIndexs(i);
			element->getPrecomputes(shapeFunctionDrivativesPrecomputed + i * drivativeEntry,
				deformationGradientPrecomputed + i * deformGradientEntry);
		}

		delete element;
		return precomputes;
	}

	void InvertibleHyperelasticMaterial::generateMatrixAndVirtualWorks(const Reference<Mesh> &mesh, const Reference<NodeIndexer> &indexer,
		const double *precomputes, const int *matrixIndices, BlockedSymSpMatrix& matrix, double *vws) const {
		InvertibleHyperelasticElement *element = dynamic_cast<InvertibleHyperelasticElement *>(mesh->getMaterialElement(type));
		const int elementCount = mesh->getElementCount();
		const int nodePerElementCount = mesh->getNodePerElementCount();
		const int drivativeEntry = element->getDirvateEntryCount();
		const int deformGradientEntry = element->getDeformGradientsPreEntryCount();
		const int quadratureCount = element->getQuadraturePointCount();
		const int subMatrixEntryCount = ((3 * nodePerElementCount + 1) * 3 * nodePerElementCount) / 2;

		const int workingEntryCount = subMatrixEntryCount + 3 * nodePerElementCount + (9 + 3 + 9 + 9 + 3 + 6 + 9) * quadratureCount;
		double *memory = new double[workingEntryCount];
		Initiation(memory, workingEntryCount);

		double *subMatrix = memory;
		double *subVirtualWorks = subMatrix + subMatrixEntryCount;
		double *gradients = subVirtualWorks + 3 * nodePerElementCount;
		double *diags = gradients + 9 * quadratureCount;
		double *leftOrthoMats = diags + 3 * quadratureCount;
		double *rightOrthoMats = leftOrthoMats + 9 * quadratureCount;
		double *energyGradient = rightOrthoMats + 9 * quadratureCount;
		double *energyHassian = energyGradient + 3 * quadratureCount;
		double *stresses = energyHassian + 6 * quadratureCount;

		const double *shapeFunctionDrivativesPrecomputed = precomputes;
		const double *deformationGradientPrecomputed = precomputes + elementCount * drivativeEntry;

		int *elementNodeIndices = (int *)alloca(3 * nodePerElementCount * sizeof(int));

		for (int elementIndex = 0; elementIndex < elementCount; elementIndex++) {
			element->setNodeIndexs(elementIndex);
			indexer->getElementNodesGlobalIndices(*element, nodePerElementCount, elementNodeIndices);
			const double *drivatePre = shapeFunctionDrivativesPrecomputed + elementIndex * drivativeEntry;
			const double *gradientPre = deformationGradientPrecomputed + elementIndex * deformGradientEntry;

			element->generateDeformationGradient(gradientPre, gradients);
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
			const int *localIndices = matrixIndices + subMatrixEntryCount * elementIndex;
			for (int subRow = 0; subRow < nodePerElementCount * 3; subRow++) {
				if (elementNodeIndices[subRow] >= 0) {
					for (int subColumn = subRow; subColumn < nodePerElementCount * 3; subColumn++) {
						if (elementNodeIndices[subColumn] >= 0) 
							matrix.addEntry(localIndices[entryIndex], subMatrix[entryIndex]);

						entryIndex += 1;
					}
				}
				else
					entryIndex += nodePerElementCount * 3 - subRow;
			}

			//generate virtual works
			element->generateNodalVirtualWorks(drivatePre, stresses, subVirtualWorks);
			for (int localIndex = 0; localIndex < nodePerElementCount * 3; localIndex++) {
				int globalIndex = elementNodeIndices[localIndex];
				if (globalIndex >= 0)
					vws[globalIndex] += subVirtualWorks[localIndex];
			}

		}

		delete element;
		delete memory;
	}

	void InvertibleHyperelasticMaterial::modifiedDeformGradient(const double *gradient, double *diags, double *leftOrthoMat, double *rightOrthoMat) const {
		double triMat[6], eigenVals[3];
		Tensor2<double> UT, VT;
		constexpr double epsilon = 1e-8;
		constexpr int diagIndices[3] = { 0, 3, 5 };

		for (int i = 0; i < 3; i++)
			for (int j = i; j < 3; j++)
				triMat[diagIndices[i] + j - i] = gradient[j * 3 + 0] * gradient[i * 3 + 0] +
				gradient[j * 3 + 1] * gradient[i * 3 + 1] + gradient[j * 3 + 2] * gradient[i * 3 + 2];

		eigenSym3x3(triMat, eigenVals, &UT(0, 0));
		
		if (UT.Determinant() < 0.0) {
			for (int i = 0; i < 3; i++)
				UT(0, i) = -UT(0, i);
		}

		for (int i = 0; i < 3; i++)
			eigenVals[i] = eigenVals[i] > 0.0 ? sqrt(eigenVals[i]) : 0.0;

		Tensor2<double> F(gradient);
		VT = UT * F;
		int condition = (eigenVals[0] < epsilon) + ((eigenVals[1] < epsilon) << 1) + ((eigenVals[2] < epsilon) << 2);
		switch (condition){
		case 0:
		{
			for (int i = 0; i < 3; i++) {
				double inverse = 1.0 / eigenVals[i];
				for (int j = 0; j < 3; j++)
					VT(i, j) *= inverse;
			}
			if (VT.Determinant() < 0) {
				int smallestIndex = eigenVals[0] < eigenVals[1] ? (eigenVals[0] < eigenVals[2] ? 0 : 2) : (eigenVals[1] < eigenVals[2] ? 1 : 2);
				for (int i = 0; i < 3; i++)
					VT(smallestIndex, i) = -VT(smallestIndex, i);
				eigenVals[smallestIndex] = -eigenVals[smallestIndex];
			}
			break;
		}
		case 1:
		{
			double inverse1 = 1.0 / eigenVals[1];
			double inverse2 = 1.0 / eigenVals[2];
			for (int j = 0; j < 3; j++) {
				VT(1, j) *= inverse1;
				VT(2, j) *= inverse2;
			}
			Vector3d another = Vector3d(VT(1, 0), VT(1, 1), VT(1, 2)) %
				Vector3d(VT(2, 0), VT(2, 1), VT(2, 2));
			memcpy(&VT(0, 0), &another[0], sizeof(double) * 3);
			break;
		}
		case 2:
		{
			double inverse0 = 1.0 / eigenVals[0];
			double inverse2 = 1.0 / eigenVals[2];
			for (int j = 0; j < 3; j++) {
				VT(0, j) *= inverse0;
				VT(2, j) *= inverse2;
			}
			Vector3d another = Vector3d(VT(2, 0), VT(2, 1), VT(2, 2)) %
				Vector3d(VT(0, 0), VT(0, 1), VT(0, 2));
			memcpy(&VT(1, 0), &another[0], sizeof(double) * 3);
			break;
		}
		case 3:
		{
			double inverse = 1.0 / eigenVals[2];
			for (int i = 0; i < 3; i++)
				VT(2, i) *= inverse;
			Vector3d v1, v2;
			coordinateSystem(Vector3d(VT(2, 0), VT(2, 1), VT(2, 2)), v1, v2);
			memcpy(&VT(0, 0), &v1[0], sizeof(double) * 3);
			memcpy(&VT(1, 0), &v2[0], sizeof(double) * 3);
			break;
		}
		case 4:
		{
			double inverse0 = 1.0 / eigenVals[0];
			double inverse1 = 1.0 / eigenVals[1];
			for (int j = 0; j < 3; j++) {
				VT(0, j) *= inverse0;
				VT(1, j) *= inverse1;
			}
			Vector3d another = Vector3d(VT(0, 0), VT(0, 1), VT(0, 2)) %
				Vector3d(VT(1, 0), VT(1, 1), VT(1, 2));
			memcpy(&VT(2, 0), &another[0], sizeof(double) * 3);
			break;
		}
		case 5:
		{
			double inverse = 1.0 / eigenVals[1];
			for (int i = 0; i < 3; i++)
				VT(1, i) *= inverse;
			Vector3d v1, v2;
			coordinateSystem(Vector3d(VT(1, 0), VT(1, 1), VT(1, 2)), v1, v2);
			memcpy(&VT(2, 0), &v1[0], sizeof(double) * 3);
			memcpy(&VT(0, 0), &v2[0], sizeof(double) * 3);
			break;
		}
		case 6:
		{
			double inverse = 1.0 / eigenVals[0];
			for (int i = 0; i < 3; i++)
				VT(0, i) *= inverse;
			Vector3d v1, v2;
			coordinateSystem(Vector3d(VT(0, 0), VT(0, 1), VT(0, 2)), v1, v2);
			memcpy(&VT(1, 0), &v1[0], sizeof(double) * 3);
			memcpy(&VT(2, 0), &v2[0], sizeof(double) * 3);
			break;
		}
		case 7:
		{
			memset(&VT(0, 0), 0, sizeof(double) * 9);
			VT(0, 0) = VT(1, 1) = VT(2, 2) = 1.0;
			break;
		}
		default:
			Severe("Unexpected condition in InvertibleHyperelasticMaterial::modifiedDeformGradient");
			break;
		}

		for (int i = 0; i < 3; i++)
			diags[i] = std::max(eigenVals[i], trashold);
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				leftOrthoMat[i * 3 + j] = UT(j, i);
		memcpy(rightOrthoMat, &VT(0, 0), sizeof(double) * 9);
	}

	void InvertibleHyperelasticMaterial::getPiolaKirchhoffStress(const double *diags, const double *leftOrthoMat, const double *rightOrthoMat,
		const double *invariants, const double *eneryGradient, double *stress) const {

		double diagEntries[3];
		for (int i = 0; i < 3; i++)
			diagEntries[i] = 2.0 * eneryGradient[0] * diags[i] + 4.0 * eneryGradient[1] * diags[i] * diags[i] * diags[i] +
			    2.0 * eneryGradient[2] * invariants[2] / diags[i];

		double rightMat[9];
		memcpy(rightMat, rightOrthoMat, sizeof(double) * 9);
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				rightMat[i * 3 + j] *= diagEntries[i];

		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				stress[i * 3 + j] = leftOrthoMat[i * 3 + 0] * rightMat[0 * 3 + j] + 
				    leftOrthoMat[i * 3 + 1] * rightMat[1 * 3 + j] + leftOrthoMat[i * 3 + 2] * rightMat[2 * 3 + j];
	}
}

