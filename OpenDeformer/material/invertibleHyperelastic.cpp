#include "stdafx.h"
#include "invertibleHyperelastic.h"
#include "mesh.h"
#include "element.h"
#include "nodeIndexer.h"
#include "numerMethod.h"
#include "tetelement.h"
#include <iostream>

namespace ODER{
	FullOrderNonlinearMaterialCache InvertibleHyperelasticMaterial::getPrecomputes(const Reference<Mesh> &mesh) const {
		InvertibleHyperelasticElement *element = dynamic_cast<InvertibleHyperelasticElement *>(mesh->getMaterialElement(type));
		const int elementCount = mesh->getElementCount();
		const int drivativeEntry = element->getDirvateEntryCount();
		const int deformGradientEntry = element->getDeformGradientsPreEntryCount();

		const int byteCount = sizeof(Scalar) * elementCount * (drivativeEntry + deformGradientEntry);
		FullOrderNonlinearMaterialCache cache(byteCount);

		Scalar *precomputes = (Scalar *)cache.getMemoryBlock();
		Scalar *shapeFunctionDrivativesPrecomputed = precomputes;
		Scalar *deformationGradientPrecomputed = precomputes + elementCount * drivativeEntry;

		for (int i = 0; i < elementCount; i++) {
			element->setNodeIndices(i);
			element->getPrecomputes(shapeFunctionDrivativesPrecomputed + i * drivativeEntry,
				deformationGradientPrecomputed + i * deformGradientEntry);
		}

		delete element;
		return cache;
	}

	void InvertibleHyperelasticMaterial::generateMatrixAndVirtualWorks(const Reference<Mesh> &mesh, const Reference<NodeIndexer> &indexer,
		const FullOrderNonlinearMaterialCache &cahce, const SparseSymMatrixIndicesPerElementCache &matrixIndices, BlockedSymSpMatrix& matrix, Scalar *vws) const {
		InvertibleHyperelasticElement *element = dynamic_cast<InvertibleHyperelasticElement *>(mesh->getMaterialElement(type));
		const int elementCount = mesh->getElementCount();
		const int nodePerElementCount = mesh->getNodePerElementCount();
		const int drivativeEntry = element->getDirvateEntryCount();
		const int deformGradientEntry = element->getDeformGradientsPreEntryCount();
		const int quadratureCount = element->getQuadraturePointCount();
		const int subMatrixEntryCount = ((3 * nodePerElementCount + 1) * 3 * nodePerElementCount) / 2;

		const int workingEntryCount = subMatrixEntryCount + 3 * nodePerElementCount + (9 + 3 + 9 + 9 + 3 + 6 + 9) * quadratureCount;
		Scalar *memory = new Scalar[workingEntryCount];
		Initiation(memory, workingEntryCount);

		Scalar *subMatrix = memory;
		Scalar *subVirtualWorks = subMatrix + subMatrixEntryCount;
		Scalar *gradients = subVirtualWorks + 3 * nodePerElementCount;
		Scalar *diags = gradients + 9 * quadratureCount;
		Scalar *leftOrthoMats = diags + 3 * quadratureCount;
		Scalar *rightOrthoMats = leftOrthoMats + 9 * quadratureCount;
		Scalar *energyGradient = rightOrthoMats + 9 * quadratureCount;
		Scalar *energyHassian = energyGradient + 3 * quadratureCount;
		Scalar *stresses = energyHassian + 6 * quadratureCount;

		const Scalar *precomputes = (const Scalar *)cahce.getMemoryBlock();
		const Scalar *shapeFunctionDrivativesPrecomputed = precomputes;
		const Scalar *deformationGradientPrecomputed = precomputes + elementCount * drivativeEntry;

		int *elementNodeIndices = (int *)alloca(3 * nodePerElementCount * sizeof(int));

		for (int elementIndex = 0; elementIndex < elementCount; elementIndex++) {
			element->setNodeIndices(elementIndex);
			indexer->getElementNodesGlobalIndices(*element, nodePerElementCount, elementNodeIndices);
			const Scalar *drivatePre = shapeFunctionDrivativesPrecomputed + elementIndex * drivativeEntry;
			const Scalar *gradientPre = deformationGradientPrecomputed + elementIndex * deformGradientEntry;

			element->generateDeformationGradient(gradientPre, gradients);
			for (int quadrature = 0; quadrature < quadratureCount; quadrature++) {
				modifiedDeformGradient(gradients + quadrature * 9, diags + quadrature * 3,
					leftOrthoMats + quadrature * 3, rightOrthoMats + quadrature * 3);

				Scalar invarients[3];
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
			const int *localIndices = matrixIndices.getElementMatIndices(elementIndex);
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
		delete[] memory;
	}

	void InvertibleHyperelasticMaterial::modifiedDeformGradient(const Scalar *gradient, Scalar *diags, Scalar *leftOrthoMat, Scalar *rightOrthoMat) const {
		Scalar triMat[6], eigenVals[3];
		Tensor2<Scalar> UT, VT;
		constexpr Scalar epsilon = Scalar(1e-8);
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
			eigenVals[i] = eigenVals[i] > Scalar(0) ? sqrt(eigenVals[i]) : Scalar(0);

		Tensor2<Scalar> F(gradient);
		VT = UT * F;
		int condition = (eigenVals[0] < epsilon) + ((eigenVals[1] < epsilon) << 1) + ((eigenVals[2] < epsilon) << 2);
		switch (condition){
		case 0:
		{
			for (int i = 0; i < 3; i++) {
				Scalar inverse = Scalar(1.0) / eigenVals[i];
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
			Scalar inverse1 = Scalar(1.0) / eigenVals[1];
			Scalar inverse2 = Scalar(1.0) / eigenVals[2];
			for (int j = 0; j < 3; j++) {
				VT(1, j) *= inverse1;
				VT(2, j) *= inverse2;
			}
			Vector3 another = Vector3(VT(1, 0), VT(1, 1), VT(1, 2)) %
				Vector3(VT(2, 0), VT(2, 1), VT(2, 2));
			memcpy(&VT(0, 0), &another[0], sizeof(Scalar) * 3);
			break;
		}
		case 2:
		{
			Scalar inverse0 = Scalar(1.0) / eigenVals[0];
			Scalar inverse2 = Scalar(1.0) / eigenVals[2];
			for (int j = 0; j < 3; j++) {
				VT(0, j) *= inverse0;
				VT(2, j) *= inverse2;
			}
			Vector3 another = Vector3(VT(2, 0), VT(2, 1), VT(2, 2)) %
				Vector3(VT(0, 0), VT(0, 1), VT(0, 2));
			memcpy(&VT(1, 0), &another[0], sizeof(Scalar) * 3);
			break;
		}
		case 3:
		{
			Scalar inverse = Scalar(1.0) / eigenVals[2];
			for (int i = 0; i < 3; i++)
				VT(2, i) *= inverse;
			Vector3 v1, v2;
			coordinateSystem(Vector3(VT(2, 0), VT(2, 1), VT(2, 2)), v1, v2);
			memcpy(&VT(0, 0), &v1[0], sizeof(Scalar) * 3);
			memcpy(&VT(1, 0), &v2[0], sizeof(Scalar) * 3);
			break;
		}
		case 4:
		{
			Scalar inverse0 = Scalar(1.0) / eigenVals[0];
			Scalar inverse1 = Scalar(1.0) / eigenVals[1];
			for (int j = 0; j < 3; j++) {
				VT(0, j) *= inverse0;
				VT(1, j) *= inverse1;
			}
			Vector3 another = Vector3(VT(0, 0), VT(0, 1), VT(0, 2)) %
				Vector3(VT(1, 0), VT(1, 1), VT(1, 2));
			memcpy(&VT(2, 0), &another[0], sizeof(Scalar) * 3);
			break;
		}
		case 5:
		{
			Scalar inverse = Scalar(1.0) / eigenVals[1];
			for (int i = 0; i < 3; i++)
				VT(1, i) *= inverse;
			Vector3 v1, v2;
			coordinateSystem(Vector3(VT(1, 0), VT(1, 1), VT(1, 2)), v1, v2);
			memcpy(&VT(2, 0), &v1[0], sizeof(Scalar) * 3);
			memcpy(&VT(0, 0), &v2[0], sizeof(Scalar) * 3);
			break;
		}
		case 6:
		{
			Scalar inverse = Scalar(1.0) / eigenVals[0];
			for (int i = 0; i < 3; i++)
				VT(0, i) *= inverse;
			Vector3 v1, v2;
			coordinateSystem(Vector3(VT(0, 0), VT(0, 1), VT(0, 2)), v1, v2);
			memcpy(&VT(1, 0), &v1[0], sizeof(Scalar) * 3);
			memcpy(&VT(2, 0), &v2[0], sizeof(Scalar) * 3);
			break;
		}
		case 7:
		{
			memset(&VT(0, 0), 0, sizeof(Scalar) * 9);
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
		memcpy(rightOrthoMat, &VT(0, 0), sizeof(Scalar) * 9);
	}

	void InvertibleHyperelasticMaterial::getPiolaKirchhoffStress(const Scalar *diags, const Scalar *leftOrthoMat, const Scalar *rightOrthoMat,
		const Scalar *invariants, const Scalar *eneryGradient, Scalar *stress) const {

		Scalar diagEntries[3];
		for (int i = 0; i < 3; i++)
			diagEntries[i] = Scalar(2.0) * eneryGradient[0] * diags[i] + Scalar(4.0) * eneryGradient[1] * diags[i] * diags[i] * diags[i] +
			    Scalar(2.0) * eneryGradient[2] * invariants[2] / diags[i];

		Scalar rightMat[9];
		memcpy(rightMat, rightOrthoMat, sizeof(Scalar) * 9);
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				rightMat[i * 3 + j] *= diagEntries[i];

		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				stress[i * 3 + j] = leftOrthoMat[i * 3 + 0] * rightMat[0 * 3 + j] + 
				    leftOrthoMat[i * 3 + 1] * rightMat[1 * 3 + j] + leftOrthoMat[i * 3 + 2] * rightMat[2 * 3 + j];
	}
}

