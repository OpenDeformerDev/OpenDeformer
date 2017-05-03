#include "stdafx.h"
#include "tetelement.h"
#include "numerMethod.h"

namespace ODER{
	void TetElement::generateSubMassMatrix(Scalar *result) const{
		Scalar volume = getTetVolume(mesh->getVertex(nodeIndexs[0]), mesh->getVertex(nodeIndexs[1]),
			mesh->getVertex(nodeIndexs[2]), mesh->getVertex(nodeIndexs[3]));
		Scalar offdiag = volume * Scalar(0.05);
		Scalar diag = volume * Scalar(0.1);
		result[0] = diag;
		result[1] = offdiag;
		result[2] = diag;
		result[3] = offdiag;
		result[4] = offdiag;
		result[5] = diag;
		result[6] = offdiag;
		result[7] = offdiag;
		result[8] = offdiag;
		result[9] = diag;
	}

	void TetElement::getBodyVirtualWorks(Scalar bodyForce[3], Scalar *result) const{
		Scalar volume = getTetVolume(mesh->getVertex(nodeIndexs[0]), mesh->getVertex(nodeIndexs[1]),
			mesh->getVertex(nodeIndexs[2]), mesh->getVertex(nodeIndexs[3]));
		Scalar factor = volume / Scalar(4.0);
		constexpr int numNodesPerElement = 4;
		for (int i = 0; i < numNodesPerElement; i++){
			for (int j = 0; j < 3; j++)
				result[3 * i + j] = factor*bodyForce[j];
		}
	}

	void TetFacet::getSurfVirtualWorks(Scalar surfForce[3], Scalar *result) const{
		Vector3 a = mesh->getVertex(vertIndexs[0]) + mesh->getVertexDisplacementConst(vertIndexs[0]);  
		Vector3 b = mesh->getVertex(vertIndexs[1]) + mesh->getVertexDisplacementConst(vertIndexs[1]);
		Vector3 c = mesh->getVertex(vertIndexs[2]) + mesh->getVertexDisplacementConst(vertIndexs[2]);
		Scalar factor = getTriArea(a, b, c) / Scalar(3.0);
		constexpr int numVertPerSur = 3;
		for (int i = 0; i < numVertPerSur; i++){
			for (int j = 0; j < 3; j++)
				result[3 * i + j] = factor*surfForce[j];
		}
	}

	void LinearIsotropicTetElement::setBMatrixs(){
		constexpr int entryCountPerB = 3;
		getTetShapeFunctionDerivatives(mesh->getVertex(nodeIndexs[0]), mesh->getVertex(nodeIndexs[1]),
			mesh->getVertex(nodeIndexs[2]), mesh->getVertex(nodeIndexs[3]),
			BMatrixs, BMatrixs + entryCountPerB, BMatrixs + 2 * entryCountPerB, BMatrixs + 3 * entryCountPerB);
	}

	void LinearIsotropicTetElement::generateSubStiffnessMatrix(int aNodeIndex, int bNodeIndex, const Scalar *D, Scalar *result) const{
		constexpr int entryCountPerB = 3;
		Scalar volume = getTetVolume(mesh->getVertex(nodeIndexs[0]), mesh->getVertex(nodeIndexs[1]),
			mesh->getVertex(nodeIndexs[2]), mesh->getVertex(nodeIndexs[3]));
		Scalar factor = Scalar(1.0) / (Scalar(36.0) * volume);
		const Scalar *Ba = BMatrixs + entryCountPerB*aNodeIndex;
		const Scalar *Bb = BMatrixs + entryCountPerB*bNodeIndex;

		//row first
		result[0] = (D[0] * Ba[0] * Bb[0] + D[2] * (Ba[1] * Bb[1] + Ba[2] * Bb[2])) * factor;
		result[1] = (D[1] * Ba[0] * Bb[1] + D[2] * Ba[1] * Bb[0]) * factor;
		result[2] = (D[1] * Ba[0] * Bb[2] + D[2] * Ba[2] * Bb[0]) * factor;
		result[3] = (D[1] * Ba[1] * Bb[0] + D[2] * Ba[0] * Bb[1]) * factor;
		result[4] = (D[0] * Ba[1] * Bb[1] + D[2] * (Ba[0] * Bb[0] + Ba[2] * Bb[2])) * factor;
		result[5] = (D[1] * Ba[1] * Bb[2] + D[2] * Ba[2] * Bb[1]) * factor;
		result[6] = (D[1] * Ba[2] * Bb[0] + D[2] * Ba[0] * Bb[2]) * factor;
		result[7] = (D[1] * Ba[2] * Bb[1] + D[2] * Ba[1] * Bb[2]) * factor;
		result[8] = (D[0] * Ba[2] * Bb[2] + D[2] * (Ba[0] * Bb[0] + Ba[1] * Bb[1])) * factor;
	}

	void ReducedIsotropicTetElement::setBMatrixs(){
		constexpr int entryCountPerB = 3;
		getTetShapeFunctionDerivatives(mesh->getVertex(nodeIndexs[0]), mesh->getVertex(nodeIndexs[1]),
			mesh->getVertex(nodeIndexs[2]), mesh->getVertex(nodeIndexs[3]),
			BMatrixs, BMatrixs + entryCountPerB, BMatrixs + 2 * entryCountPerB, BMatrixs + 3 * entryCountPerB);
	}

	void ReducedIsotropicTetElement::Intergration(const Scalar *C, Scalar *nlpart, Scalar *nnpart) const{
		constexpr int entryCountPerB = 3;
		constexpr int numNodesPerElement = 4;
		constexpr int commonEntryNum = numNodesPerElement*numNodesPerElement*numNodesPerElement;
		memset(nlpart, 0, commonEntryNum * 3 * sizeof(Scalar));
		memset(nnpart, 0, commonEntryNum * numNodesPerElement * sizeof(Scalar));

		const Scalar determinant = Scalar(6.0) * getTetVolume(mesh->getVertex(nodeIndexs[0]), mesh->getVertex(nodeIndexs[1]),
			mesh->getVertex(nodeIndexs[2]), mesh->getVertex(nodeIndexs[3]));
		const Scalar determinant2 = determinant * determinant;
		const Scalar nlFactor = Scalar(1.0) / (Scalar(6.0) * determinant2);
		const Scalar nnFactor = Scalar(1.0) / (Scalar(12.0) * determinant * determinant2);

		const Scalar lambda = C[0], mu2 = Scalar(2.0) * C[1];

		for (int a = 0; a < numNodesPerElement; a++){
			const Scalar *dNas = BMatrixs + entryCountPerB*a;
			Vector3 dNa(dNas[0], dNas[1], dNas[2]);

			for (int b = 0; b < numNodesPerElement; b++){
				const Scalar *dNbs = BMatrixs + entryCountPerB*b;

				Tensor2<Scalar> t = dNa ^ Vector3(dNbs[0], dNbs[1], dNbs[2]);
				//Scalar contraction in isotropic form
				Scalar digCommmon = lambda*(t(0, 0) + t(1, 1) + t(2, 2));
				t(0, 0) = digCommmon + mu2*t(0, 0);
				t(0, 1) *= mu2;
				t(0, 2) *= mu2;
				t(1, 0) *= mu2;
				t(1, 1) = digCommmon + mu2*t(1, 1);
				t(1, 2) *= mu2;
				t(2, 0) *= mu2;
				t(2, 1) *= mu2;
				t(2, 2) = digCommmon + mu2*t(2, 2);

				int indexOffset = a * 16 + b * 4;
				for (int c = 0; c < numNodesPerElement; c++){
					//nlpart
					const Scalar *dNcs = BMatrixs + entryCountPerB*c;
					Vector3 dNc(dNcs[0], dNcs[1], dNcs[2]);
					Vector3 result = t * dNc;
					int nlsubIndex = (indexOffset + c) * 3;
					nlpart[nlsubIndex + 0] = result[0] * nlFactor;
					nlpart[nlsubIndex + 1] = result[1] * nlFactor;
					nlpart[nlsubIndex + 2] = result[2] * nlFactor;
					//nnpart
					for (int d = 0; d < numNodesPerElement; d++){
						const Scalar *dNds = BMatrixs + entryCountPerB*d;
						Tensor2<Scalar> left = Vector3(dNds[0], dNds[1], dNds[2]) ^ dNc;
						int nnsubIndex = (indexOffset + c) * 4;
						nnpart[nnsubIndex + d] = (left & t) * nnFactor;
					}
				}
			}
		}	
	}

	void ReducedIsotropicTetElement::generateSubStiffnessMatrix(int aNodeIndex, int bNodeIndex, const Scalar *D, Scalar *result) const{
		constexpr int entryCountPerB = 3;
		Scalar volume = getTetVolume(mesh->getVertex(nodeIndexs[0]), mesh->getVertex(nodeIndexs[1]),
			mesh->getVertex(nodeIndexs[2]), mesh->getVertex(nodeIndexs[3]));
		Scalar factor = Scalar(1.0) / (Scalar(36.0) * volume);
		const Scalar *Ba = BMatrixs + entryCountPerB*aNodeIndex;
		const Scalar *Bb = BMatrixs + entryCountPerB*bNodeIndex;

		//row first
		result[0] = (D[0] * Ba[0] * Bb[0] + D[2] * (Ba[1] * Bb[1] + Ba[2] * Bb[2])) * factor;
		result[1] = (D[1] * Ba[0] * Bb[1] + D[2] * Ba[1] * Bb[0]) * factor;
		result[2] = (D[1] * Ba[0] * Bb[2] + D[2] * Ba[2] * Bb[0]) * factor;
		result[3] = (D[1] * Ba[1] * Bb[0] + D[2] * Ba[0] * Bb[1]) * factor;
		result[4] = (D[0] * Ba[1] * Bb[1] + D[2] * (Ba[0] * Bb[0] + Ba[2] * Bb[2])) * factor;
		result[5] = (D[1] * Ba[1] * Bb[2] + D[2] * Ba[2] * Bb[1]) * factor;
		result[6] = (D[1] * Ba[2] * Bb[0] + D[2] * Ba[0] * Bb[2]) * factor;
		result[7] = (D[1] * Ba[2] * Bb[1] + D[2] * Ba[1] * Bb[2]) * factor;
		result[8] = (D[0] * Ba[2] * Bb[2] + D[2] * (Ba[0] * Bb[0] + Ba[1] * Bb[1])) * factor;
	}

	void LinearAnisortropicTetElement::setBMatrixs(){
		constexpr int entryCountPerB = 3;
		getTetShapeFunctionDerivatives(mesh->getVertex(nodeIndexs[0]), mesh->getVertex(nodeIndexs[1]),
			mesh->getVertex(nodeIndexs[2]), mesh->getVertex(nodeIndexs[3]), 
			BMatrixs, BMatrixs + entryCountPerB, BMatrixs + 2 * entryCountPerB, BMatrixs + 3 * entryCountPerB);
	}

	void LinearAnisortropicTetElement::generateSubStiffnessMatrix(int aNodeIndex, int bNodeIndex, const Scalar *D, Scalar *result) const{
		constexpr int entryCountPerB = 3;
		Scalar volume = getTetVolume(mesh->getVertex(nodeIndexs[0]), mesh->getVertex(nodeIndexs[1]),
			mesh->getVertex(nodeIndexs[2]), mesh->getVertex(nodeIndexs[3]));
		Scalar factor = Scalar(1.0) / (Scalar(36.0) * volume);
		const Scalar *Ba = BMatrixs + entryCountPerB*aNodeIndex;
		const Scalar *Bb = BMatrixs + entryCountPerB*bNodeIndex;

		Scalar mr[6][3];
		mr[0][0] = D[0] * Bb[0] + D[3] * Bb[1] + D[5] * Bb[2];
		mr[0][1] = D[1] * Bb[1] + D[3] * Bb[0] + D[4] * Bb[2];
		mr[0][2] = D[2] * Bb[2] + D[4] * Bb[1] + D[5] * Bb[0];

		mr[1][0] = D[1] * Bb[0] + D[8] * Bb[1] + D[10] * Bb[2];
		mr[1][1] = D[6] * Bb[1] + D[8] * Bb[0] + D[9] * Bb[2];
		mr[1][2] = D[7] * Bb[2] + D[9] * Bb[1] * D[10] * Bb[0];

		mr[2][0] = D[2] * Bb[0] + D[12] * Bb[1] + D[14] * Bb[2];
		mr[2][1] = D[7] * Bb[1] + D[12] * Bb[0] + D[13] * Bb[2];
		mr[2][2] = D[11] * Bb[2] + D[13] * Bb[1] + D[14] * Bb[0];

		mr[3][0] = D[3] * Bb[0] + D[15] * Bb[1] + D[17] * Bb[2];
		mr[3][1] = D[8] * Bb[1] + D[15] * Bb[0] + D[16] * Bb[2];
		mr[3][2] = D[12] * Bb[2] + D[16] * Bb[1] + D[17] * Bb[0];

		mr[4][0] = D[4] * Bb[0] + D[16] * Bb[1] + D[19] * Bb[2];
		mr[4][1] = D[9] * Bb[1] + D[16] * Bb[0] + D[18] * Bb[2];
		mr[4][2] = D[13] * Bb[2] + D[18] * Bb[1] + D[19] * Bb[0];

		mr[5][0] = D[5] * Bb[0] + D[17] * Bb[1] + D[20] * Bb[2];
		mr[5][1] = D[10] * Bb[1] + D[17] * Bb[0] + D[19] * Bb[2];
		mr[5][2] = D[14] * Bb[2] + D[19] * Bb[1] + D[20] * Bb[0];

		result[0] = (Ba[0] * mr[0][0] + Ba[1] * mr[3][0] + Ba[2] * mr[5][0])*factor;
		result[1] = (Ba[0] * mr[0][1] + Ba[1] * mr[3][1] + Ba[2] * mr[5][1])*factor;
		result[2] = (Ba[0] * mr[0][2] + Ba[1] * mr[3][2] + Ba[2] * mr[5][2])*factor;
		result[3] = (Ba[1] * mr[1][0] + Ba[0] * mr[3][0] + Ba[2] * mr[4][0])*factor;
		result[4] = (Ba[1] * mr[1][1] + Ba[0] * mr[3][1] + Ba[2] * mr[4][1])*factor;
		result[5] = (Ba[1] * mr[1][2] + Ba[0] * mr[3][2] + Ba[2] * mr[4][2])*factor;
		result[6] = (Ba[2] * mr[2][0] + Ba[1] * mr[4][0] + Ba[0] * mr[5][0])*factor;
		result[7] = (Ba[2] * mr[2][1] + Ba[1] * mr[4][1] + Ba[0] * mr[5][1])*factor;
		result[8] = (Ba[2] * mr[2][2] + Ba[1] * mr[4][2] + Ba[0] * mr[5][2])*factor;
	}

	void InvertibleHyperelasticTetElement::getPrecomputes(Scalar *drivates, Scalar *deforamtionGradients) const {
		getTetShapeFunctionDerivatives(mesh->getVertex(nodeIndexs[0]), mesh->getVertex(nodeIndexs[1]),
			mesh->getVertex(nodeIndexs[2]), mesh->getVertex(nodeIndexs[3]), drivates, drivates + 3, drivates + 6, drivates + 9);
		
		Vector3 ab = mesh->getVertex(nodeIndexs[1]) - mesh->getVertex(nodeIndexs[0]);
		Vector3 ac = mesh->getVertex(nodeIndexs[2]) - mesh->getVertex(nodeIndexs[0]);
		Vector3 ad = mesh->getVertex(nodeIndexs[3]) - mesh->getVertex(nodeIndexs[0]);

		Tensor2<Scalar> D(ab.x, ac.x, ad.x,
			ab.y, ac.y, ad.y,
			ab.z, ac.z, ad.z);

		memcpy(deforamtionGradients, &(Inverse(D)(0, 0)), sizeof(Scalar) * 9);
	}

	void InvertibleHyperelasticTetElement::generateDeformationGradient(const Scalar *precompute, Scalar *gradients) const{
		Initiation(gradients, 9);

		Vector3 v0 = mesh->getVertex(nodeIndexs[0]) + mesh->getVertexDisplacementConst(nodeIndexs[0]);
		for (int i = 1; i < 4; i++) {
			Vector3 ax = (mesh->getVertex(nodeIndexs[i]) + mesh->getVertexDisplacementConst(nodeIndexs[i])) - v0;
			for (int j = 0; j < 3; j++) {
				Scalar d = ax[j];
				for (int k = 0; k < 3; k++)
					gradients[j * 3 + k] += d * precompute[(i - 1) * 3 + k];
			}
		}
	}

	void InvertibleHyperelasticTetElement::generateSubStiffnessMatrix(const Scalar *drivates, const Scalar *diags, const Scalar *leftOrthoMats,
		const Scalar *rightOrthoMats, const Scalar *energyGradients, const Scalar *energyHassians, Scalar *result) const{
		constexpr int diagIndices[12] = { 0, 12, 23, 33, 42, 50, 57, 63, 68, 72, 75, 77 };
		constexpr int nodePerElementCount = 4;
		Scalar dPdF[81];
		Initiation(dPdF, 81);
		getdPdF(diags, energyGradients, energyHassians, dPdF);

		Scalar volume = getTetVolume(mesh->getVertex(nodeIndexs[0]), mesh->getVertex(nodeIndexs[1]),
			mesh->getVertex(nodeIndexs[2]), mesh->getVertex(nodeIndexs[3]));
		Scalar factor = Scalar(1.0) / (Scalar(36.0) * volume);

		const Tensor2<Scalar> VT(rightOrthoMats);
		Vector3 dN[4];
		for (int i = 0; i < nodePerElementCount; i++) {
			for (int j = 0; j < 3; j++)
				dN[i][j] = drivates[i * 3 + j];
			dN[i] = VT * dN[i];
		}
	
		const Tensor2<Scalar> U(leftOrthoMats);
		const Tensor2<Scalar> UT = Transpose(U);
		Tensor2<Scalar> subMat;
		for (int aNodeIndex = 0; aNodeIndex < nodePerElementCount; aNodeIndex++) {
			const Vector3 dNa = dN[aNodeIndex];
			//generate diag subMat
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					Scalar entry = 0.0;
					for (int k = 0; k < 3; k++)
						for (int l = 0; l < 3; l++)
							entry += dPdF[tensorIndex(i, k, j, l)] * dNa[k] * dNa[l];
					subMat(i, j) = entry;
				}
			}
			subMat = U * subMat * UT;
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3 - i; j++)
					result[diagIndices[aNodeIndex * 3 + i] + j] = subMat(i, i + j) * factor;
			}
			
			//generate offdiag subMat
			for (int bNodeIndex = aNodeIndex + 1; bNodeIndex < nodePerElementCount; bNodeIndex++) {
				const Vector3 dNb = dN[bNodeIndex];
				for (int i = 0; i < 3; i++) {
					for (int j = 0; j < 3; j++) {
						Scalar entry = 0.0;
						for (int k = 0; k < 3; k++)
							for (int l = 0; l < 3; l++)
								entry += dPdF[tensorIndex(i, k, j, l)] * dNa[k] * dNb[l];
						subMat(i, j) = entry;
					}
				}
				subMat = U * subMat * UT;
				for (int i = 0; i < 3; i++) {
					for (int j = 0; j < 3; j++)
						result[diagIndices[aNodeIndex * 3 + i] + ((bNodeIndex - aNodeIndex) * 3 - i) + j] = subMat(i, j) * factor;
				}
			}
		}
	}

	void InvertibleHyperelasticTetElement::generateNodalVirtualWorks(const Scalar *derivate, const Scalar *stress, Scalar *result) const{
		constexpr int nodePerElementCount = 4;
		constexpr Scalar factor = Scalar(1.0) / Scalar(6.0);
		Tensor2<Scalar> piolaSress(stress);
		for (int i = 0; i < nodePerElementCount - 1; i++) {
			Vector3 vw = piolaSress * Vector3(derivate[i * 3 + 0], derivate[i * 3 + 1], derivate[i * 3 + 2]) * factor;
			memcpy(result + i * 3, &vw[0], sizeof(Scalar) * 3);
		}

		for (int i = 0; i < 3; i++)
			result[9 + i] = -(result[0 + i] + result[3 + i] + result[6 + i]);
	}

	void InvertibleHyperelasticTetElement::getdPdF(const Scalar *diag, const Scalar *energyGradient, 
		const Scalar *energyHassian, Scalar dPdF[81]) const {

		const Scalar diagSquare[3] = { diag[0] * diag[0],  diag[1] * diag[1],  diag[2] * diag[2] };
		const Scalar invariant = diagSquare[0] * diagSquare[1] * diagSquare[2];
		const Scalar diagCross[3] = { diag[0] * diag[1], diag[0] * diag[2], diag[1] * diag[2] };

		Scalar alpha00 = Scalar(2.0) * energyGradient[0] + Scalar(8.0) * diagSquare[0] * energyGradient[1];
		Scalar alpha01 = Scalar(2.0) * energyGradient[0] + Scalar(4.0) * (diagSquare[0] + diagSquare[1]) * energyGradient[1];
		Scalar alpha02 = Scalar(2.0) * energyGradient[0] + Scalar(4.0) * (diagSquare[0] + diagSquare[2]) * energyGradient[1];
		Scalar alpha11 = Scalar(2.0) * energyGradient[0] + Scalar(8.0) * diagSquare[1] * energyGradient[1];
		Scalar alpha12 = Scalar(2.0) * energyGradient[0] + Scalar(4.0) * (diagSquare[1] + diagSquare[2]) * energyGradient[1];
		Scalar alpha22 = Scalar(2.0) * energyGradient[0] + Scalar(8.0) * diagSquare[2] * energyGradient[1];

		Scalar delta00 = invariant * energyGradient[2] / diagSquare[0];
		Scalar delta01 = invariant * energyGradient[2] / diagCross[0];
		Scalar delta02 = invariant * energyGradient[2] / diagCross[1];
		Scalar delta11 = invariant * energyGradient[2] / diagSquare[1];
		Scalar delta12 = invariant * energyGradient[2] / diagCross[2];
		Scalar delta22 = invariant * energyGradient[2] / diagSquare[2];

		Scalar beta00 = Scalar(4.0) * diagSquare[0] * energyGradient[1] - Scalar(2.0) * delta00;
		Scalar beta01 = Scalar(4.0) * diagCross[0] * energyGradient[1] - Scalar(2.0) * delta01;
		Scalar beta02 = Scalar(4.0) * diagCross[1] * energyGradient[1] - Scalar(2.0) * delta02;
		Scalar beta11 = Scalar(4.0) * diagSquare[1] * energyGradient[1] - Scalar(2.0) * delta11;
		Scalar beta12 = Scalar(4.0) * diagCross[2] * energyGradient[1] - Scalar(2.0) * delta12;
		Scalar beta22 = Scalar(4.0) * diagSquare[2] * energyGradient[1] - Scalar(2.0) * delta22;

		const Tensor2<Scalar> hassian(energyHassian[0], energyHassian[1], energyHassian[2],
			                          energyHassian[1], energyHassian[3], energyHassian[4],
			                          energyHassian[2], energyHassian[4], energyHassian[5]);

		const Vector3 vec0(Scalar(2.0) * diag[0], Scalar(4.0) * diag[0] * diagSquare[0], Scalar(2.0) * invariant / diag[0]);
		const Vector3 vec1(Scalar(2.0) * diag[1], Scalar(4.0) * diag[1] * diagSquare[1], Scalar(2.0) * invariant / diag[1]);
		const Vector3 vec2(Scalar(2.0) * diag[2], Scalar(4.0) * diag[2] * diagSquare[2], Scalar(2.0) * invariant / diag[2]);

		Scalar gamma00 = vec0 * (hassian * vec0) + Scalar(4.0) * delta00;
		Scalar gamma01 = vec1 * (hassian * vec0) + Scalar(4.0) * delta01;
		Scalar gamma02 = vec2 * (hassian * vec0) + Scalar(4.0) * delta02;
		Scalar gamma11 = vec1 * (hassian * vec1) + Scalar(4.0) * delta11;
		Scalar gamma12 = vec2 * (hassian * vec1) + Scalar(4.0) * delta12;
		Scalar gamma22 = vec2 * (hassian * vec2) + Scalar(4.0) * delta22;

		//ensure (semi-)definiteness
	    Scalar trimat[6];
		trimat[0] = alpha00 + beta00 + gamma00;
		trimat[1] = gamma01;
		trimat[2] = gamma02;
		trimat[3] = alpha11 + beta11 + gamma11;
		trimat[4] = gamma12;
		trimat[5] = alpha22 + beta22 + gamma22;
		forceSemidefinite3x3(trimat);
		forceSemidefinite2x2(alpha01, beta01);
		forceSemidefinite2x2(alpha02, beta02);
		forceSemidefinite2x2(alpha12, beta12);

		dPdF[tensorIndex(0, 0, 0, 0)] = trimat[0];
		dPdF[tensorIndex(0, 0, 1, 1)] = trimat[1];
		dPdF[tensorIndex(0, 0, 2, 2)] = trimat[2];
		dPdF[tensorIndex(1, 1, 0, 0)] = trimat[1];
		dPdF[tensorIndex(1, 1, 1, 1)] = trimat[3];
		dPdF[tensorIndex(1, 1, 2, 2)] = trimat[4];
		dPdF[tensorIndex(2, 2, 0, 0)] = trimat[2];
		dPdF[tensorIndex(2, 2, 1, 1)] = trimat[4];
		dPdF[tensorIndex(2, 2, 2, 2)] = trimat[5];

		dPdF[tensorIndex(0, 1, 0, 1)] = alpha01;
		dPdF[tensorIndex(0, 1, 1, 0)] = beta01;
		dPdF[tensorIndex(1, 0, 0, 1)] = beta01;
		dPdF[tensorIndex(1, 0, 1, 0)] = alpha01;

		dPdF[tensorIndex(0, 2, 0, 2)] = alpha02;
		dPdF[tensorIndex(0, 2, 2, 0)] = beta02;
		dPdF[tensorIndex(2, 0, 0, 2)] = beta02;
		dPdF[tensorIndex(2, 0, 2, 0)] = alpha02;

		dPdF[tensorIndex(1, 2, 1, 2)] = alpha12;
		dPdF[tensorIndex(1, 2, 2, 1)] = beta12;
		dPdF[tensorIndex(2, 1, 1, 2)] = beta12;
		dPdF[tensorIndex(2, 1, 2, 1)] = alpha12;
	}

	void InvertibleHyperelasticTetElement::forceSemidefinite3x3(Scalar mat[6]) const {
		Tensor2<Scalar> QT;
		Scalar eigenVals[3];
		eigenSym3x3(mat, eigenVals, &(QT(0, 0)));
		bool negative = false;
		for (int i = 0; i < 3; i++) {
			if (eigenVals[i] < 0) {
				negative = true;
				eigenVals[i] = 0.0;
			}
		}

		if (negative) {
			const Tensor2<Scalar> Q = Transpose(QT);
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++)
					QT(i, j) *= eigenVals[i];
			}
			const Tensor2<Scalar> A = Q * QT;
			mat[0] = A(0, 0);
			mat[1] = A(0, 1);
			mat[2] = A(0, 2);
			mat[3] = A(1, 1);
			mat[4] = A(1, 2);
			mat[5] = A(2, 2);
		}
	}

	void InvertibleHyperelasticTetElement::forceSemidefinite2x2(Scalar& diag, Scalar& offdiag) const {
		Scalar eigenVal0 = diag + offdiag;
		Scalar eigenVal1 = diag - offdiag;

		int condtion = (eigenVal0 < 0) + ((eigenVal1 < 0) << 1);
		switch (condtion) {
		case 0:
			break;
		case 1:
		{
			eigenVal1 *= Scalar(0.5);
			diag = eigenVal1;
			offdiag = -eigenVal1;
			break;
		}
		case 2:
		{
			eigenVal0 *= Scalar(0.5);
			diag = eigenVal0;
			offdiag = eigenVal0;
			break;
		}
		case 3:
		{
			diag = 0.0;
			offdiag = 0.0;
			break;
		}
		default:
			Severe("Unexpected condition in InvertibleHyperelasticTetElement::forceSemidefinite2x2");
			break;
		}
	}

	void getTetShapeFunctionDerivatives(const Vector3& a, const Vector3& b, const Vector3& c, const Vector3& d,
		Scalar *dn0, Scalar *dn1, Scalar *dn2, Scalar *dn3) {

		Vector3 da = a - d;
		Vector3 db = b - d;
		Vector3 dc = c - d;
		Vector3 ca = a - c;
		Vector3 cb = b - c;

		Vector3 dNa, dNb, dNc, dNd;

		if (da * (db % dc) > Scalar(0)) {
			dNa = db % dc;
			dNb = dc % da;
			dNc = da % db;
			dNd = cb % ca;
		}
		else {
			dNa = dc % db;
			dNb = da % dc;
			dNc = db % da;
			dNd = ca % cb;
		}

		for (int i = 0; i < 3; i++) {
			dn0[i] = dNa[i];
			dn1[i] = dNb[i];
			dn2[i] = dNc[i];
			dn3[i] = dNd[i];
		}
	}
}