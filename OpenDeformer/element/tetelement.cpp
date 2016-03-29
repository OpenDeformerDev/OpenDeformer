#include "stdafx.h"
#include "tetelement.h"
#include "numerMethod.h"

namespace ODER{
	void TetElement::generateSubMassMatrix(double *result) const{
		double volume = getTetVolume(mesh->getVertex(nodeIndexs[0]), mesh->getVertex(nodeIndexs[1]),
			mesh->getVertex(nodeIndexs[2]), mesh->getVertex(nodeIndexs[3]));
		double offdiag = volume*0.05;
		double diag = volume*0.1;
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

	void TetElement::getBodyVirtualWorks(double bodyForce[3], double *result) const{
		double volume = getTetVolume(mesh->getVertex(nodeIndexs[0]), mesh->getVertex(nodeIndexs[1]),
			mesh->getVertex(nodeIndexs[2]), mesh->getVertex(nodeIndexs[3]));
		double factor = volume / 4.0;
		constexpr int numNodesPerElement = 4;
		for (int i = 0; i < numNodesPerElement; i++){
			for (int j = 0; j < 3; j++)
				result[3 * i + j] = factor*bodyForce[j];
		}
	}

	void TetFacet::getSurfVirtualWorks(double surfForce[3], double *result) const{
		double factor = getTriArea(mesh->getVertex(vertIndexs[0]), mesh->getVertex(vertIndexs[1]), mesh->getVertex(vertIndexs[2])) / 3.0;
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

	void LinearIsotropicTetElement::generateSubStiffnessMatrix(int aNodeIndex, int bNodeIndex, const double *D, double *result) const{
		constexpr int entryCountPerB = 3;
		double volume = getTetVolume(mesh->getVertex(nodeIndexs[0]), mesh->getVertex(nodeIndexs[1]),
			mesh->getVertex(nodeIndexs[2]), mesh->getVertex(nodeIndexs[3]));
		double factor = 1.0 / (36.0*volume);
		const double *Ba = BMatrixs + entryCountPerB*aNodeIndex;
		const double *Bb = BMatrixs + entryCountPerB*bNodeIndex;

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

	void ReducedIsotropicTetElement::Intergration(const double *C, double *nlpart, double *nnpart) const{
		constexpr int entryCountPerB = 3;
		constexpr int numNodesPerElement = 4;
		constexpr int commonEntryNum = numNodesPerElement*numNodesPerElement*numNodesPerElement;
		memset(nlpart, 0, commonEntryNum * 3 * sizeof(double));
		memset(nnpart, 0, commonEntryNum * numNodesPerElement * sizeof(double));

		const double determinant = 6.0 * getTetVolume(mesh->getVertex(nodeIndexs[0]), mesh->getVertex(nodeIndexs[1]),
			mesh->getVertex(nodeIndexs[2]), mesh->getVertex(nodeIndexs[3]));
		const double determinant2 = determinant * determinant;
		const double nlFactor = 1.0 / (6.0 * determinant2);
		const double nnFactor = 1.0 / (12.0 * determinant * determinant2);

		const double lambda = C[0], mu2 = 2.0*C[1];

		for (int a = 0; a < numNodesPerElement; a++){
			const double *dNas = BMatrixs + entryCountPerB*a;
			VectorBase<double> dNa(dNas[0], dNas[1], dNas[2]);

			for (int b = 0; b < numNodesPerElement; b++){
				const double *dNbs = BMatrixs + entryCountPerB*b;

				Tensor2<double> t = dNa ^ VectorBase<double>(dNbs[0], dNbs[1], dNbs[2]);
				//double contraction in isotropic form
				double digCommmon = lambda*(t(0, 0) + t(1, 1) + t(2, 2));
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
					const double *dNcs = BMatrixs + entryCountPerB*c;
					VectorBase<double> dNc(dNcs[0], dNcs[1], dNcs[2]);
					VectorBase<double> result = t * dNc;
					int nlsubIndex = (indexOffset + c) * 3;
					nlpart[nlsubIndex + 0] = result[0] * nlFactor;
					nlpart[nlsubIndex + 1] = result[1] * nlFactor;
					nlpart[nlsubIndex + 2] = result[2] * nlFactor;
					//nnpart
					for (int d = 0; d < numNodesPerElement; d++){
						const double *dNds = BMatrixs + entryCountPerB*d;
						Tensor2<double> left = VectorBase<double>(dNds[0], dNds[1], dNds[2]) ^ dNc;
						int nnsubIndex = (indexOffset + c) * 4;
						nnpart[nnsubIndex + d] = (left & t) * nnFactor;
					}
				}
			}
		}	
	}

	void ReducedIsotropicTetElement::generateSubStiffnessMatrix(int aNodeIndex, int bNodeIndex, const double *D, double *result) const{
		constexpr int entryCountPerB = 3;
		double volume = getTetVolume(mesh->getVertex(nodeIndexs[0]), mesh->getVertex(nodeIndexs[1]),
			mesh->getVertex(nodeIndexs[2]), mesh->getVertex(nodeIndexs[3]));
		double factor = 1.0 / (36.0*volume);
		const double *Ba = BMatrixs + entryCountPerB*aNodeIndex;
		const double *Bb = BMatrixs + entryCountPerB*bNodeIndex;

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

	void LinearAnisortropicTetElement::generateSubStiffnessMatrix(int aNodeIndex, int bNodeIndex, const double *D, double *result) const{
		constexpr int entryCountPerB = 3;
		double volume = getTetVolume(mesh->getVertex(nodeIndexs[0]), mesh->getVertex(nodeIndexs[1]),
			mesh->getVertex(nodeIndexs[2]), mesh->getVertex(nodeIndexs[3]));
		double factor = 1.0 / (36.0*volume);
		const double *Ba = BMatrixs + entryCountPerB*aNodeIndex;
		const double *Bb = BMatrixs + entryCountPerB*bNodeIndex;

		double mr[6][3];
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

	void InvertibleHyperelasticTetElement::getPrecomputes(double *drivates, double *deforamtionGradients) const {
		getTetShapeFunctionDerivatives(mesh->getVertex(nodeIndexs[0]), mesh->getVertex(nodeIndexs[1]),
			mesh->getVertex(nodeIndexs[2]), mesh->getVertex(nodeIndexs[3]), drivates, drivates + 3, drivates + 6, drivates + 9);
		
		Vector ab = mesh->getVertex(nodeIndexs[1]) - mesh->getVertex(nodeIndexs[0]);
		Vector ac = mesh->getVertex(nodeIndexs[2]) - mesh->getVertex(nodeIndexs[0]);
		Vector ad = mesh->getVertex(nodeIndexs[3]) - mesh->getVertex(nodeIndexs[0]);

		Tensor2<double> D(ab.x, ac.x, ad.x,
			ab.y, ac.y, ad.y,
			ab.z, ac.z, ad.z);

		memcpy(deforamtionGradients, &(Inverse(D)(0, 0)), sizeof(double) * 9);
	}

	void InvertibleHyperelasticTetElement::generateDeformationGradient(const double *precompute, const double *u, double *gradients) const{
		Initiation(gradients, 9);

		for (int i = 1; i < 4; i++) {
			Vector ax = mesh->getVertex(nodeIndexs[i]) - mesh->getVertex(nodeIndexs[0]);
			for (int j = 0; j < 3; j++) {
				double d = ax[j] + (u[i * 3 + j] - u[j]);
				for (int k = 0; k < 3; k++)
					gradients[j * 3 + k] += d * precompute[(i - 1) * 3 + k];
			}
		}
	}

	void InvertibleHyperelasticTetElement::generateSubStiffnessMatrix(const double *drivates, const double *diags, const double *leftOrthoMats,
		const double *rightOrthoMats, const double *energyGradients, const double *energyHassians, double *result) const{
		constexpr int diagIndices[12] = { 0, 12, 23, 33, 42, 50, 57, 63, 68, 72, 75, 77 };
		constexpr int nodePerElementCount = 4;
		double dPdF[81];
		Initiation(dPdF, 81);
		getdPdF(diags, energyGradients, energyHassians, dPdF);

		double volume = getTetVolume(mesh->getVertex(nodeIndexs[0]), mesh->getVertex(nodeIndexs[1]),
			mesh->getVertex(nodeIndexs[2]), mesh->getVertex(nodeIndexs[3]));
		double factor = 1.0 / (36.0 * volume);

		const Tensor2<double> VT(rightOrthoMats);
		VectorBase<double> dN[4];
		for (int i = 0; i < nodePerElementCount; i++) {
			for (int j = 0; j < 3; j++)
				dN[i][j] = drivates[i * 3 + j];
			dN[i] = VT * dN[i];
		}
	
		const Tensor2<double> U(leftOrthoMats);
		const Tensor2<double> UT = Transpose(U);
		Tensor2<double> subMat;
		for (int aNodeIndex = 0; aNodeIndex < nodePerElementCount; aNodeIndex++) {
			const VectorBase<double> dNa = dN[aNodeIndex];
			//generate diag subMat
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					double entry = 0.0;
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
				const VectorBase<double> dNb = dN[bNodeIndex];
				for (int i = 0; i < 3; i++) {
					for (int j = 0; j < 3; j++) {
						double entry = 0.0;
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

	void InvertibleHyperelasticTetElement::generateNodalVirtualWorks(const double *derivate, const double *stress, double *result) const{
		constexpr int nodePerElementCount = 4;
		constexpr double factor = 1.0 / 6.0;
		Tensor2<double> piolaSress(stress);
		for (int i = 0; i < nodePerElementCount - 1; i++) {
			VectorBase<double> vw = piolaSress * VectorBase<double>(derivate[i * 3 + 0], derivate[i * 3 + 1], derivate[i * 3 + 2]) * factor;
			memcpy(result + i * 3, &vw[0], sizeof(double) * 3);
		}

		for (int i = 0; i < 3; i++)
			result[9 + i] = -(result[0 + i] + result[3 + i] + result[6 + i]);
	}

	void InvertibleHyperelasticTetElement::getdPdF(const double *diag, const double *energyGradient, 
		const double *energyHassian, double dPdF[81]) const {

		const double diagSquare[3] = { diag[0] * diag[0],  diag[1] * diag[1],  diag[2] * diag[2] };
		const double invariant = diagSquare[0] * diagSquare[1] * diagSquare[2];
		const double diagCross[3] = { diag[0] * diag[1], diag[0] * diag[2], diag[1] * diag[2] };

		double alpha00 = 2.0 * energyGradient[0] + 8.0 * diagSquare[0] * energyGradient[1];
		double alpha01 = 2.0 * energyGradient[0] + 4.0 * (diagSquare[0] + diagSquare[1]) * energyGradient[1];
		double alpha02 = 2.0 * energyGradient[0] + 4.0 * (diagSquare[0] + diagSquare[2]) * energyGradient[1];
		double alpha11 = 2.0 * energyGradient[0] + 8.0 * diagSquare[1] * energyGradient[1];
		double alpha12 = 2.0 * energyGradient[0] + 4.0 * (diagSquare[1] + diagSquare[2]) * energyGradient[1];
		double alpha22 = 2.0 * energyGradient[0] + 8.0 * diagSquare[2] * energyGradient[1];

		double delta00 = invariant * energyGradient[2] / diagSquare[0];
		double delta01 = invariant * energyGradient[2] / diagCross[0];
		double delta02 = invariant * energyGradient[2] / diagCross[1];
		double delta11 = invariant * energyGradient[2] / diagSquare[1];
		double delta12 = invariant * energyGradient[2] / diagCross[2];
		double delta22 = invariant * energyGradient[2] / diagSquare[2];

		double beta00 = 4.0 * diagSquare[0] * energyGradient[1] - 2.0 * delta00;
		double beta01 = 4.0 * diagCross[0] * energyGradient[1] - 2.0 * delta01;
		double beta02 = 4.0 * diagCross[1] * energyGradient[1] - 2.0 * delta02;
		double beta11 = 4.0 * diagSquare[1] * energyGradient[1] - 2.0 * delta11;
		double beta12 = 4.0 * diagCross[2] * energyGradient[1] - 2.0 * delta12;
		double beta22 = 4.0 * diagSquare[2] * energyGradient[1] - 2.0 * delta22;

		const Tensor2<double> hassian(energyHassian[0], energyHassian[1], energyHassian[2],
			                          energyHassian[1], energyHassian[3], energyHassian[4],
			                          energyHassian[2], energyHassian[4], energyHassian[5]);

		const VectorBase<double> vec0(2.0 * diag[0], 4.0 * diag[0] * diagSquare[0], 2.0 * invariant / diag[0]);
		const VectorBase<double> vec1(2.0 * diag[1], 4.0 * diag[1] * diagSquare[1], 2.0 * invariant / diag[1]);
		const VectorBase<double> vec2(2.0 * diag[2], 4.0 * diag[2] * diagSquare[2], 2.0 * invariant / diag[2]);

		double gamma00 = vec0 * (hassian * vec0) + 4.0 * delta00;
		double gamma01 = vec1 * (hassian * vec0) + 4.0 * delta01;
		double gamma02 = vec2 * (hassian * vec0) + 4.0 * delta02;
		double gamma11 = vec1 * (hassian * vec1) + 4.0 * delta11;
		double gamma12 = vec2 * (hassian * vec1) + 4.0 * delta12;
		double gamma22 = vec2 * (hassian * vec2) + 4.0 * delta22;

		//ensure (semi-)definiteness
	    double trimat[6];
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

	void InvertibleHyperelasticTetElement::forceSemidefinite3x3(double mat[6]) const {
		Tensor2<double> QT;
		double eigenVals[3];
		eigenSym3x3(mat, eigenVals, &(QT(0, 0)));
		bool negative = false;
		for (int i = 0; i < 3; i++) {
			if (eigenVals[i] < 0) {
				negative = true;
				eigenVals[i] = 0.0;
			}
		}

		if(negative) {
			const Tensor2<double> Q = Transpose(QT);
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++)
					QT(i, j) *= eigenVals[i];
			}
			const Tensor2<double> A = Q * QT;
			mat[0] = A(0, 0);
			mat[1] = A(0, 1);
			mat[2] = A(0, 2);
			mat[3] = A(1, 1);
			mat[4] = A(1, 2);
			mat[5] = A(2, 2);
		}
	}

	void InvertibleHyperelasticTetElement::forceSemidefinite2x2(double& diag, double& offdiag) const {
		double eigenVal0 = diag + offdiag;
		double eigenVal1 = diag - offdiag;

		int condtion = (eigenVal0 < 0) + ((eigenVal1 < 0) << 1);
		switch (condtion) {
		case 0:
			break;
		case 1:
		{
			eigenVal1 *= 0.5;
			diag = eigenVal1;
			offdiag = -eigenVal1;
			break;
		}
		case 2:
		{
			eigenVal0 *= 0.5;
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

	void getTetShapeFunctionDerivatives(const Vector& a, const Vector& b, const Vector& c, const Vector& d,
		double *dn0, double *dn1, double *dn2, double *dn3) {

		Vector da = a - d;
		Vector db = b - d;
		Vector dc = c - d;
		Vector ca = a - c;
		Vector cb = b - c;

		Vector dNa, dNb, dNc, dNd;

		Vector dTest = da % dc + db % da + dc % db;
		if (da * (db % dc) > 0.0) {
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