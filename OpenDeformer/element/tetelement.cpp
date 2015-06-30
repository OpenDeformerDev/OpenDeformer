#include "stdafx.h"
#include "tetelement.h"

namespace ODER{
	void TetElement::getShapeFunctionDerivatives(double *dndx, double *dndy, double *dndz) const{
		Vector a = mesh->getVertex(nodeIndexs[0]);
		Vector b = mesh->getVertex(nodeIndexs[1]);
		Vector c = mesh->getVertex(nodeIndexs[2]);
		Vector d = mesh->getVertex(nodeIndexs[3]);

		Vector da = a - d;
		Vector db = b - d;
		Vector dc = c - d;
		Vector ca = a - c;
		Vector cb = b - c;

		dndx[0] = -(db.y*dc.z - db.z*dc.y);
		dndx[1] = (da.y*dc.z - da.z*dc.y);
		dndx[2] = -(da.y*db.z - da.z*db.y);
		dndx[3] = (ca.y*cb.z - ca.z*cb.y);

		dndy[0] = (db.x*dc.z - db.z*dc.x);
		dndy[1] = -(da.x*dc.z - da.z*dc.x);
		dndy[2] = (da.x*db.z - da.z*db.x);
		dndy[3] = -(ca.x*cb.z - ca.z*cb.x);

		dndz[0] = -(db.x*dc.y - db.y*dc.x);
		dndz[1] = (da.x*dc.y - da.y*dc.x);
		dndz[2] = -(da.x*db.y - da.y*db.x);
		dndz[3] = (ca.x*cb.y - ca.y*cb.x);
	}


	void TetElement::generateSubMassMatrix(double *result) const{
		double volume = getVolume();
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


	float TetElement::getVolume() const{
		Vector a = mesh->getVertex(nodeIndexs[0]);
		Vector b = mesh->getVertex(nodeIndexs[1]);
		Vector c = mesh->getVertex(nodeIndexs[2]);
		Vector d = mesh->getVertex(nodeIndexs[3]);

		Vector ab = b - a;
		Vector ac = c - a;
		Vector ad = d - a;

		return fabsf(ab*(ac%ad)) / 6.f;
	}

	void TetElement::getBodyVirtualWorks(double bodyForce[3], double *result) const{
		double volume = getVolume();
		double factor = volume / 4.0;
		constexpr int numNodesPerElement = 4;
		for (int i = 0; i < numNodesPerElement; i++){
			for (int j = 0; j < 3; j++)
				result[3 * i + j] = factor*bodyForce[j];
		}
	}

	float TetFacet::getArea() const{
		Vector a = mesh->getVertex(vertIndexs[0]);
		Vector b = mesh->getVertex(vertIndexs[1]);
		Vector c = mesh->getVertex(vertIndexs[2]);

		return ((b - a) % (c - a)).length() * 0.5f;
	}

	void TetFacet::getSurfVirtualWorks(double surfForce[3], double *result) const{
		double factor = getArea() / 3.0;
		constexpr int numVertPerSur = 3;
		for (int i = 0; i < numVertPerSur; i++){
			for (int j = 0; j < 3; j++)
				result[3 * i + j] = factor*surfForce[j];
		}
	}

	void LinearIsotropicTetElement::setBMatrixs(){
		constexpr int numNodesPerElement = 4;
		constexpr int entryCountPerB = 3;
		double dndx[numNodesPerElement], dndy[numNodesPerElement], dndz[numNodesPerElement];
		getShapeFunctionDerivatives(dndx, dndy, dndz);
		for (int i = 0; i < numNodesPerElement; i++){
			BMatrixs[entryCountPerB * i + 0] = dndx[i];
			BMatrixs[entryCountPerB * i + 1] = dndy[i];
			BMatrixs[entryCountPerB * i + 2] = dndz[i];
		}
	}

	void LinearIsotropicTetElement::generateSubStiffnessMatrix(int aNodeIndex, int bNodeIndex, const double *D, double *result) const{
		constexpr int entryCountPerB = 3;
		double volume = getVolume();
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
		constexpr int numNodesPerElement = 4;
		constexpr int entryCountPerB = 3;
		double dndx[numNodesPerElement], dndy[numNodesPerElement], dndz[numNodesPerElement];
		getShapeFunctionDerivatives(dndx, dndy, dndz);
		for (int i = 0; i < numNodesPerElement; i++){
			BMatrixs[entryCountPerB * i + 0] = dndx[i];
			BMatrixs[entryCountPerB * i + 1] = dndy[i];
			BMatrixs[entryCountPerB * i + 2] = dndz[i];
		}
	}

	void ReducedIsotropicTetElement::Intergration(const double *C, double *nlpart, double *nnpart) const{
		constexpr int entryCountPerB = 3;
		constexpr int numNodesPerElement = 4;
		constexpr int commonEntryNum = numNodesPerElement*numNodesPerElement*numNodesPerElement;
		memset(nlpart, 0, commonEntryNum * 3 * sizeof(double));
		memset(nnpart, 0, commonEntryNum * numNodesPerElement * sizeof(double));

		const double determinant = 6.0 * getVolume();
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
		double volume = getVolume();
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
		constexpr int numNodesPerElement = 4;
		constexpr int entryCountPerB = 3;
		double dndx[numNodesPerElement], dndy[numNodesPerElement], dndz[numNodesPerElement];
		getShapeFunctionDerivatives(dndx, dndy, dndz);
		for (int i = 0; i < numNodesPerElement; i++){
			BMatrixs[entryCountPerB * i + 0] = dndx[i];
			BMatrixs[entryCountPerB * i + 1] = dndy[i];
			BMatrixs[entryCountPerB * i + 2] = dndz[i];
		}
	}

	void LinearAnisortropicTetElement::generateSubStiffnessMatrix(int aNodeIndex, int bNodeIndex, const double *D, double *result) const{
		constexpr int entryCountPerB = 3;
		double volume = getVolume();
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
}