#include "stdafx.h"
#include "tetelement.h"

using namespace ODER;

TetElement::TetElement(TetMesh *m, int n, MarterialType t) :Element(m, n, t){
	if (matchMaterialFlag(type, Marterial_Linear)){
		const int entrys = 4 * 3;
		BMatrixs = new double[entrys];
	}
	else{
		Severe("this feature in TetElement Constructor is not implemented");
	}
}

TetElement::TetElement(TetMesh *m, MarterialType t) : Element(m, t){
	if (matchMaterialFlag(type, Marterial_Linear)){
		const int entrys = 4 * 3;
		BMatrixs = new double[entrys];
	}
	else{
		Severe("this feature in TetElement Constructor is not implemented");
	}
}

void TetElement::getShapeFunctionDerivatives(double *dndx, double *dndy, double *dndz) const{
	Vector a = mesh->vertices[nodeIndexs[0]];
	Vector b = mesh->vertices[nodeIndexs[1]];
	Vector c = mesh->vertices[nodeIndexs[2]];
	Vector d = mesh->vertices[nodeIndexs[3]];

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

void TetElement::setBMatrixs(){
	if (matchMaterialFlag(type, Marterial_Linear)){
		const int numNodesPerElement = 4;
		const int entryCountPerB = 3;
		double dndx[numNodesPerElement], dndy[numNodesPerElement], dndz[numNodesPerElement];
		getShapeFunctionDerivatives(dndx, dndy, dndz);
		for (int i = 0; i < numNodesPerElement; i++){
			BMatrixs[entryCountPerB * i + 0] = dndx[i];
			BMatrixs[entryCountPerB * i + 1] = dndy[i];
			BMatrixs[entryCountPerB * i + 2] = dndz[i];
		}
	}
	else{
		Severe("this feature in TetElement::setBMatrixs is not implemented");
	}
}

void TetElement::generateSubStiffnessMatrix(int aNodeIndex, int bNodeIndex, const double *D, double *result) const{
	if (matchMaterialFlag(type, MarterialType(Marterial_Linear | Marterial_Isotropic))){
		const int entryCountPerB = 3;
		const double *Ba = BMatrixs + entryCountPerB*aNodeIndex;
		const double *Bb = BMatrixs + entryCountPerB*bNodeIndex;

		double factor = 1.0 / (36.0*volume);
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
	else{
		Severe("this feature in TetElement::generateSubStiffnessMatrix is not implemented");
	}
}

void TetElement::setVolume(){
	Vector a = mesh->vertices[nodeIndexs[0]];
	Vector b = mesh->vertices[nodeIndexs[1]];
	Vector c = mesh->vertices[nodeIndexs[2]];
	Vector d = mesh->vertices[nodeIndexs[3]];

	Vector ab = b - a;
	Vector ac = c - a;
	Vector ad = d - a;

	volume = fabsf(ab*(ac%ad)) / 6.f;
}

void TetElement::getBodyVirtualWorks(double bodyForce[3], double *result) const{
	double factor = volume / 4.0;
	const int numNodesPerElement = 4;
	for (int i = 0; i < numNodesPerElement; i++){
		for (int j = 0; j < 3; j++)
			result[3 * i + j] = factor*bodyForce[j];
	}
}

float TetFacet::getArea() const{
	Vector a = mesh->vertices[vertIndexs[0]];
	Vector b = mesh->vertices[vertIndexs[1]];
	Vector c = mesh->vertices[vertIndexs[2]];

	return ((b - a) % (c - a)).length() * 0.5f;
}

void TetFacet::getSurfVirtualWorks(double surfForce[3], double *result) const{
	double factor = getArea() / 3.0;
	const int numVertPerSur = 3;
	for (int i = 0; i < numVertPerSur; i++){
		for (int j = 0; j < 3; j++)
			result[3 * i + j] = factor*surfForce[j];
	}
}