#include "stdafx.h"
#include "newmark.h"
#include "nodeIndexer.h"

using namespace ODER;

NewmarkIntergrator::NewmarkIntergrator(int totalDOFS, double beta, double gamma, const SparseMatrix &M, const SparseMatrix &K, int DOFS, double massDamp, double stiffDamp, double ts)
    :Intergrator(DOFS, massDamp, stiffDamp, ts){
	totalDofs = totalDOFS;
	betaDeltaT2 = beta*timeStep*timeStep;
	gammaDeltaT = gamma*timeStep;
	minusBetaDeltaT2 = timeStep*timeStep*0.5 - betaDeltaT2;
	minusGammaDeltaT = timeStep - gammaDeltaT;
	frequencies2 = allocAligned<double>(dofs);
	basises = new double[dofs*totalDofs];

	EigenSlover slover;
	slover.getEigenValVectors(dofs, M, K, frequencies2, basises);

#ifdef ODER_DEBUG
	for (int i = 0; i < dofs; i++){
		if (frequencies2[i] < 0.0)
			Severe("Stiffness Matrix is not semi-define");
	}
#endif
}

void NewmarkIntergrator::setExternalVirtualWork(const Forcer& forcer){
	forcer.getVirtualWorks(dofs, totalDofs, basises, externalVirtualWork);
	memcpy(a, externalVirtualWork, dofs*sizeof(double));
}

void NewmarkIntergrator::runOneTimeStep(){
	memcpy(pre_d, d, dofs*sizeof(double));
	memcpy(pre_v, v, dofs*sizeof(double));
	memcpy(pre_a, a, dofs*sizeof(double));

	for (int i = 0; i < dofs; i++){
		double d_pre = pre_d[i];
		double v_pre = pre_v[i];
		double a_pre = pre_a[i];
		double frequency2 = frequencies2[i];
		//predict d and v
		double d_predict = d_pre + timeStep*v_pre + minusBetaDeltaT2*a_pre;
		double v_predict = v_pre + minusGammaDeltaT*a_pre;
		//caculating a 
		double xi = massDamping + stiffnessDamping*frequency2;
		double righthand = externalVirtualWork[i] - xi*v_predict - frequency2*d_predict;
		double lefthand = 1.0 + gammaDeltaT*xi + betaDeltaT2*frequency2;
		a[i] = righthand / lefthand;
		//uptate d and v
		v[i] = v_predict + gammaDeltaT*a[i];
		d[i] = d_predict + betaDeltaT2*a[i];
	}
}

void NewmarkIntergrator::getDisplacements(const NodeIndexer &indexer, int displacementCount, double *displacements) const{
	memset(displacements, 0, displacementCount*sizeof(double));
	const int* constrainIndices = NULL;
	int constrainSize = indexer.getConstrainIndices(&constrainIndices);

	double *basis = basises;
	for (int i = 0; i < dofs; i++){
		int constrainIndex = 0;
		double displacement = d[i];
		for (int j = 0; j < displacementCount; j++){
			if (constrainIndex < constrainSize && j == constrainIndices[constrainIndex])
				constrainIndex++;
			else
				displacements[j] += displacement * basis[j - constrainIndex];
		}
		basis += totalDofs;
	}
}

void NewmarkIntergrator::getVertexPositions(const Reference<Mesh> &mesh, const NodeIndexer &indexer, int verticesCount, Vector *vertices, double *displacementBuffer) const{
	const int* constrainIndices = NULL;
	int constrainSize = indexer.getConstrainIndices(&constrainIndices);
	double *displacements = displacementBuffer;
	if (displacementBuffer == NULL)
		displacements = new double[totalDofs];

	getRawDisplacements(displacements);

	int constrainIndex = 0;
	int displacementIndex = 0;
	for (int vertIndex = 0; vertIndex < verticesCount; vertIndex++){
		for (int axis = 0; axis < 3; axis++){
			if (constrainIndex < constrainSize && (3 * vertIndex + axis) == constrainIndices[constrainIndex]){
				vertices[vertIndex][axis] = mesh->vertices[vertIndex][axis];
				constrainIndex++;
			}
			else
				vertices[vertIndex][axis] = mesh->vertices[vertIndex][axis] + float(displacements[displacementIndex++]);
		}
	}

	if (displacementBuffer == NULL)
	    delete[] displacements;
}

void NewmarkIntergrator::getRawDisplacements(double *displacements) const{
	memset(displacements, 0, totalDofs*sizeof(double));
	double *basis = basises;
	for (int i = 0; i < dofs; i++){
		double displacement = d[i];
		for (int j = 0; j < totalDofs; j++){
			displacements[j] += displacement * basis[j];
		}
		basis += totalDofs;
	}
}

NewmarkIntergrator::~NewmarkIntergrator(){
	freeAligned(frequencies2);
	delete[] basises;
}
