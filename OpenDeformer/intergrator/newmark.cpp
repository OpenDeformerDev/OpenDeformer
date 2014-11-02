#include "stdafx.h"
#include "newmark.h"
#include "nodeIndexer.h"
#include "forcer.h"

namespace ODER{
	NewmarkIntergrator::NewmarkIntergrator(int totalDOFS, double beta, double gamma, int DOFS, double massDamp, double stiffDamp, double ts,
		const Reference<Mesh> m, const Reference<NodeIndexer>& nodeIndexer, const Reference<MecMaterial> &mater)
		:Intergrator(DOFS, massDamp, stiffDamp, ts, m, nodeIndexer, mater){
		totalDofs = totalDOFS;
		betaDeltaT2 = beta*timeStep*timeStep;
		gammaDeltaT = gamma*timeStep;
		minusBetaDeltaT2 = timeStep*timeStep*0.5 - betaDeltaT2;
		minusGammaDeltaT = timeStep - gammaDeltaT;
		frequencies2 = allocAligned<double>(dofs);
		basises = new double[dofs*totalDofs];

		int matrixOrder = indexer->getMatrixOrder(mesh);
		SparseMatrixAssembler M(matrixOrder);
		SparseMatrixAssembler K(matrixOrder);
		material->generateMassMatrix(mesh, indexer, M);
		material->generateStiffnessMatrix(mesh, indexer, K);

		EigenSlover slover;
		slover.getEigenValVectors(dofs, M, K, frequencies2, basises);

		loadFactors = NULL;
		if (matchMaterialFlag(material->getMaterialType(), Marterial_NonLinear))
			loadFactors = allocAligned<double>(material->getNonlinearAsymptoticOrder());


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

	void NewmarkIntergrator::getDisplacements(int displacementCount, double *displacements) const{
		memset(displacements, 0, displacementCount*sizeof(double));
		const int* constrainIndices = NULL;
		int constrainSize = indexer->getConstrainIndices(&constrainIndices);

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

	void NewmarkIntergrator::getVertexPositions(Vector *vertices, double *displacementBuffer) const{
		const int* constrainIndices = NULL;
		int constrainSize = indexer->getConstrainIndices(&constrainIndices);
		double *displacements = displacementBuffer;
		if (displacementBuffer == NULL)
			displacements = new double[totalDofs];

		getRawDisplacements(displacements);

		int constrainIndex = 0;
		int displacementIndex = 0;
		int vertCount = mesh->numNodes;
		for (int vertIndex = 0; vertIndex < vertCount; vertIndex++){
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
		if (loadFactors)
			freeAligned(loadFactors);
	}
}
