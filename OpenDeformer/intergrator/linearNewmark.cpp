#include "stdafx.h"
#include "linearNewmark.h"
#include "nodeIndexer.h"
#include "mechMaterial.h"
#include "forcer.h"

namespace ODER{
	LinearNewmark::LinearNewmark(Scalar beta, Scalar gamma, int DOFS, Scalar massDamp, Scalar stiffDamp, Scalar ts,
		const Reference<Mesh>& m, const Reference<NodeIndexer>& nodeIndexer, const HookeMaterial& mater)
		:Intergrator(DOFS, massDamp, stiffDamp, ts){
		d = allocAligned<Scalar>(dofs);
		v = allocAligned<Scalar>(dofs);
		a = allocAligned<Scalar>(dofs);
		pre_d = allocAligned<Scalar>(dofs);
		pre_v = allocAligned<Scalar>(dofs);
		pre_a = allocAligned<Scalar>(dofs);
		externalVirtualWork = allocAligned<Scalar>(dofs);

		memset(d, 0, dofs*sizeof(Scalar));
		memset(v, 0, dofs*sizeof(Scalar));
		memset(a, 0, dofs*sizeof(Scalar));
		memset(pre_d, 0, dofs*sizeof(Scalar));
		memset(pre_v, 0, dofs*sizeof(Scalar));
		memset(pre_a, 0, dofs*sizeof(Scalar));
		memset(externalVirtualWork, 0, dofs*sizeof(Scalar));

		betaDeltaT2 = beta*timeStep*timeStep;
		gammaDeltaT = gamma*timeStep;
		minusBetaDeltaT2 = timeStep*timeStep*Scalar(0.5) - betaDeltaT2;
		minusGammaDeltaT = timeStep - gammaDeltaT;
		totalDofs = nodeIndexer->getMatrixOrder(m);
		frequencies2 = allocAligned<Scalar>(dofs);
		basises = new Scalar[dofs*totalDofs];

		SparseMatrixAssembler M(totalDofs);
		SparseMatrixAssembler K(totalDofs);
		mater.generateMassMatrix(m, nodeIndexer, M);
		mater.generateStiffnessMatrix(m, nodeIndexer, K);

		EigenSolver slover;
		slover.getEigenValVectors(dofs, M, K, frequencies2, basises);


#ifdef ODER_DEBUG
		for (int i = 0; i < dofs; i++){
			if (frequencies2[i] < 0.0)
				Severe("Stiffness Matrix is not semi-define");
		}
#endif
	}

	void LinearNewmark::setExternalVirtualWork(const Forcer& forcer){
		forcer.getVirtualWorks(dofs, totalDofs, basises, externalVirtualWork);
		memcpy(a, externalVirtualWork, dofs*sizeof(Scalar));
	}

	void LinearNewmark::runOneTimeStep(){
		memcpy(pre_d, d, dofs*sizeof(Scalar));
		memcpy(pre_v, v, dofs*sizeof(Scalar));
		memcpy(pre_a, a, dofs*sizeof(Scalar));

		for (int i = 0; i < dofs; i++){
			Scalar d_pre = pre_d[i];
			Scalar v_pre = pre_v[i];
			Scalar a_pre = pre_a[i];
			Scalar frequency2 = frequencies2[i];
			//predict d and v
			Scalar d_predict = d_pre + timeStep*v_pre + minusBetaDeltaT2*a_pre;
			Scalar v_predict = v_pre + minusGammaDeltaT*a_pre;
			//caculating a 
			Scalar xi = massDamping + stiffnessDamping*frequency2;
			Scalar righthand = externalVirtualWork[i] - xi*v_predict - frequency2*d_predict;
			Scalar lefthand = Scalar(1.0) + gammaDeltaT*xi + betaDeltaT2*frequency2;
			a[i] = righthand / lefthand;
			//uptate d and v
			v[i] = v_predict + gammaDeltaT*a[i];
			d[i] = d_predict + betaDeltaT2*a[i];
		}
	}

	void LinearNewmark::getRawDisplacements(Scalar *displacements) const{
		memset(displacements, 0, totalDofs*sizeof(Scalar));
		Scalar *basis = basises;
		for (int i = 0; i < dofs; i++){
			Scalar displacement = d[i];
			for (int j = 0; j < totalDofs; j++){
				displacements[j] += displacement * basis[j];
			}
			basis += totalDofs;
		}
	}

	void LinearNewmark::updateMeshVerticesDisplacements(const Reference<NodeIndexer> &indexer, Reference<Mesh> &mesh) const {
		Scalar *displacements = new Scalar[totalDofs];
		getRawDisplacements(displacements);

		auto constrainIter = indexer->getConstrainIterBegin();
		auto constrainEnd = indexer->getConstrainIterEnd();
		int displacementIndex = 0;
		int vertCount = mesh->getNodeCount();
		for (int vertIndex = 0; vertIndex < vertCount; vertIndex++) {
			for (int axis = 0; axis < 3; axis++) {
				if (constrainIter != constrainEnd && (3 * vertIndex + axis) == *constrainIter) 
					constrainIter++;
				else
					mesh->getVertexDisplacement(vertIndex)[axis] = displacements[displacementIndex++];
			}
		}

		delete[] displacements;
	}

	LinearNewmark::~LinearNewmark(){
		freeAligned(a);
		freeAligned(d);
		freeAligned(v);
		freeAligned(pre_a);
		freeAligned(pre_d);
		freeAligned(pre_v);
		freeAligned(externalVirtualWork);
		freeAligned(frequencies2);
		delete[] basises;
	}
}
