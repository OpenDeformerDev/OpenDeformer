#include "stdafx.h"
#include "asymptoticNewmark.h"
#include "nodeIndexer.h"
#include "forcer.h"
#include "numerMethod.h"

namespace ODER{
	AsymptoticNewmark::AsymptoticNewmark(Scalar beta, Scalar gamma, int DOFS, Scalar massDamp, Scalar stiffDamp, Scalar ts,
		const Reference<Mesh>& m, const Reference<NodeIndexer>& nodeIndexer, ReducedHyperelasticMaterial* mater)
		:Intergrator(DOFS, massDamp, stiffDamp, ts){
		int orderCount = mater->getNonlinearAsymptoticOrder();
		int entrys = orderCount*dofs;
		d = allocAligned<Scalar>(entrys);
		v = allocAligned<Scalar>(entrys);
		a = allocAligned<Scalar>(entrys);
		pre_d = allocAligned<Scalar>(entrys);
		pre_v = allocAligned<Scalar>(entrys);
		pre_a = allocAligned<Scalar>(entrys);
		externalVirtualWork = allocAligned<Scalar>(entrys);

		memset(d, 0, entrys*sizeof(Scalar));
		memset(v, 0, entrys*sizeof(Scalar));
		memset(a, 0, entrys*sizeof(Scalar));
		memset(pre_d, 0, entrys*sizeof(Scalar));
		memset(pre_v, 0, entrys*sizeof(Scalar));
		memset(pre_a, 0, entrys*sizeof(Scalar));
		memset(externalVirtualWork, 0, entrys*sizeof(Scalar));

		mesh = m;
		indexer = nodeIndexer;
		material = mater;
		betaDeltaT2 = beta*timeStep*timeStep;
		gammaDeltaT = gamma*timeStep;
		minusBetaDeltaT2 = timeStep*timeStep*Scalar(0.5) - betaDeltaT2;
		minusGammaDeltaT = timeStep - gammaDeltaT;
		totalDofs = indexer->getMatrixOrder(mesh);
		frequencies2 = allocAligned<Scalar>(dofs);
		basises = allocAligned<Scalar>(dofs*totalDofs);

		SparseMatrixAssembler M(totalDofs);
		SparseMatrixAssembler K(totalDofs);
		material->generateMassMatrix(mesh, indexer, M);
		material->generateStiffnessMatrix(mesh, indexer, K); 

		EigenSolver slover;
		slover.getEigenValVectors(dofs, M, K, frequencies2, basises);

		material->preprocessWithReduction(mesh, indexer);
		loadFactors = allocAligned<Scalar>(orderCount);
		fullDisplacements = allocAligned<Scalar>(orderCount*totalDofs);
		vwBuffer = allocAligned<Scalar>(totalDofs);

#ifdef ODER_DEBUG
		for (int i = 0; i < dofs; i++){
			if (frequencies2[i] < 0.0)
				Severe("Stiffness Matrix is not semi-define");
		}
#endif
	}

	void AsymptoticNewmark::setExternalVirtualWork(const Forcer& forcer){
		forcer.getVirtualWorks(dofs, totalDofs, basises, externalVirtualWork);
		memcpy(a, externalVirtualWork, dofs*sizeof(Scalar));
	}

	void AsymptoticNewmark::runOneTimeStep(){
		int orderCount = material->getNonlinearAsymptoticOrder();
		memset(fullDisplacements, 0, orderCount*totalDofs*sizeof(Scalar));
		memcpy(pre_d, d, orderCount*dofs*sizeof(Scalar));
		memcpy(pre_v, v, orderCount*dofs*sizeof(Scalar));
		memcpy(pre_a, a, orderCount*dofs*sizeof(Scalar));

		//order 1
		Scalar factor = 0.0;
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
			factor += d[i] * d[i];
		}

		//calculate load factor 1
		Scalar loadFactor0 = Scalar(1.0) / sqrt(Scalar(1.0) + factor);
		loadFactors[0] = loadFactor0;
		//set full order 1 displacement
		const Scalar *basis = basises;
		for (int i = 0; i < dofs; i++){
			Scalar entry = d[i] * loadFactor0;
			for (int j = 0; j < totalDofs; j++){
				fullDisplacements[j] += entry * basis[j];
			}
			basis += totalDofs;
		}

		Scalar loadFactor0_second = loadFactor0*loadFactor0;
		for (int order = 1; order < orderCount; order++){
			material->getNodeForces(mesh, indexer, order, totalDofs, fullDisplacements, vwBuffer);
			//set reduced virtualworks
			int offset = order*dofs;
			const Scalar *basis = basises;
			Scalar *reducedVW = externalVirtualWork + order*dofs;
			for (int i = 0; i < dofs; i++){
				Scalar entry = 0.0;
				for (int j = 0; j < totalDofs; j++)
					entry += basis[j] * vwBuffer[j];
				externalVirtualWork[offset + i] = entry;
				basis += totalDofs;
			}

			factor = 0.0;
			for (int i = 0; i < dofs; i++){
				Scalar d_pre = pre_d[offset + i];
				Scalar v_pre = pre_v[offset + i];
				Scalar a_pre = pre_a[offset + i];
				Scalar frequency2 = frequencies2[i];
				//predict d and v
				Scalar d_predict = d_pre + timeStep*v_pre + minusBetaDeltaT2*a_pre;
				Scalar v_predict = v_pre + minusGammaDeltaT*a_pre;
				//caculating a 
				Scalar xi = massDamping + stiffnessDamping*frequency2;
				Scalar righthand = externalVirtualWork[offset + i] - xi*v_predict - frequency2*d_predict;
				Scalar lefthand = Scalar(1.0) + gammaDeltaT*xi + betaDeltaT2*frequency2;
				a[offset + i] = righthand / lefthand;
				//uptate d and v
				v[offset + i] = v_predict + gammaDeltaT*a[offset + i];
				d[offset + i] = d_predict + betaDeltaT2*a[offset + i];
				factor += d[offset + i] * d[i];
			}
			Scalar loadFactorPerOrder = -loadFactor0_second * factor;
			loadFactors[order] = loadFactorPerOrder;
			//set full displacement per order
			Scalar *displacement = fullDisplacements + order*totalDofs;
			basis = basises;
			for (int i = 0; i < dofs; i++){
				Scalar entry = d[offset + i] + loadFactorPerOrder * d[i];
				for (int j = 0; j < totalDofs; j++){
					displacement[j] += entry * basis[j];
				}
				basis += totalDofs;
			}
		}
	}

	void AsymptoticNewmark::getRawDisplacements(Scalar *displacements) const{
		int orderCount = material->getNonlinearAsymptoticOrder();
		const Scalar *factors = loadFactors;
		Scalar root = findRoot([orderCount, factors](Scalar x)->Scalar{
			Scalar ret = 0.0, para = x;
			for (int i = 0; i < orderCount; i++){
				ret += para*factors[i];
				para *= x;
			}
			return ret - Scalar(1.0);
		}, Scalar(0.0), Scalar(2.0), Scalar(2e-6));

		memset(displacements, 0, totalDofs*sizeof(Scalar));

		const Scalar *displacementPerOrder = fullDisplacements;
		Scalar factor = root;
		for (int order = 0; order < orderCount; order++){
			for (int i = 0; i < totalDofs; i++){
				displacements[i] += factor*displacementPerOrder[i];
			}
			displacementPerOrder += totalDofs;
			factor *= root;
		}
	}

	void AsymptoticNewmark::updateMeshVerticesDisplacements(const Reference<NodeIndexer> &indexer, Reference<Mesh> &mesh) const {
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

	AsymptoticNewmark::~AsymptoticNewmark(){
		freeAligned(a);
		freeAligned(d);
		freeAligned(v);
		freeAligned(pre_a);
		freeAligned(pre_d);
		freeAligned(pre_v);
		freeAligned(externalVirtualWork);
		freeAligned(frequencies2);

		delete material;
		freeAligned(basises);
		freeAligned(frequencies2);
		freeAligned(loadFactors);
		freeAligned(fullDisplacements);
		freeAligned(vwBuffer);
	}
}