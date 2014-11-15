#include "stdafx.h"
#include "asymptoticNewmark.h"
#include "nodeIndexer.h"
#include "forcer.h"
#include "latool.h"

namespace ODER{
	AsymptoticNewmark::AsymptoticNewmark(double beta, double gamma, int DOFS, double massDamp, double stiffDamp, double ts,
		const Reference<Mesh> m, const Reference<NodeIndexer>& nodeIndexer, HyperelasticMaterial* mater)
		:Intergrator(mater->getNonlinearAsymptoticOrder() * DOFS, massDamp, stiffDamp, ts){
		mesh = m;
		indexer = nodeIndexer;
		material = mater;
		betaDeltaT2 = beta*timeStep*timeStep;
		gammaDeltaT = gamma*timeStep;
		minusBetaDeltaT2 = timeStep*timeStep*0.5 - betaDeltaT2;
		minusGammaDeltaT = timeStep - gammaDeltaT;
		totalDofs = indexer->getMatrixOrder(mesh);
		frequencies2 = allocAligned<double>(dofs);
		basises = new double[dofs*totalDofs];

		SparseMatrixAssembler M(totalDofs);
		SparseMatrixAssembler K(totalDofs);
		mater->generateMassMatrix(mesh, indexer, M);
		mater->generateStiffnessMatrix(mesh, indexer, K);

		EigenSlover slover;
		slover.getEigenValVectors(dofs, M, K, frequencies2, basises);

		int orderCount = material->getNonlinearAsymptoticOrder();
		material->preprocessWithReduction(mesh, indexer);
		loadFactors = allocAligned<double>(orderCount);
		fullDisplacements = allocAligned<double>(orderCount*totalDofs);
		vwBuffer = allocAligned<double>(totalDofs);

#ifdef ODER_DEBUG
		for (int i = 0; i < dofs; i++){
			if (frequencies2[i] < 0.0)
				Severe("Stiffness Matrix is not semi-define");
		}
#endif
	}

	void AsymptoticNewmark::setExternalVirtualWork(const Forcer& forcer){
		forcer.getVirtualWorks(dofs, totalDofs, basises, externalVirtualWork);
		memcpy(a, externalVirtualWork, dofs*sizeof(double));
	}

	void AsymptoticNewmark::runOneTimeStep(){
		int orderCount = material->getNonlinearAsymptoticOrder();
		memcpy(pre_d, d, orderCount*dofs*sizeof(double));
		memcpy(pre_v, v, orderCount*dofs*sizeof(double));
		memcpy(pre_a, a, orderCount*dofs*sizeof(double));

		//order 0
		double factor = 0.0;
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
			factor += d[i] * d[i];
		}
		setFullDisplacement(0);
		//calculate load factor 0
		loadFactors[0] = 1.0 / (factor + 1);

		for (int order = 1; order < orderCount; order++){
			material->getNodeForces(mesh, indexer, order, fullDisplacements, vwBuffer);
			setReducedVirtualWorks(order);
			factor = 0.0;
			int offset = order*totalDofs;
			for (int i = 0; i < dofs; i++){
				double d_pre = pre_d[offset + i];
				double v_pre = pre_v[offset + i];
				double a_pre = pre_a[offset + i];
				double frequency2 = frequencies2[i];
				//predict d and v
				double d_predict = d_pre + timeStep*v_pre + minusBetaDeltaT2*a_pre;
				double v_predict = v_pre + minusGammaDeltaT*a_pre;
				//caculating a 
				double xi = massDamping + stiffnessDamping*frequency2;
				double righthand = externalVirtualWork[offset + i] - xi*v_predict - frequency2*d_predict;
				double lefthand = 1.0 + gammaDeltaT*xi + betaDeltaT2*frequency2;
				a[offset + i] = righthand / lefthand;
				//uptate d and v
				v[offset + i] = v_predict + gammaDeltaT*a[offset + i];
				d[offset + i] = d_predict + betaDeltaT2*a[offset + i];
				factor = d[offset + i] * d[i];
			}
			loadFactors[order] = -loadFactors[0] * factor;
			for (int i = 0; i < dofs; i++){
				a[offset + i] -= factor*a[i];
				v[offset + i] -= factor*v[i];
				d[offset + i] -= factor*d[i];
			}
			setFullDisplacement(order);
		}
	}

	void AsymptoticNewmark::getRawDisplacements(double *displacements) const{
		int orderCount = material->getNonlinearAsymptoticOrder();
		const double *factors = loadFactors;
		double a = findRoot([orderCount, factors](double x)->double{
			double ret = 0.0, para = x;
			for (int i = 0; i < orderCount; i++){
				ret += para*factors[i];
				para *= x;
			}
			return ret;
		}, 0.0, 2.0, 2e-6);

		memset(displacements, 0, totalDofs*sizeof(double));
		const double *displacementPerOrder = fullDisplacements;
		double factor = 1.0;
		for (int order = 0; order < orderCount; order++){
			for (int i = 0; i < totalDofs; i++){
				displacements[i] = factor*displacementPerOrder[i];
			}
			displacementPerOrder += totalDofs;
			factor *= a;
		}
	}

	void AsymptoticNewmark::setFullDisplacement(int order){
		const double *reduced = d + order*dofs;
		double *displacement = fullDisplacements + order*totalDofs;
		const double *basis = basises;
		for (int i = 0; i < dofs; i++){
			double entry = reduced[i];
			for (int j = 0; j < totalDofs; j++){
				displacement[j] += entry * basis[j];
			}
			basis += totalDofs;
		}
	}

	void AsymptoticNewmark::setReducedVirtualWorks(int order){
		const double *basis = basises;
		double *reducedVW = externalVirtualWork + order*dofs;
		for (int i = 0; i < dofs; i++){
			double entry = 0.0;
			for (int j = 0; j < totalDofs; j++)
				entry += basis[j] * vwBuffer[j];
			reducedVW[i] = entry;
		}
	}

	AsymptoticNewmark::~AsymptoticNewmark(){
		delete[] basises;
		delete material;
		freeAligned(frequencies2);
		freeAligned(loadFactors);
		freeAligned(fullDisplacements);
	}
}