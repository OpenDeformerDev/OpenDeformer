#include "stdafx.h"
#include "asymptoticNewmark.h"
#include "nodeIndexer.h"
#include "forcer.h"
#include "numerMethod.h"

namespace ODER{
	AsymptoticNewmark::AsymptoticNewmark(double beta, double gamma, int DOFS, double massDamp, double stiffDamp, double ts,
		const Reference<Mesh> m, const Reference<NodeIndexer>& nodeIndexer, ReducedHyperelasticMaterial* mater)
		:Intergrator(DOFS, massDamp, stiffDamp, ts){
		int orderCount = mater->getNonlinearAsymptoticOrder();
		int entrys = orderCount*dofs;
		d = allocAligned<double>(entrys);
		v = allocAligned<double>(entrys);
		a = allocAligned<double>(entrys);
		pre_d = allocAligned<double>(entrys);
		pre_v = allocAligned<double>(entrys);
		pre_a = allocAligned<double>(entrys);
		externalVirtualWork = allocAligned<double>(entrys);

		memset(d, 0, entrys*sizeof(double));
		memset(v, 0, entrys*sizeof(double));
		memset(a, 0, entrys*sizeof(double));
		memset(pre_d, 0, entrys*sizeof(double));
		memset(pre_v, 0, entrys*sizeof(double));
		memset(pre_a, 0, entrys*sizeof(double));
		memset(externalVirtualWork, 0, entrys*sizeof(double));

		mesh = m;
		indexer = nodeIndexer;
		material = mater;
		betaDeltaT2 = beta*timeStep*timeStep;
		gammaDeltaT = gamma*timeStep;
		minusBetaDeltaT2 = timeStep*timeStep*0.5 - betaDeltaT2;
		minusGammaDeltaT = timeStep - gammaDeltaT;
		totalDofs = indexer->getMatrixOrder(mesh);
		frequencies2 = allocAligned<double>(dofs);
		basises = allocAligned<double>(dofs*totalDofs);

		SparseMatrixAssembler M(totalDofs);
		SparseMatrixAssembler K(totalDofs);
		material->generateMassMatrix(mesh, indexer, M);
		material->generateStiffnessMatrix(mesh, indexer, K); 

		EigenSolver slover;
		slover.getEigenValVectors(dofs, M, K, frequencies2, basises);

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
		memset(fullDisplacements, 0, orderCount*totalDofs*sizeof(double));
		memcpy(pre_d, d, orderCount*dofs*sizeof(double));
		memcpy(pre_v, v, orderCount*dofs*sizeof(double));
		memcpy(pre_a, a, orderCount*dofs*sizeof(double));

		//order 1
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

		//calculate load factor 1
		double loadFactor0 = 1.0 / sqrt(1.0 + factor);
		loadFactors[0] = loadFactor0;
		//set full order 1 displacement
		const double *basis = basises;
		for (int i = 0; i < dofs; i++){
			double entry = d[i] * loadFactor0;
			for (int j = 0; j < totalDofs; j++){
				fullDisplacements[j] += entry * basis[j];
			}
			basis += totalDofs;
		}

		double loadFactor0_second = loadFactor0*loadFactor0;
		for (int order = 1; order < orderCount; order++){
			material->getNodeForces(mesh, indexer, order, totalDofs, fullDisplacements, vwBuffer);
			//set reduced virtualworks
			int offset = order*dofs;
			const double *basis = basises;
			double *reducedVW = externalVirtualWork + order*dofs;
			for (int i = 0; i < dofs; i++){
				double entry = 0.0;
				for (int j = 0; j < totalDofs; j++)
					entry += basis[j] * vwBuffer[j];
				externalVirtualWork[offset + i] = entry;
				basis += totalDofs;
			}

			factor = 0.0;
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
				factor += d[offset + i] * d[i];
			}
			double loadFactorPerOrder = -loadFactor0_second * factor;
			loadFactors[order] = loadFactorPerOrder;
			//set full displacement per order
			double *displacement = fullDisplacements + order*totalDofs;
			basis = basises;
			for (int i = 0; i < dofs; i++){
				double entry = d[offset + i] + loadFactorPerOrder * d[i];
				for (int j = 0; j < totalDofs; j++){
					displacement[j] += entry * basis[j];
				}
				basis += totalDofs;
			}
		}
	}

	void AsymptoticNewmark::getRawDisplacements(double *displacements) const{
		int orderCount = material->getNonlinearAsymptoticOrder();
		const double *factors = loadFactors;
		double root = findRoot([orderCount, factors](double x)->double{
			double ret = 0.0, para = x;
			for (int i = 0; i < orderCount; i++){
				ret += para*factors[i];
				para *= x;
			}
			return ret - 1.0;
		}, 0.0, 2.0, 2e-6);

		memset(displacements, 0, totalDofs*sizeof(double));

		const double *displacementPerOrder = fullDisplacements;
		double factor = root;
		for (int order = 0; order < orderCount; order++){
			for (int i = 0; i < totalDofs; i++){
				displacements[i] += factor*displacementPerOrder[i];
			}
			displacementPerOrder += totalDofs;
			factor *= root;
		}
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