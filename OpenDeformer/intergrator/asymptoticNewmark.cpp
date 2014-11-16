#include "stdafx.h"
#include "asymptoticNewmark.h"
#include "nodeIndexer.h"
#include "forcer.h"
#include "latool.h"

namespace ODER{
	AsymptoticNewmark::AsymptoticNewmark(double beta, double gamma, int DOFS, double massDamp, double stiffDamp, double ts,
		const Reference<Mesh> m, const Reference<NodeIndexer>& nodeIndexer, HyperelasticMaterial* mater)
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

		EigenSlover slover;
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

		//calculate load factor 0
		loadFactors[0] = 1.0 / sqrt(factor + 1);
		setFullDisplacement(0, loadFactors[0]);

		for (int order = 1; order < orderCount; order++){
			material->getNodeForces(mesh, indexer, order, totalDofs, fullDisplacements, vwBuffer);
			setReducedVirtualWorks(order);
			factor = 0.0;
			int offset = order*dofs;
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
			setFullDisplacement(order, loadFactors[order]);
		}
	}

	void AsymptoticNewmark::getRawDisplacements(double *displacements) const{
		int orderCount = material->getNonlinearAsymptoticOrder();
		const double *factors = loadFactors;
		double a = findRoot([orderCount, factors](double x)->double{
			double ret = 0.0, para = 1.0;
			for (int i = 0; i < orderCount; i++){
				para *= x;
				ret += para*factors[i];
			}
			return ret - 1.0;
		}, 0.0, 2.0, 2e-6);
		memset(displacements, 0, totalDofs*sizeof(double));

		const double *displacementPerOrder = fullDisplacements;
		double factor = 1.0;
		for (int order = 0; order < orderCount; order++){
			factor *= a;
			for (int i = 0; i < totalDofs; i++){
				displacements[i] += factor*displacementPerOrder[i];
			}
			displacementPerOrder += totalDofs;
		}
	}

	void AsymptoticNewmark::setFullDisplacement(int order, double loadFactor){
		const double *reduced = d + order*dofs;
		double *displacement = fullDisplacements + order*totalDofs;
		memset(displacement, 0, totalDofs * sizeof(double));
		if (order == 0){
			const double *basis = basises;
			for (int i = 0; i < dofs; i++){
				double entry = loadFactor*reduced[i];
				for (int j = 0; j < totalDofs; j++){
					displacement[j] += entry * basis[j];
				}
				basis += totalDofs;
			}
		}
		else{
			const double *basis = basises;
			for (int i = 0; i < dofs; i++){
				double entry = loadFactor * d[i] + reduced[i];
				for (int j = 0; j < totalDofs; j++){
					displacement[j] += entry * basis[j];
				}
				basis += totalDofs;
			}
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
			basis += totalDofs;
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