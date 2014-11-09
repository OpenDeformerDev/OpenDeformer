#include "stdafx.h"
#include "asymptoticNewmark.h"
#include "nodeIndexer.h"
#include "forcer.h"
#include "latool.h"

namespace ODER{
	AsymptoticNewmark::AsymptoticNewmark(int totalDOFS, double beta, double gamma, int DOFS, double massDamp, double stiffDamp, double ts,
		const Reference<Mesh> m, const Reference<NodeIndexer>& nodeIndexer, HyperelasticMaterial* mater)
		:Intergrator(mater->getNonlinearAsymptoticOrder() * DOFS, massDamp, stiffDamp, ts){
		mesh = m;
		indexer = nodeIndexer;
		material = mater;
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
		mater->generateMassMatrix(mesh, indexer, M);
		mater->generateStiffnessMatrix(mesh, indexer, K);

		EigenSlover slover;
		slover.getEigenValVectors(dofs, M, K, frequencies2, basises);

		material->preprocessWithReduction(mesh, indexer);
		loadFactors = allocAligned<double>(material->getNonlinearAsymptoticOrder());

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
		int matrixOrder = indexer->getMatrixOrder(mesh);
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
		loadFactors[0] = 1.0 / (factor + 1);

		for (int order = 1; order < orderCount; order++){
			material->getNodeForces(mesh, indexer, order, d, externalVirtualWork + order*matrixOrder);
			factor = 0.0;
			for (int i = 0; i < dofs; i++){
				double d_pre = pre_d[order*matrixOrder + i];
				double v_pre = pre_v[order*matrixOrder + i];
				double a_pre = pre_a[order*matrixOrder + i];
				double frequency2 = frequencies2[i];
				//predict d and v
				double d_predict = d_pre + timeStep*v_pre + minusBetaDeltaT2*a_pre;
				double v_predict = v_pre + minusGammaDeltaT*a_pre;
				//caculating a 
				double xi = massDamping + stiffnessDamping*frequency2;
				double righthand = externalVirtualWork[order*matrixOrder + i] - xi*v_predict - frequency2*d_predict;
				double lefthand = 1.0 + gammaDeltaT*xi + betaDeltaT2*frequency2;
				a[order*matrixOrder + i] = righthand / lefthand;
				//uptate d and v
				v[order*matrixOrder + i] = v_predict + gammaDeltaT*a[order*matrixOrder + i];
				d[order*matrixOrder + i] = d_predict + betaDeltaT2*a[order*matrixOrder + i];
				factor = d[order*matrixOrder + i] * d[i];
			}
			loadFactors[order] = -loadFactors[0] * factor;
			for (int i = 0; i < dofs; i++){
				a[order*matrixOrder + i] -= factor*a[i];
				v[order*matrixOrder + i] -= factor*v[i];
				d[order*matrixOrder + i] -= factor*d[i];
			}
		}
	}

	void AsymptoticNewmark::getRawDisplacements(double *displacements) const{
		int orderCount = material->getNonlinearAsymptoticOrder();
		const double *factors = loadFactors;
		double a = findRoot([orderCount, factors](double x)->double{
			double ret = 0.0, para = 1.0;
			for (int i = 0; i < orderCount; i++){
				ret += para*factors[i];
				para *= x;
			}
			return ret;
		}, 0.0, 2.0, 2e-8);

		memset(displacements, 0, totalDofs*sizeof(double));
		for (int order = 0; order < orderCount; order++){
			double *basis = basises;
			double factor = 1.0;
			for (int i = 0; i < dofs; i++){
				double displacement = d[i];
				for (int j = 0; j < totalDofs; j++){
					displacements[j] += factor * displacement * basis[j];
				}
				basis += totalDofs;
			}
			factor *= a;
		}
	}

	AsymptoticNewmark::~AsymptoticNewmark(){
		delete[] basises;
		delete material;
		freeAligned(frequencies2);
		freeAligned(loadFactors);
	}
}