#include "stdafx.h"
#include "linearNewmark.h"
#include "nodeIndexer.h"
#include "forcer.h"

namespace ODER{
	LinearNewmark::LinearNewmark(double beta, double gamma, int DOFS, double massDamp, double stiffDamp, double ts,
		const Reference<Mesh> m, const Reference<NodeIndexer>& nodeIndexer, const MecMaterial& mater)
		:Intergrator(DOFS, massDamp, stiffDamp, ts){
		betaDeltaT2 = beta*timeStep*timeStep;
		gammaDeltaT = gamma*timeStep;
		minusBetaDeltaT2 = timeStep*timeStep*0.5 - betaDeltaT2;
		minusGammaDeltaT = timeStep - gammaDeltaT;
		totalDofs = nodeIndexer->getMatrixOrder(m);
		frequencies2 = allocAligned<double>(dofs);
		basises = new double[dofs*totalDofs];

		SparseMatrixAssembler M(totalDofs);
		SparseMatrixAssembler K(totalDofs);
		mater.generateMassMatrix(m, nodeIndexer, M);
		mater.generateStiffnessMatrix(m, nodeIndexer, K);

		EigenSlover slover;
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
		memcpy(a, externalVirtualWork, dofs*sizeof(double));
	}

	void LinearNewmark::runOneTimeStep(){
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

	void LinearNewmark::getRawDisplacements(double *displacements) const{
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

	LinearNewmark::~LinearNewmark(){
		freeAligned(frequencies2);
		delete[] basises;
	}
}
