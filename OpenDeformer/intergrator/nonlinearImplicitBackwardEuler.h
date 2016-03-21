#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_INTERGRATOR_NONLINEARIMPLICITBACKWARDEULER_H
#define ODER_INTERGRATOR_NONLINEARIMPLICITBACKWARDEULER_H

#include "oder.h"
#include "intergrator.h"
#include "fullOrderNonlinear.h"
#include "sparseMatrix.h"
#include "linearSolver.h"
#include "forcer.h"
#include "numerMethod.h"

namespace ODER {
	template<class SpMatrix> class NonlinearImplicitBackwardEuler : public Intergrator {
	public:
		NonlinearImplicitBackwardEuler(int DOFS, double massDamp, double stiffDamp, double ts,
			const Reference<Mesh>& m, const Reference<NodeIndexer>& nodeIndexer, const FullOrderNonlinearMaterial<SpMatrix> * mater, LinearSolver<SpMatrix>* linearSolver);
		void setExternalVirtualWork(const Forcer& forcer);
		void runOneTimeStep();
		void getRawDisplacements(double *displacements) const;
		~NonlinearImplicitBackwardEuler();
	private:
		double *d;
		double *v;
		double *pre_v;
		double *internalVirtualWork;
		double *externalVirtualWork;
		double *rhs;

		double *memory;
		Reference<Mesh> mesh;
		Reference<NodeIndexer> indexer;
		const FullOrderNonlinearMaterial<SpMatrix> *material;
		LinearSolver<SpMatrix> *solver;
		SpMatrix *tagentMatrix;
		SpMatrix *massMatrix;
		std::vector<std::unordered_map<int, int>> matrixIndices;
	};

	template<class SpMatrix> NonlinearImplicitBackwardEuler<SpMatrix>::NonlinearImplicitBackwardEuler(int DOFS, double massDamp, double stiffDamp, double ts,
		const Reference<Mesh>& m, const Reference<NodeIndexer>& nodeIndexer, const FullOrderNonlinearMaterial<SpMatrix> *mater, LinearSolver<SpMatrix> *linearSolver) 
		: Intergrator(DOFS, massDamp, stiffDamp, ts), mesh(m), indexer(nodeIndexer), material(mater), solver(linearSolver){
		memory = new double[6 * dofs];
		Initiation(memory, 6 * dofs);
		d = memory;
		v = memory + dofs;
		pre_v = memory + 2 * dofs;
		internalVirtualWork = memory + 3 * dofs;
		externalVirtualWork = memory + 4 * dofs;
		rhs = memory + 5 * dofs;

		using AssemblerType = std::conditional_t<std::is_same<SpMatrix, BlockedSymSpMatrix>::value, BlockedSymSpMatrixAssembler, SparseMatrixAssembler>;
		AssemblerType structureAssembler(dofs);
		material->getMatrixStructure(mesh, indexer, structureAssembler);

		tagentMatrix = new SpMatrix(structureAssembler);
		massMatrix = new SpMatrix(structureAssembler);
		matrixIndices = tagentMatrix->getIndices();
		material->generateMassMatrix(mesh, indexer, matrixIndices, *massMatrix);
	}

	template<class SpMatrix> NonlinearImplicitBackwardEuler<SpMatrix>::~NonlinearImplicitBackwardEuler() {
		delete[] memory;
		delete tagentMatrix;
		delete massMatrix;
	}

	template<class SpMatrix> void NonlinearImplicitBackwardEuler<SpMatrix>::runOneTimeStep() {
		memcpy(pre_v, v, sizeof(double) * dofs);

		tagentMatrix->setZeros();
		Initiation(internalVirtualWork, dofs);
		material->generateMatrixAndVirtualWorks(mesh, indexer, d, matrixIndices, *tagentMatrix, internalVirtualWork);
		tagentMatrix->Scale(timeStep + stiffnessDamping);
		tagentMatrix->Add(massDamping, *massMatrix);
		SpMDV(*tagentMatrix, pre_v, rhs);
		for (int i = 0; i < dofs; i++)
			rhs[i] = timeStep * (externalVirtualWork[i] - internalVirtualWork[i] - rhs[i]);

		tagentMatrix->Scale(timeStep);
		tagentMatrix->Add(1.0, *massMatrix);

		solver->resetLinearSystem(tagentMatrix);
		solver->solveLinearSystem(rhs, v);
		for (int i = 0; i < dofs; i++) {
			v[i] += pre_v[i];
			d[i] += v[i] * timeStep;
		}
	}

	template<class SpMatrix> void NonlinearImplicitBackwardEuler<SpMatrix>::setExternalVirtualWork(const Forcer& forcer) {
		forcer.getVirtualWorks(dofs, dofs, NULL, externalVirtualWork);
	}

	template<class SpMatrix> void NonlinearImplicitBackwardEuler<SpMatrix>::getRawDisplacements(double *displacements) const {
		memcpy(displacements, d, sizeof(double) * dofs);
	}
}

#endif