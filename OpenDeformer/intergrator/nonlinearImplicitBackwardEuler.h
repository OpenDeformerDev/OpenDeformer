#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_INTERGRATOR_NONLINEARIMPLICITBACKWARDEULER_H
#define ODER_INTERGRATOR_NONLINEARIMPLICITBACKWARDEULER_H

#include "oder.h"
#include "intergrator.h"
#include "fullOrderNonlinear.h"
#include "linearSolver.h"
#include "forcer.h"

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
		double *pre_d;
		double *pre_v;
		double *internalVirtualWork;
		double *externalVirtualWork;

		double *memory;
		Reference<Mesh> mesh;
		Reference<NodeIndexer> indexer;
		const FullOrderNonlinearMaterial<SpMatrix> *material
		LinearSolver<SpMatrix> *solver;
		SpMatrix *tagentMatrix;
		SpMatrix *massMatrix;
	};

	template<class SpMatrix> NonlinearImplicitBackwardEuler<SpMatrix>::NonlinearImplicitBackwardEuler(int DOFS, double massDamp, double stiffDamp, double ts,
		const Reference<Mesh>& m, const Reference<NodeIndexer>& nodeIndexer, const FullOrderNonlinearMaterial<SpMatrix> *mater, LinearSolver<SpMatrix> *linearSolver) 
		: Intergrator(DOFS, massDamp, stiffDamp, ts), mesh(m), indexer(nodeIndexer), material(mater), solver(linearSolver){
		memory = new double[6 * dofs];
		Initiation(memory, 6 * dofs);
		d = memory;
		v = memory + dofs;
		pre_d = memory + 2 * dofs;
		pre_v = memory + 3 * dofs;
		internalVirtualWork = memory + 4 * dofs;
		externalVirtualWork = memory + 5 * dofs;

		using AssemblerType = std::conditional_t<std::is_same<SpMatrix, BlockedSymSpMatrix>::value, BlockedSymSpMatrixAssembler, SparseMatrixAssembler>;
		AssemblerType structureAssembler, massAssembler;
		material->generateMassMatrix(mesh, indexer, massAssembler);
		material->getMatrixStructure(mesh, indexer, structureAssembler);

		tagentMatrix = new SpMatrix(structureAssembler);
		massMatrix = new SpMatrix(massAssembler);
		linearSolver->resetLinerSystem(tagentMatrix);
	}

	template<class SpMatrix> NonlinearImplicitBackwardEuler<SpMatrix>::~NonlinearImplicitBackwardEuler() {
		delete[] memory;
		delete tagentMatrix;
		delete massMatrix;
	}

	template<class SpMatrix> void NonlinearImplicitBackwardEuler<SpMatrix>::setExternalVirtualWork(const Forcer& forcer) {
		forcer.getVirtualWorks(dofs, dofs, NULL, externalVirtualWork);
	}

	template<class SpMatrix> void NonlinearImplicitBackwardEuler<SpMatrix>::getRawDisplacements(double *displacements) const {
		memcpy(displacements, d, sizeof(double) * dofs);
	}

}

#endif