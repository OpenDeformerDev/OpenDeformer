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
		NonlinearImplicitBackwardEuler(int DOFS, Scalar massDamp, Scalar stiffDamp, Scalar ts,
			const Reference<Mesh>& m, const Reference<NodeIndexer>& nodeIndexer, const FullOrderNonlinearMaterial<SpMatrix> * mater, LinearSolver<SpMatrix>* linearSolver);
		void setExternalVirtualWork(const Forcer& forcer);
		void runOneTimeStep();
		void updateMeshVerticesDisplacements(const Reference<NodeIndexer> &indexer, Reference<Mesh> &mesh) const {}
		~NonlinearImplicitBackwardEuler();
	private:
		Scalar *v;
		Scalar *delta_v;
		Scalar *internalVirtualWork;
		Scalar *externalVirtualWork;
		Scalar *rhs;

		Scalar *memory;
		Reference<Mesh> mesh;
		Reference<NodeIndexer> indexer;
		const FullOrderNonlinearMaterial<SpMatrix> *material;
		Scalar *precomputes;
		LinearSolver<SpMatrix> *solver;
		SpMatrix *tagentMatrix;
		SpMatrix *massMatrix;
		int *matrixIndices;
	};

	template<class SpMatrix> NonlinearImplicitBackwardEuler<SpMatrix>::NonlinearImplicitBackwardEuler(int DOFS, Scalar massDamp, Scalar stiffDamp, Scalar ts,
		const Reference<Mesh>& m, const Reference<NodeIndexer>& nodeIndexer, const FullOrderNonlinearMaterial<SpMatrix> *mater, LinearSolver<SpMatrix> *linearSolver) 
		: Intergrator(DOFS, massDamp, stiffDamp, ts), mesh(m), indexer(nodeIndexer), material(mater), solver(linearSolver){
		memory = new Scalar[5 * dofs];
		Initiation(memory, 5 * dofs);
		v = memory;
		delta_v = memory + dofs;
		internalVirtualWork = memory + 2 * dofs;
		externalVirtualWork = memory + 3 * dofs;
		rhs = memory + 4 * dofs;

		using AssemblerType = SpMatrix::Assembler;
		AssemblerType structureAssembler(dofs);
		material->getMatrixStructure(mesh, indexer, structureAssembler);

		tagentMatrix = new SpMatrix(structureAssembler);
		massMatrix = new SpMatrix(structureAssembler);

		material->getMatrixIndicesPerElement(mesh, indexer, *tagentMatrix, &matrixIndices);

		solver->resetLinearSystem(tagentMatrix);
		material->generateMassMatrix(mesh, indexer, matrixIndices, *massMatrix);
		precomputes = material->getPrecomputes(mesh);
	}

	template<class SpMatrix> NonlinearImplicitBackwardEuler<SpMatrix>::~NonlinearImplicitBackwardEuler() {
		delete[] memory;
		delete[] matrixIndices;
		delete[] precomputes;
		delete tagentMatrix;
		delete massMatrix;
	}

	template<class SpMatrix> void NonlinearImplicitBackwardEuler<SpMatrix>::runOneTimeStep() {
		//memcpy(pre_v, v, sizeof(Scalar) * dofs);

		tagentMatrix->setZeros();
		Initiation(internalVirtualWork, dofs);
		Initiation(rhs, dofs);
		material->generateMatrixAndVirtualWorks(mesh, indexer, precomputes, matrixIndices, *tagentMatrix, internalVirtualWork);
		tagentMatrix->Scale(timeStep + stiffnessDamping);
		tagentMatrix->Add(massDamping, *massMatrix);
		SpMDV(*tagentMatrix, v, rhs);
		for (int i = 0; i < dofs; i++)
			rhs[i] = timeStep * (externalVirtualWork[i] - internalVirtualWork[i] - rhs[i]);

		tagentMatrix->Scale(timeStep);
		tagentMatrix->Add(1.0, *massMatrix);

		solver->solveLinearSystem(rhs, delta_v);

		for (int i = 0; i < dofs; i++) 
			v[i] += delta_v[i];

		//update mesh displacements
		auto constrainIter = indexer->getConstrainIterBegin();
		auto constrainEnd = indexer->getConstrainIterEnd();
		int dofIndex = 0;
		int vertCount = mesh->getNodeCount();
		for (int vertIndex = 0; vertIndex < vertCount; vertIndex++) {
			for (int axis = 0; axis < 3; axis++) {
				if (constrainIter != constrainEnd && (3 * vertIndex + axis) == *constrainIter) 
					constrainIter++;
				else
					mesh->getVertexDisplacement(vertIndex)[axis] += (v[dofIndex++] * timeStep);
			}
		}
	}

	template<class SpMatrix> void NonlinearImplicitBackwardEuler<SpMatrix>::setExternalVirtualWork(const Forcer& forcer) {
		forcer.getVirtualWorks(dofs, dofs, NULL, externalVirtualWork);
	}
}

#endif