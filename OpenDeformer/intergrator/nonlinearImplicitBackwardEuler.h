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
			const Reference<Mesh>& m, const Reference<NodeIndexer>& nodeIndexer, 
			const Reference<FullOrderNonlinearMaterial<SpMatrix>>& mater, const Reference<LinearSolver<SpMatrix>>& linearSolver);
		void setExternalVirtualWork(const Forcer& forcer);
		void addExternalDampingForceMatrix(const DampingForcer& forcer);
		void clearExternalDampingForceMatrix();
		void runOneTimeStep();
		void updateMeshVerticesDisplacements(const Reference<NodeIndexer> &indexer, Reference<Mesh> &mesh) const;
		void getMeshVerticesVelocities(const Reference<NodeIndexer> &indexer, const Reference<Mesh> &mesh, Vector3 *velocities) const;
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
		const Reference<FullOrderNonlinearMaterial<SpMatrix>> material;
		FullOrderNonlinearMaterialCache materialCache;
		Reference<LinearSolver<SpMatrix>> solver;
		SpMatrix *tagentMatrix;
		SpMatrix *massMatrix;
		SpMatrix *dampMatrix;
		SparseSymMatrixIndicesPerElementCache matrixIndices;
	};

	template<class SpMatrix> NonlinearImplicitBackwardEuler<SpMatrix>::NonlinearImplicitBackwardEuler(int DOFS, Scalar massDamp, Scalar stiffDamp, Scalar ts,
		const Reference<Mesh>& m, const Reference<NodeIndexer>& nodeIndexer, 
		const Reference<FullOrderNonlinearMaterial<SpMatrix>>& mater, const Reference<LinearSolver<SpMatrix>>& linearSolver)
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
		dampMatrix = new SpMatrix(structureAssembler);

		matrixIndices = material->getMatrixIndicesPerElement(mesh, indexer, *tagentMatrix);

		solver->resetLinearSystem(tagentMatrix);
		material->generateMassMatrix(mesh, indexer, matrixIndices, *massMatrix);
		materialCache = material->getPrecomputes(mesh);
	}

	template<class SpMatrix> NonlinearImplicitBackwardEuler<SpMatrix>::~NonlinearImplicitBackwardEuler() {
		delete[] memory;
		delete tagentMatrix;
		delete massMatrix;
		delete dampMatrix;
	}

	template<class SpMatrix> void NonlinearImplicitBackwardEuler<SpMatrix>::runOneTimeStep() {
		tagentMatrix->setZeros();
		Initiation(internalVirtualWork, dofs);
		Initiation(rhs, dofs);
		material->generateMatrixAndVirtualWorks(mesh, indexer, materialCache, matrixIndices, *tagentMatrix, internalVirtualWork);
		tagentMatrix->Scale(timeStep + stiffnessDamping);
		tagentMatrix->Add(massDamping, *massMatrix);
		tagentMatrix->Add(Scalar(1), *dampMatrix);
		SpMDV(*tagentMatrix, v, rhs);
		for (int i = 0; i < dofs; i++)
			rhs[i] = timeStep * (externalVirtualWork[i] - internalVirtualWork[i] - rhs[i]);

		tagentMatrix->Scale(timeStep);
		tagentMatrix->Add(Scalar(1), *massMatrix);

		solver->solveLinearSystem(rhs, delta_v);

		for (int i = 0; i < dofs; i++) 
			v[i] += delta_v[i];
	}

	template<class SpMatrix> void NonlinearImplicitBackwardEuler<SpMatrix>::updateMeshVerticesDisplacements(const Reference<NodeIndexer> &indexer, Reference<Mesh> &mesh) const {
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

	template<class SpMatrix> void NonlinearImplicitBackwardEuler<SpMatrix>::addExternalDampingForceMatrix(const DampingForcer& forcer) {
		forcer.addDampingMatrix(matrixIndices, *dampMatrix);
	}

	template<class SpMatrix> void NonlinearImplicitBackwardEuler<SpMatrix>::clearExternalDampingForceMatrix() {
		dampMatrix->setZeros();
	}

	template<class SpMatrix> void NonlinearImplicitBackwardEuler<SpMatrix>::getMeshVerticesVelocities(const Reference<NodeIndexer> &indexer, 
		const Reference<Mesh> &mesh, Vector3 *velocities) const {
		auto constrainIter = indexer->getConstrainIterBegin();
		auto constrainEnd = indexer->getConstrainIterEnd();
		int dofIndex = 0;
		int vertCount = mesh->getNodeCount();
		for (int vertIndex = 0; vertIndex < vertCount; vertIndex++) {
			for (int axis = 0; axis < 3; axis++) {
				if (constrainIter != constrainEnd && (3 * vertIndex + axis) == *constrainIter) {
					velocities[vertIndex][axis] = Scalar(0);
					constrainIter++;
				}
				else
					velocities[vertIndex][axis] = v[dofIndex++];
			}
		}
	}
}

#endif