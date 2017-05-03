#include "stdafx.h"
#include "simulator.h"
#include "mesh.h"
#include "intergrator.h"
#include "mechMaterial.h"
#include "memory.h"

namespace ODER{
	Simulator::Simulator(Simulator &&simulator){
		mesh = std::move(simulator.mesh);
		indexer = std::move(simulator.indexer);
		intergrator = simulator.intergrator;

		simulator.intergrator = NULL;
	}

	Simulator& Simulator::operator=(Simulator &&simulator){
		std::swap(mesh, simulator.mesh);
		std::swap(indexer, simulator.indexer);
		std::swap(intergrator, simulator.intergrator);
		return *this;
	}

	void Simulator::getVertexPositions(Vector3f *vertices) {
		intergrator->updateMeshVerticesDisplacements(indexer, mesh);
		int vertCount = mesh->getNodeCount();
		for (int vertIndex = 0; vertIndex < vertCount; vertIndex++) {
			Vector3 vert = mesh->getVertex(vertIndex) + mesh->getVertexDisplacement(vertIndex);
			for (int axis = 0; axis < 3; axis++) {
				vertices[vertIndex][axis] = float(vert[axis]);
			}
		}
	}

	void Simulator::runTimeSteps(int stepCount){
		for (int step = 0; step < stepCount; step++)
			intergrator->runOneTimeStep();
	}

	Simulator::~Simulator(){
		delete intergrator;
	}
}