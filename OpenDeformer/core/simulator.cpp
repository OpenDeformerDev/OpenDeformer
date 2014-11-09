#include "stdafx.h"
#include "simulator.h"
#include "mesh.h"
#include "intergrator.h"
#include "mecMaterial.h"
#include "memory.h"

namespace ODER{
	Simulator::Simulator(const Reference<Mesh> m, const Reference<NodeIndexer>& nodeIndexer, Intergrator *intergrate)
		:mesh(m), indexer(nodeIndexer), intergrator(intergrate){
		rawDisplacements = allocAligned<double>(indexer->getMatrixOrder(mesh));
	};
	void Simulator::getVertexPositions(Vector *vertices) const{
		const int* constrainIndices = NULL;
		int constrainSize = indexer->getConstrainIndices(&constrainIndices);
		intergrator->getRawDisplacements(rawDisplacements);

		int constrainIndex = 0;
		int displacementIndex = 0;
		int vertCount = mesh->getNodeCount();
		for (int vertIndex = 0; vertIndex < vertCount; vertIndex++){
			for (int axis = 0; axis < 3; axis++){
				if (constrainIndex < constrainSize && (3 * vertIndex + axis) == constrainIndices[constrainIndex]){
					vertices[vertIndex][axis] = mesh->getVertex(vertIndex)[axis];
					constrainIndex++;
				}
				else
					vertices[vertIndex][axis] = mesh->getVertex(vertIndex)[axis] + float(rawDisplacements[displacementIndex++]);
			}
		}
	}
	void Simulator::runTimeSteps(int stepCount){
		for (int step = 0; step < stepCount; step++)
			intergrator->runOneTimeStep();
	}
	Simulator::~Simulator(){
		delete intergrator;
		freeAligned(rawDisplacements);
	}
}