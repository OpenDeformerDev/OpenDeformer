#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_CORE_SIMULATOR_H
#define ODER_CORE_SIMULATOR_H

#include "oder.h"
#include "memory.h"

namespace ODER{
	class Simulator{
	public:
		Simulator():intergrator(NULL), rawDisplacements(NULL){}
		Simulator(const Reference<Mesh> m, const Reference<NodeIndexer>& nodeIndexer, Intergrator *intergrate);
		Simulator(Simulator &&simulator);
		Simulator& operator=(Simulator &&simulator);
		Simulator(const Simulator &simulator) = delete;
		Simulator& operator=(const Simulator &simulator) = delete;
		void getVertexPositions(Vector *vertices) const;
		void runTimeSteps(int stepCount);
		~Simulator();
	private:
		Reference<Mesh> mesh;
		Reference<NodeIndexer> indexer;
		Intergrator *intergrator;
		double *rawDisplacements;
	};
}

#endif