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
		Simulator():intergrator(NULL){}
		Simulator(const Reference<Mesh> m, const Reference<NodeIndexer>& nodeIndexer, Intergrator *intergrate)
			: mesh(m), indexer(nodeIndexer), intergrator(intergrate) {}
		Simulator(Simulator &&simulator);
		Simulator& operator=(Simulator &&simulator);
		Simulator(const Simulator &simulator) = delete;
		Simulator& operator=(const Simulator &simulator) = delete;
		void getVertexPositions(Vector3f *vertices);
		void runTimeSteps(int stepCount);
		~Simulator();
	private:
		Reference<Mesh> mesh;
		Reference<NodeIndexer> indexer;
		Intergrator *intergrator;
	};
}

#endif