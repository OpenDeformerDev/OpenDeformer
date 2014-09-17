#pragma once
#include "mesh.h"
#include "oder.h"

namespace ODER{
	class Mesher{
	public:
		Mesher(){}
		virtual Reference<Mesh> generateMesh() = 0;
	};
}