#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_CORE_MESHER_H
#define ODER_CORE_MESHER_H

#include "mesh.h"
#include "oder.h"

namespace ODER{
	class Mesher{
	public:
		virtual Reference<Mesh> generateMesh(int *vertexLableMap = NULL) = 0;
		virtual ~Mesher() = default;
	};
}

#endif