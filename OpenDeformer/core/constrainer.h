#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_CORE_CONSTRAINER_H
#define ODER_CORE_CONSTRAINER_H

#include "oder.h"
#include "memory.h"
#include "latool.h"
#include "sparseMatrix.h"
#include "mesh.h"
#include "forcer.h"
#include "element.h"
#include <set>
#include <vector>

namespace ODER{
	enum ConstrainType{
		Constrain_X,
		Constrain_Y,
		Constrain_Z,
		Constrain_XY,
		Constrain_YZ,
		Constrain_XZ,
		Constrain_XYZ
	};

	class Constrainer{
	public:
		Constrainer(const Reference<Mesh> &m) :mesh(m){}
		void setConstrains(int nodeCounts, int *nodeIndices, ConstrainType *constrainTypes = NULL);
		void setConstrains(const Vector& origin, const Vector& normal, ConstrainType constrainType = Constrain_XYZ);
		void outputConstrains(std::vector<int> &out) const;
	private:
		void setConstrain(int nodeIndex, ConstrainType constrainType);
		std::set<int> constrains;
		Reference<Mesh> mesh;
	};
}

#endif