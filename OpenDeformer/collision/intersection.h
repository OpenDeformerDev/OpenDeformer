#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_COLLISION_INTERSECTION_H
#define ODER_COLLISION_INTERSECTION_H

#include "collisionShapes.h"

namespace ODER {
	namespace Collision {
		void Intersection(const StaticPlaneShape& plane, const DynamicTetrahedronShape& tet, Vector3& centroid, Vector3& dir);
		void computeBarycentricCoordinates(const DynamicTetrahedronShape& tet, const Vector3& p, Scalar coords[4]);
	}

}

#endif
