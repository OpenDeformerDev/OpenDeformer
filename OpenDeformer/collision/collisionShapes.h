#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_COLLISION_COLLISIONSHAPES_H
#define ODER_COLLISION_COLLISIONSHAPES_H

#include "oder.h"
#include "latool.h"
#include "element.h"

namespace ODER {
	namespace Collision {

		struct CollisionShape {};

		struct StaticPlaneShape : public CollisionShape {
			Vector3 planeNormal;
			Scalar planeConstant;
		};

		struct DynamicTetrahedronShape : public CollisionShape {
			Vector3 points[4];
		};
	}
}

#endif