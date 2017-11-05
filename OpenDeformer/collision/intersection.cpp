#include "stdafx.h"
#include "intersection.h"

namespace ODER {
	namespace Collision {
		namespace {
			Vector3 Intersection(const StaticPlaneShape& plane, const Vector3& pa, const Vector3& pb) {
				Vector3 ab = pb - pa;
				Scalar t = (-plane.planeConstant - plane.planeNormal * pa) / (plane.planeNormal * ab);

				return pa + Clamp(t, Scalar(0), Scalar(1)) * ab;
			}
		}

		void Intersection(const StaticPlaneShape& plane, const DynamicTetrahedronShape& tet, Vector3& centroid, Vector3& dir) {
			constexpr Scalar eps = Scalar(1e-6);

			Vector3 negPoints[4], posiPoints[4];
			int negPointsCount = 0, posiPointsCount = 0;

			for (int i = 0; i < 4; i++) {
				if (plane.planeNormal * tet.points[i] + plane.planeConstant < eps)
					negPoints[negPointsCount++] = tet.points[i];
				else
					posiPoints[posiPointsCount++] = tet.points[i];
			}

			if (negPointsCount == 0) {
				Scalar minDis = std::numeric_limits<Scalar>::max();
				for (int i = 0; i < 4; i++) {
					Scalar dis = plane.planeNormal * tet.points[i] + plane.planeConstant;
					if (dis < minDis) {
						minDis = dis;
						centroid = tet.points[i];
					}
				}
				dir = plane.planeNormal;
			}
			else if (negPointsCount == 1) {
				Vector3 p0 = Intersection(plane, negPoints[0], posiPoints[0]);
				Vector3 p1 = Intersection(plane, negPoints[0], posiPoints[1]);
				Vector3 p2 = Intersection(plane, negPoints[0], posiPoints[2]);

				centroid = (negPoints[0] + p0 + p1 + p2) * Scalar(0.25);

				Vector3 ab = p0 - negPoints[0];
				Vector3 ac = p1 - negPoints[1];
				Vector3 ad = p2 - negPoints[2];

				dir = Normalize((ab % ac) + (ac % ad) + (ad % ab));

				if (ab * (ac % ad) < Scalar(0))
					dir = -dir;
			}
			else if (negPointsCount == 2) {
				Vector3 p0 = Intersection(plane, negPoints[0], posiPoints[0]);
				Vector3 p1 = Intersection(plane, negPoints[0], posiPoints[1]);
				Vector3 p2 = Intersection(plane, negPoints[1], posiPoints[0]);
				Vector3 p3 = Intersection(plane, negPoints[1], posiPoints[1]);

				Vector3 ab = negPoints[1] - negPoints[0];

				Vector3 c0 = negPoints[0] + negPoints[1] + p2 + p3;
				Scalar v0 = fabs((negPoints[1] - negPoints[0]) * ((p2 - negPoints[0]) % (p3 - negPoints[1])));

				Vector3 c1 = negPoints[0] + p0 + p1 + p2;
				Scalar v1 = fabs((p0 - negPoints[0]) * ((p1 - negPoints[0]) % (p2 - negPoints[0])));

				Vector3 c2 = negPoints[0] + p1 + p2 + p3;
				Scalar v2 = fabs((p1 - negPoints[0]) * ((p2 - negPoints[0]) % (p3 - negPoints[0])));

				centroid = (c0 * v0 + c1 * v1 + c2 * v2) / (v0 + v1 + v2) * Scalar(0.25);

				dir = (p0 - negPoints[0]) % (p1 - negPoints[0]) + (p3 - negPoints[1]) % (p2 - negPoints[1]) +
					ab % (p2 - negPoints[0]) + (p2 - negPoints[0]) % (p0 - negPoints[0]) +
					(p3 - negPoints[0]) % ab + (p1 - negPoints[0]) % (p3 - negPoints[0]);
				dir = Normalize(dir);

				if ((negPoints[1] - negPoints[0]) *
					((posiPoints[0] - negPoints[0]) % (posiPoints[1] - negPoints[0])) < Scalar(0))
					dir = -dir;
			}
			else if (negPointsCount == 3) {
				Vector3 p0 = Intersection(plane, negPoints[0], posiPoints[0]);
				Vector3 p1 = Intersection(plane, negPoints[1], posiPoints[0]);
				Vector3 p2 = Intersection(plane, negPoints[2], posiPoints[0]);

				Vector3 c0 = negPoints[0] + p0 + p1 + p2;
				Scalar v0 = fabs((p0 - negPoints[0]) * ((p1 - negPoints[1]) % (p2 - negPoints[2])));

				Vector3 c1 = negPoints[0] + negPoints[1] + negPoints[2] + p1;
				Scalar v1 = fabs((negPoints[1] - negPoints[0]) * ((negPoints[2] - negPoints[0]) % (p1 - negPoints[0])));

				Vector3 c2 = negPoints[0] + negPoints[2] + p1 + p2;
				Scalar v2 = fabs((negPoints[2] - negPoints[0]) * ((p1 - negPoints[0]) % (p2 - negPoints[0])));

				centroid = (c0 * v0 + c1 * v1 + c2 * v2) / (v0 + v1 + v2) * Scalar(0.25);

				dir = (negPoints[2] - negPoints[0]) % (negPoints[1] - negPoints[0]) + (p0 - negPoints[0]) % (p2 - negPoints[0]) +
					(p2 - negPoints[0]) % (negPoints[2] - negPoints[0]) + (negPoints[1] - negPoints[0]) % (p2 - negPoints[0]) +
					(p2 - negPoints[0]) % (p0 - negPoints[0]) + (negPoints[2] - negPoints[1]) % (p2 - negPoints[1]) +
					(p2 - negPoints[1]) % (p1 - negPoints[1]);

				if ((negPoints[0] - posiPoints[0]) *
					((negPoints[1] - posiPoints[0]) % (negPoints[2] - posiPoints[0])) < Scalar(0))
					dir = -dir;
			}
			else {
				centroid = (tet.points[0] + tet.points[1] + tet.points[2] + tet.points[3]) * Scalar(0.25);
				dir = plane.planeNormal;
			}
		}

	}
}