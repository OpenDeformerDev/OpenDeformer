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
				dir = Vector3(0, 0, 0);
			}
			else if (negPointsCount == 1) {
				Vector3 p0 = Intersection(plane, negPoints[0], posiPoints[0]);
				Vector3 p1 = Intersection(plane, negPoints[0], posiPoints[1]);
				Vector3 p2 = Intersection(plane, negPoints[0], posiPoints[2]);

				centroid = (negPoints[0] + p0 + p1 + p2) * Scalar(0.25);

				Vector3 ab = p0 - negPoints[0];
				Vector3 ac = p1 - negPoints[0];
				Vector3 ad = p2 - negPoints[0];

				Scalar svol = ab * (ac % ad);
				dir = fabs(svol) * Normalize((ab % ac) + (ac % ad) + (ad % ab));

				if (svol < Scalar(0))
					dir = -dir;
			}
			else if (negPointsCount == 2) {
				Vector3 p0 = Intersection(plane, negPoints[0], posiPoints[0]);
				Vector3 p1 = Intersection(plane, negPoints[0], posiPoints[1]);
				Vector3 p2 = Intersection(plane, negPoints[1], posiPoints[0]);
				Vector3 p3 = Intersection(plane, negPoints[1], posiPoints[1]);

				Vector3 ab = negPoints[1] - negPoints[0];

				Vector3 c0 = negPoints[0] + negPoints[1] + p2 + p3;
				Scalar v0 = fabs((negPoints[1] - negPoints[0]) * ((p2 - negPoints[0]) % (p3 - negPoints[0])));

				Vector3 c1 = negPoints[0] + p0 + p1 + p2;
				Scalar v1 = fabs((p0 - negPoints[0]) * ((p1 - negPoints[0]) % (p2 - negPoints[0])));

				Vector3 c2 = negPoints[0] + p1 + p2 + p3;
				Scalar v2 = fabs((p1 - negPoints[0]) * ((p2 - negPoints[0]) % (p3 - negPoints[0])));

				centroid = (c0 * v0 + c1 * v1 + c2 * v2) / (v0 + v1 + v2) * Scalar(0.25);

				dir = (p0 - negPoints[0]) % (p1 - negPoints[0]) + (p3 - negPoints[1]) % (p2 - negPoints[1]) +
					ab % (p2 - negPoints[0]) + (p2 - negPoints[0]) % (p0 - negPoints[0]) +
					(p3 - negPoints[0]) % ab + (p1 - negPoints[0]) % (p3 - negPoints[0]);
				dir = (v0 + v1 + v2) * Normalize(dir);

				if ((negPoints[1] - negPoints[0]) *
					((posiPoints[0] - negPoints[0]) % (posiPoints[1] - negPoints[0])) < Scalar(0))
					dir = -dir;
			}
			else if (negPointsCount == 3) {
				Vector3 p0 = Intersection(plane, negPoints[0], posiPoints[0]);
				Vector3 p1 = Intersection(plane, negPoints[1], posiPoints[0]);
				Vector3 p2 = Intersection(plane, negPoints[2], posiPoints[0]);

				Vector3 ab = p0 - posiPoints[0];
				Vector3 ac = p1 - posiPoints[0];
				Vector3 ad = p2 - posiPoints[0];

				Vector3 c = tet.points[0] + tet.points[1] + tet.points[2] + tet.points[3];
				Scalar sv = (negPoints[0] - posiPoints[0]) * ((negPoints[1] - posiPoints[0]) % (negPoints[2] - posiPoints[0]));
				Scalar v = fabs(sv);

				Vector3 c0 = posiPoints[0] + p0 + p1 + p2;
				Scalar v0 = fabs(ab * (ac % ad));

				centroid = (v * c - v0 * c0) / (v - v0) * Scalar(0.25);
				
				dir = (v - v0) * Normalize((ab % ac) + (ac % ad) + (ad % ab));

				if (sv > Scalar(0))
					dir = -dir;
			}
			else {
				centroid = (tet.points[0] + tet.points[1] + tet.points[2] + tet.points[3]) * Scalar(0.25);
				Scalar vol = fabs((tet.points[1] - tet.points[0]) * ((tet.points[2] - tet.points[0]) % (tet.points[3] - tet.points[0])));
				dir = vol * plane.planeNormal;
			}
		}

		void computeBarycentricCoordinates(const DynamicTetrahedronShape& tet, const Vector3& p, Scalar coords[4]) {
			constexpr Scalar eps = Scalar(1e-6);

			Vector3 ab = tet.points[1] - tet.points[0];
			Vector3 ac = tet.points[2] - tet.points[0];
			Vector3 ad = tet.points[3] - tet.points[0];
			Vector3 bc = tet.points[2] - tet.points[1];
			Vector3 bd = tet.points[3] - tet.points[1];

			Vector3 ap = p - tet.points[0];
			Vector3 bp = p - tet.points[1];

			Scalar va = fabs(bp * (bc % bd));
			Scalar vb = fabs(ap * (ac % ad));
			Scalar vc = fabs(ap * (ab % ad));
			Scalar vd = fabs(ap * (ab % ac));

			Scalar vol = fabs(ab * (ac % ad));

			if (vol <= eps) {
				coords[0] = coords[1] = coords[2] = coords[3] = Scalar(0.25);
				return;
			}

			Scalar invVolume = 1 / vol;

			coords[0] = va * invVolume;
			coords[1] = vb * invVolume;
			coords[2] = vc * invVolume;
			coords[3] = vd * invVolume;
		}

	}
}