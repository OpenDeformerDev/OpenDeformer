#include "stdafx.h"
#include "intersection.h"
#include "delPrimitive.h"
#include "predicate.h"
#include "geometer.h"
#include <queue>
#include <random>

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

		void digCavity(Vertex *u, const Segment& s, const Predicator<Scalar>& predicator, 
			TriMeshDataStructure& mesh, std::vector<int>& deletedIndices, 
			std::vector<std::pair<Scalar, std::vector<Vertex *>>>& conflictVertsVecs,
			std::vector<Vertex *> conflictVerts,
			std::vector<SegmentWithIndex>& horizon) {
			Vertex *w = s.v[0], *v = s.v[1], *x = NULL;
			int index = -1;
			if (!mesh.Adjacent(s, &x, &index)) return;

			if (predicator.orient3d(u->point, w->point, v->point, x->point) > Scalar(0)) {
				mesh.deleteTriangle(w, v, x);
				if (index > 0) {
					deletedIndices.push_back(index);
					for (auto vert : conflictVertsVecs[index].second) conflictVerts.push_back(vert);

					conflictVertsVecs[index].first = Scalar(0);
					conflictVertsVecs[index].second.clear();
				}

				digCavity(u, Segment(x, v), predicator, mesh, deletedIndices, conflictVertsVecs, conflictVerts, horizon);
				digCavity(u, Segment(w, x), predicator, mesh, deletedIndices, conflictVertsVecs, conflictVerts, horizon);
			}
			else {
				horizon.push_back(SegmentWithIndex(v, w, -1));
				mesh.addTriangle(u, v, w);
			}
		}

		void computeConvexHull(std::vector<Vertex *>& verts, TriMeshDataStructure& mesh) {
			if (verts.size() < 4) return;

			Predicator<Scalar> predicator;

			Vertex *p0 = verts[0];
			Vertex *p1 = verts[1];
			Vertex *p2 = NULL;

			for (int i = 2; i < verts.size(); i++) {
				if (!predicator.coLine(p0->point, p1->point, verts[i]->point)) {
					p2 = verts[2];
					std::swap(verts[2], verts[i]);
					break;
				}
			}

			Vertex *p3 = NULL;
			if (p2) {
				Scalar maxOrient = Scalar(0);
				int maxIndex = -1;
				for (int i = 3; i < verts.size(); i++) {
					Scalar orient = predicator.orient3d(p0->point, p1->point, p2->point, verts[i]->point);
					if (fabs(orient) > fabs(maxOrient)) {
						maxIndex = i;
						maxOrient = orient;
					}
				}
				if (maxIndex >= 0) {
					p3 = verts[maxIndex];
					std::swap(verts[3], verts[maxIndex]);

					if (maxOrient < Scalar(0)) std::swap(p0, p1);
				}
			}

			if (p3) {
				std::vector<std::pair<Scalar, std::vector<Vertex *>>> conflictVertsVecs;
				std::vector<int> reusableIndices;
				std::vector<SegmentWithIndex> horizon; 
				std::vector<Vertex *> conflictVerts;
				conflictVertsVecs.reserve(verts.size());
				reusableIndices.reserve(verts.size());
				horizon.reserve(verts.size());
				conflictVerts.reserve(verts.size());


				TriangleWithIndex facets[4];
				facets[0] = TriangleWithIndex(p0, p1, p2, -1);
				facets[1] = TriangleWithIndex(p0, p2, p3, -1);
				facets[2] = TriangleWithIndex(p1, p3, p2, -1);
				facets[3] = TriangleWithIndex(p0, p3, p1, -1);

				for (int i = 4; i < verts.size(); i++) {
					Vertex *v = verts[i];
					for (int j = 0; j < 4; j++) {
						Scalar orient = predicator.orient3d(v->point, facets[j].v[0]->point, facets[j].v[1]->point, facets[j].v[2]->point);
						if (orient > Scalar(0)) {
							if (facets[j].index < 0) {
								facets[j].index = conflictVertsVecs.size();
								std::vector<Vertex *> conflictVerts;
								conflictVerts.push_back(v);
								conflictVertsVecs.push_back(std::make_pair(orient, std::move(conflictVerts)));
							}
							else {
								conflictVertsVecs[facets[j].index].second.push_back(v);
								if (orient > conflictVertsVecs[facets[j].index].first) {
									conflictVertsVecs[facets[j].index].first = orient;
									std::swap(conflictVertsVecs[facets[j].index].second.back(), conflictVertsVecs[facets[j].index].second.front());
								}
							}
							break;
						}
					}
				}
				for (int i = 0; i < 4; i++) mesh.addTriangle(facets[i].v[0], facets[i].v[1], facets[i].v[2], facets[i].index);
				
				std::queue<TriangleWithIndex> conflictedTriangles;
				for (int i = 0; i < 4; i++)
					if (facets[i].index >= 0)
						conflictedTriangles.push(facets[i]);

				while (!conflictedTriangles.empty()) {
					TriangleWithIndex f = conflictedTriangles.front();
					conflictedTriangles.pop();

					Vertex *conflict = conflictVertsVecs[f.index].second.front();

					for (int i = 1; i < conflictVertsVecs[f.index].second.size(); i++)
						conflictVerts.push_back(conflictVertsVecs[f.index].second[i]);
					conflictVertsVecs[f.index].first = Scalar(0);
					conflictVertsVecs[f.index].second.clear();
					reusableIndices.push_back(f.index);

					mesh.deleteTriangle(f.v[0], f.v[1], f.v[2]);	

					digCavity(conflict, Segment(f.v[1], f.v[0]), predicator, mesh, reusableIndices, conflictVertsVecs, conflictVerts, horizon);
					digCavity(conflict, Segment(f.v[0], f.v[2]), predicator, mesh, reusableIndices, conflictVertsVecs, conflictVerts, horizon);
					digCavity(conflict, Segment(f.v[2], f.v[1]), predicator, mesh, reusableIndices, conflictVertsVecs, conflictVerts, horizon);

					for (auto vert : conflictVerts) {
						for (auto& seg : horizon) {
							Scalar orient = predicator.orient3d(vert->point, conflict->point, seg.v[0]->point, seg.v[1]->point);
							if (orient > Scalar(0)) {
								if (seg.index < 0) {
									if (!reusableIndices.empty()) {
										int newIndex = reusableIndices.back();
										reusableIndices.pop_back();

										conflictVertsVecs[newIndex].first = orient;
										conflictVertsVecs[newIndex].second.push_back(vert);
										seg.index = newIndex;
									}
									else {
										int newIndex = reusableIndices.size();
										std::vector<Vertex *> conflictVertsVec;
										conflictVertsVec.push_back(vert);
										conflictVertsVecs.push_back(std::make_pair(orient, std::move(conflictVertsVec)));
										seg.index = newIndex;
									}
								}
								else {
									conflictVertsVecs[seg.index].second.push_back(vert);
									if (orient > conflictVertsVecs[seg.index].first) {
										conflictVertsVecs[seg.index].first = orient;
										std::swap(conflictVertsVecs[seg.index].second.front(), conflictVertsVecs[seg.index].second.back());
									}
								}
								break;
							}
						}
					}

					for (auto seg : horizon)
						if (seg.index >= 0)
							conflictedTriangles.push(TriangleWithIndex(conflict, seg.v[0], seg.v[1], seg.index));

					horizon.clear();
					conflictVerts.clear();
				}
			}

		}

		template<int d> bool solveLinearProgramming(int n, Scalar *c, Scalar *A, Scalar *u, Scalar *v) {
			std::random_device rng;
			std::default_random_engine randomEngine(rng());
			for (int i = n - 1; i > 0; i--) {
				std::uniform_int_distribution<int> indexDistrubution(0, i);
				int ni = indexDistrubution(randomEngine);
				for (int j = 0; j < d + 2; j++)
					std::swap(A[i * (d + 2) + j], A[ni * (d + 2) + j]);
			}

			Scalar uu[d], vv[d];
			for (int i = 0; i < d; i++) {
				uu[i] = Scalar(0);
				vv[i] = (c[i] >= Scalar(0) ? Scalar(1) : Scalar(-1));
			}

			Scalar *AA = new Scalar[(d + 1) * (n + 1)];

			for (int i = 0; i < n; i++) {
				Scalar a[d + 2];
				for (int j = 0; j < d + 2; j++) a[j] = A[i * (d + 2) + j];

				Scalar auu = Scalar(0), avv = Scalar(0);
				for (int j = 0; j < d; j++) {
					auu += a[j] * uu[j];
					avv += a[j] * vv[j];
				}

				if ((avv > a[d + 1]) || (avv == a[d + 1] && auu > a[d])) {
					Scalar ak = a[0];
					int k = 0;
					for (int j = 1; j < d; j++) {
						if (fabs(a[j]) > fabs(ak)) {
							ak = a[j];
							k = j;
						}
					}

					if (ak == Scalar(0)) {
						delete[] AA;
						return false;
					}

					//eliminate kth variable from constraints
					for (int j = 0; j < i; j++) {
						Scalar scale = A[j * (d + 2) + k] / ak;
						int ind = 0;
						Scalar aaaa = A[j * (d + 2) + k] - scale * a[k];
						for (int l = 0; l < d + 2; l++) {
							if (l != k) {
								AA[j * (d + 1) + ind] = A[j * (d + 2) + l] - scale * a[l];
								ind += 1;
							}
						}
					}
					Scalar cc[d - 1];
					{
						Scalar scale = c[k] / ak;
						int ind = 0;
						for (int j = 0; j < d; j++) {
							if (j != k) {
								cc[ind] = c[j] - scale * a[j];
								ind += 1;
							}
						}
					}
					//incorporate the constraint on kth variable
					{
						int ind = 0;
						Scalar scale = Scalar(1) / ak;
						for (int j = 0; j < d + 1; j++) {
							if (j != k) {
								AA[i * (d + 1) + ind] = -scale * a[j];
								AA[(i + 1) * (d + 1) + ind] = scale * a[j];
								ind += 1;
							}
						}
						AA[i * (d + 1) + d] = Scalar(1) - scale * a[d + 1];
						AA[(i + 1) * (d + 1) + d] = Scalar(1) + scale * a[d + 1];
					}

 					if (!solveLinearProgramming<d - 1>(i + 2, cc, AA, uu, vv)) {
						delete[] AA;
						return false;
					}

					//"lift" the solution
					for (int j = d - 2; j >= k; j--) {
						uu[j + 1] = uu[j];
						vv[j + 1] = vv[j];
					}
					uu[k] = Scalar(0);
					vv[k] = Scalar(0);

					auu = Scalar(0); avv = Scalar(0);
					for (int j = 0; j < d; j++) {
						auu += a[j] * uu[j];
						avv += a[j] * vv[j];
					}
					uu[k] = (a[d] - auu) / ak;
					vv[k] = (a[d + 1] - avv) / ak;
				}
			}

			for (int i = 0; i < d; i++) {
				u[i] = uu[i];
				v[i] = vv[i];
			}

			delete[] AA;
			return true;
		}

		template<> bool solveLinearProgramming<1>(int n, Scalar *c, Scalar *A, Scalar *u, Scalar *v) {
			Scalar h = Scalar(0), h_lambda = Scalar(1);
			Scalar l = Scalar(0), l_lambda = Scalar(-1);
			Scalar z = Scalar(0), z_lambda = Scalar(0);

			for (int i = 0; i < n; i++) {
				Scalar a1 = A[i * 3];
				Scalar a2 = A[i * 3 + 1];
				Scalar a3 = A[i * 3 + 2];
				if (a1 != Scalar(0)) {
					a2 /= a1;
					a3 /= a1;

					if (a1 > Scalar(0)) {
						if (a3 < h_lambda) {
							h_lambda = a3;
							h = a2;
						}
						else if (a3 == h_lambda)
							h = std::min(h, a2);
					}
					else {
						if (a3 > l_lambda) {
							l_lambda = a3;
							l = a2;
						}
						else if (a3 == l_lambda)
							l = std::max(l, a2);
					}
				}
				else {
					if (a3 < z_lambda) {
						z_lambda = a3;
						z = a2;
					}
					else if (a3 == z_lambda)
						z = std::min(z, a2);
				}
			}

			if (z_lambda < Scalar(0) || (z_lambda == Scalar(0) && z < Scalar(0)) || 
				l_lambda > h_lambda || (l_lambda == h_lambda && l > h)) 
				return false;

			if (*c >= Scalar(0)) {
				*u = h; *v = h_lambda;
			}
			else {
				*u = l; *v = l_lambda;
			}
			return true;
		}

		void Intersection(const DynamicTetrahedronShape& teta, const DynamicTetrahedronShape& tetb, Vector3& centroid, Vector3& dir) {
			Vector3 tetaPoints[4] = { teta.points[0], teta.points[1], teta.points[2], teta.points[3] };
			Vector3 tetbPoints[4] = { tetb.points[0], tetb.points[1], tetb.points[2], tetb.points[3] };

			Predicator<Scalar> predicator;
			if (predicator.orient3d(tetaPoints[0], tetaPoints[1], tetaPoints[2], tetaPoints[3]) < Scalar(0)) std::swap(tetaPoints[2], tetaPoints[3]);
			if (predicator.orient3d(tetbPoints[0], tetbPoints[1], tetbPoints[2], tetbPoints[3]) < Scalar(0)) std::swap(tetbPoints[2], tetbPoints[3]);

			Vector3 normals[8];
			Scalar constants[8];

			normals[0] = Normalize(Geometer::triangleNormal(tetaPoints[1], tetaPoints[2], tetaPoints[3]));
			normals[1] = Normalize(Geometer::triangleNormal(tetaPoints[0], tetaPoints[3], tetaPoints[2]));
			normals[2] = Normalize(Geometer::triangleNormal(tetaPoints[0], tetaPoints[1], tetaPoints[3]));
			normals[3] = Normalize(Geometer::triangleNormal(tetaPoints[0], tetaPoints[2], tetaPoints[1]));

			normals[4] = Normalize(Geometer::triangleNormal(tetbPoints[1], tetbPoints[2], tetbPoints[3]));
			normals[5] = Normalize(Geometer::triangleNormal(tetbPoints[0], tetbPoints[3], tetbPoints[2]));
			normals[6] = Normalize(Geometer::triangleNormal(tetbPoints[0], tetbPoints[1], tetbPoints[3]));
			normals[7] = Normalize(Geometer::triangleNormal(tetbPoints[0], tetbPoints[2], tetbPoints[1]));

			constants[0] = -(normals[0] * tetaPoints[1]);
			constants[1] = -(normals[1] * tetaPoints[0]);
			constants[2] = -(normals[2] * tetaPoints[0]);
			constants[3] = -(normals[3] * tetaPoints[0]);

			constants[4] = -(normals[4] * tetbPoints[1]);
			constants[5] = -(normals[5] * tetbPoints[0]);
			constants[6] = -(normals[6] * tetbPoints[0]);
			constants[7] = -(normals[7] * tetbPoints[0]);

			Scalar A[8 * 6];
			for (int i = 0; i < 8; i++) {
				A[i * 6 + 0] = -normals[i].x;
				A[i * 6 + 1] = -normals[i].y;
				A[i * 6 + 2] = -normals[i].z;
				A[i * 6 + 3] = Scalar(1);
				A[i * 6 + 4] = constants[i];
				A[i * 6 + 5] = Scalar(0);
			}
			Scalar c[4] = { 0, 0, 0, 1 };
			Scalar u[4], v[4];

			bool success = solveLinearProgramming<4>(8, c, A, u, v);

			if (!success || u[3] <= Scalar(0)) {
				dir = Vector3(0, 0, 0);
				return;
			}

			std::vector<Vertex *> vertices;
			vertices.reserve(8);
			TriMeshDataStructure mesh;

			constexpr Scalar eps = 1e-8;
			Vector3 origin(u[0], u[1], u[2]);
			for (int i = 0; i < 8; i++) {
				Vector3 dualPoint = normals[i] / (constants[i] + normals[i] * origin);

				bool merged = false;
				for (auto vert : vertices) {
					if ((vert->point - dualPoint).length2() <= eps) {
						merged = true;
						break;
					}
				}
				if (!merged) vertices.push_back(mesh.allocVertex(dualPoint, 0));
			}

			computeConvexHull(vertices, mesh);

			std::vector<std::vector<Vector3>> polyhedron;
			for (auto vert : vertices) {
				auto begin = mesh.getIncidentTriangles(vert);
				std::vector<Vector3> polygon;
				polygon.reserve(8);

				if (begin != NULL) {
					auto circulator = begin;

					Triangle triangle = *circulator;
					Vector3 normal = Geometer::triangleNormal(triangle.v[0]->point, triangle.v[1]->point, triangle.v[2]->point);
					Scalar constant = -(normal * triangle.v[0]->point);
					Vector3 primalPoint = normal / -constant + origin;
					polygon.push_back(primalPoint);

					while (++circulator != begin) {
						triangle = *circulator;
						normal = Geometer::triangleNormal(triangle.v[0]->point, triangle.v[1]->point, triangle.v[2]->point);
						constant = -(normal * vert->point);
						Vector3 nextPrimalPoint = normal / -constant + origin;

						Scalar e = (nextPrimalPoint - polygon.back()).length2();

						if ((nextPrimalPoint - polygon.back()).length2() > eps)
							polygon.push_back(nextPrimalPoint);
					}
					
					if (polygon.size() >= 3) polyhedron.push_back(std::move(polygon));
				}
			}
			
			Vector3 cent = Vector3(0, 0, 0), direction = Vector3(0, 0, 0);
			Scalar volume = 0;
			Vector3 oppoPoint = polyhedron[0][0];

			for (std::vector<std::vector<Vector3>>::size_type i = 1; i < polyhedron.size(); i++) {
				Vector3 orgPoint = polyhedron[i][0];

				for (std::vector<Vector3>::size_type j = 1; j < polyhedron[i].size() - 1; j++) {
					Vector3 destPoint = polyhedron[i][j];
					Vector3 apexPoint = polyhedron[i][j + 1];


					Scalar tetVolume = fabs(Geometer::triangleNormal(orgPoint, destPoint, apexPoint) * (orgPoint - oppoPoint));
					cent += tetVolume * Scalar(0.25) * (oppoPoint + orgPoint + destPoint + apexPoint);
					volume += tetVolume;
				}
			}
			centroid = cent / volume;

			for (auto &polygon : polyhedron) {
				Vector3 orgPoint = polygon[0];
				for (std::vector<Vector3>::size_type i = 1; i < polygon.size() - 1; i++) {
					Vector3 destPoint = polygon[i];
					Vector3 apexPoint = polygon[i + 1];

					Vector3 normal = Geometer::triangleNormal(orgPoint, destPoint, apexPoint);
					direction += normal;
				}
			}

			dir = volume * Normalize(direction);
		}
	}
}