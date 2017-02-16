#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_MESHER_DELMESHER_H
#define ODER_MESHER_DELMESHER_H

#include "predicate.h"
#include "delPrimitive.h"
#include "mesher.h"
#include "aabb.h"
#include <vector>
#include <deque>

namespace ODER{
	class DelTriangulator{
	public:
		DelTriangulator() = default;
		~DelTriangulator() = default;
		void generateSubPolygons(Vertex **vertices, int *segments, int vertexCount, int segmentCount, 
			const DelVector& above, bool boundaryOnly, std::vector<Triangle>& output);
		void generateSubPolygons(Vertex **vertices, Segment *segments, int vertexCount, int segmentCount, 
			const DelVector& above, bool boundaryOnly, std::vector<Triangle>& output);
		void insertSegments(const Triangle *triangles, const Segment *segments, int triangleCount, int segmentCount, 
			const DelVector& above, std::vector<Triangle>& output);
	private:
		void insertSegment(const Segment& s);
		void triangulateHalfHole(const std::vector<Vertex *>& vertices);
		void findCavity(const Segment& s, std::vector<Vertex *>& positive, std::vector<Vertex *>& negtive);
		void insertVertexToCavity(Vertex *u, Vertex *v, Vertex *w,
			bool mark, bool oriTest, std::deque<Vertex *>& marked, int depth);
		void insertVertexToCavity(Vertex *u, Vertex *v, Vertex *w, int depth);
		void triangulateConvexPoly(Vertex *u, const std::vector<Vertex *>& convexPoly);
		Triangle findPosition(Vertex *u, const Triangle &f) const;
		void digCavity(Vertex *u, const Segment &s, Triangle *rf, int depth);
		void propagateClean(const Segment& s, int depth);
		inline int findGhost(const Triangle &f){
			for (int i = 0; i < 3; i++){
				if (f.v[i]->isGhost())
					return i;
			}
			return -1;
		}
		inline int findGhost(Vertex **v, int n){
			for (int i = 0; i < n; i++){
				if (v[i]->isGhost())
					return i;
			}
			return -1;
		}

		static Predicator<REAL> predicator;
		TriMeshDataStructure meshRep;
		Vertex *ghost;

		std::vector<Segment> segments;

		TriMeshDataStructure cavityRep;
	};

	class DelMesher : public Mesher{
	public:
		struct DelMesherCriteria {
			DelMesherCriteria() {
				constexpr REAL maxVal = std::numeric_limits<REAL>::max();
				maxTetRadius = maxVal; maxFacetRadius = maxVal; maxSegLength = maxVal;
				maxTetRatio = REAL(2.0);
				facetAngleTol = REAL(179.0);
			}
			REAL maxTetRadius;
			REAL maxFacetRadius;
			REAL maxSegLength;
			REAL maxTetRatio;
			REAL facetAngleTol;
		};

		DelMesher(Vector *surfvs, int *triangls, int numv, int numtri, DelMesherCriteria criteria = DelMesherCriteria());
		DelMesher(const DelMesher& mesher) = delete;
		DelMesher& operator=(const DelMesher& mesher) = delete;
		~DelMesher();
		Reference<Mesh> generateMesh(int *vertexLableMap);

	private:
		bool findSegment(const Segment& s) const{ return meshRep.isSegment(s); }
		Tetrahedron findPosition(Vertex *u, const Tetrahedron& t, const TetMeshDataStructure& meshRep) const;
		Tetrahedron findPositionWithOcclusion(Vertex *u, const Tetrahedron& t, 
			const TetMeshDataStructure& meshRep, const TriMeshDataStructure& occluderRep) const;
		Triangle findPosition(Vertex *u, const Triangle& f) const;
		Triangle findPositionWithOcclusion(Vertex *u, const Triangle& f) const;

		DelVector calculateAbovePoint(const std::vector<Vertex *> &vertices, const Triangle& ref);
		void detectAcuteVertices() const;
		void constrainedTriangulation();
		void triangulation3D(std::vector<Vertex *>& vertices, TetMeshDataStructure& meshRep);

		void constrainedRefinment();
		void refineSubSegment(const Segment& s, Vertex *parent, bool entwinement, bool refineFaceTestEnable, bool skinnyTetTestEnable);
		void refineSubPolygon(TriangleWithGeometry& f, bool skinnyTetTestEnable);
		void refineTetrahedron(const Tetrahedron& tet, bool& processed);
		 
		//segment recovery
		void splitSubSegment(const Segment& s, Vertex* ref, bool missingFacetTest);
		//face recovery
		bool facetRecovery(Triangle& f, std::vector<Vertex *>& regionVertices,
			std::vector<Segment>& regionBoundaries, std::vector<Triangle>& regionFaces,
			std::vector<Vertex *>& positiveVertices, std::vector<Triangle>& positiveFaces,
			std::vector<Vertex *>& negativeVertices, std::vector<Triangle>& negativeFaces,
			std::vector<Tetrahedron>& deleted, std::vector<Tetrahedron>& inserted, 
			bool refineFacetTest = false, bool skinnyTest = false);
		void findMissingRegion(const Triangle& missed, std::vector<Vertex *> &regionVertices, 
			std::vector<Segment>& regionBoundary, std::vector<Triangle>& regionFaces);
		void propagateFindRegion(const Segment& edge, std::vector<Vertex *> &regionVertices,
			std::vector<Segment>& regionBoundary, std::vector<Triangle>& regionFaces, int depth);
		bool findCavity(const std::vector<Segment>& regionBoundary, const std::vector<Triangle>& regionFaces,
			std::vector<Vertex *>& positiveVertices, std::vector<Triangle>& positiveFaces,
			std::vector<Vertex *>& negativeVertices, std::vector<Triangle>& negativeFaces,
			std::vector<Tetrahedron>& deleted, bool refineFacetTest);
		bool findCrossEdge(const Segment& boundary, const std::vector<Triangle>& regionFaces, Segment& cross) const;
		bool triangulateCavity(const std::vector<Triangle>& regionFaces, std::vector<Triangle>& boundaryFaces,
			std::vector<Vertex *>& cavityVertices, std::vector<Tetrahedron>& deleted, std::vector<Tetrahedron>& inserted, Triangle& encroached);
		void propagateCleanCavity(const Triangle& f, TetMeshDataStructure& cavityRep, int depth);
		void refineRegion(const Triangle& regionFacet, bool missingSegTest, bool refineSegTest, bool setInsertionRadius = false, bool skinnyTest = false);

		Vertex* allocVertex(const DelVector &vert, REAL weight, VertexType type = VertexType(VertexType::Vertex_Undefined));
		Vertex* allocVertex(const Vertex &vert);
		void deallocVertex(Vertex *vert);

		struct VolumeVertexInsertionFlags {
			VolumeVertexInsertionFlags() {
				cdt = false;
				skinnyTetTest = false;
				trueInsertion = true;
				refineSegTest = false;
				missingSegTest = false;
				refineFacetTest = false;
				missingFacetTest = false;
				rejected = false;
				insertRadiusTest = false;
				parent = NULL;
			}
			bool cdt, skinnyTetTest, trueInsertion, refineSegTest, missingSegTest, refineFacetTest, missingFacetTest, insertRadiusTest;
			mutable bool rejected;
			mutable Vertex *parent;
		};
		struct SurfaceVertexInsertionFlags {
			SurfaceVertexInsertionFlags() {
				trueInsertion = true;
				deletionRecord = false;
				refineTest = false;
				missingTest = false;
			}
			bool trueInsertion, deletionRecord, refineTest, missingTest;
		};
		void insertVertex(Vertex *u, const Tetrahedron& tet, TetMeshDataStructure& meshRep, 
			const VolumeVertexInsertionFlags& vifs = VolumeVertexInsertionFlags(), Tetrahedron *rt = NULL);
		void insertVertexOnSegment(Vertex *u, const Segment& s, const Tetrahedron& tet, 
			TetMeshDataStructure& meshRep, const VolumeVertexInsertionFlags& vifs);
		void insertVertexOnSurface(Vertex *u, const Triangle& f, const Tetrahedron& tet,
			TetMeshDataStructure& meshRep, const VolumeVertexInsertionFlags& vifs);
		//Caution: Do not dealloc the deleted vertex for memory reuse
		void deleteVertex(Vertex *u, TetMeshDataStructure& meshRep, const VolumeVertexInsertionFlags& vifs = VolumeVertexInsertionFlags());
		void insertSurfaceVertex(Vertex *u, const Triangle &f, const SurfaceVertexInsertionFlags& vifs);
		void insertSurfaceSegmentVertex(Vertex *u, const Segment &s, const SurfaceVertexInsertionFlags& vifs);

		void digCavity(Vertex *u, const Triangle& f, TetMeshDataStructure& meshRep,
			const VolumeVertexInsertionFlags& vifs, int depth, Tetrahedron *rt = NULL);
		void digCavity(Vertex *u, const Segment &s, const DelVector& above, int index, const SurfaceVertexInsertionFlags& vifs);
		void triangulateCavity(Vertex *u, std::vector<Triangle>& boundaries, TetMeshDataStructure& meshRep, const VolumeVertexInsertionFlags& vifs);

		bool Encroached(const Segment &s, Vertex **encroachedVert = NULL) const;
		bool Encroached(TriangleWithGeometry& f, Vertex **encroachedVert) const;
		bool Encroached(const Segment &s, Vertex *v) const;
		bool Encroached(TriangleWithGeometry &f, Vertex *v) const;
		bool Encroached(const DelVector &orthocenter, REAL radius, Vertex *v) const;

		size_t getPolygonVertices(int facetIndex, Vertex ***verts) const;
		bool Adjacent(const Segment &s, Vertex *v) const;
		bool Adjacent(const Segment &s, int facetIndex) const;
		bool Adjacent(int facetIndex, Vertex *v) const;
		bool mayIneligible(int facetIndex, Vertex *v) const;

		bool skinnyTetTest(Tetrahedron& t) const;

		bool findIntersectedTetrahedron(Vertex *a, const DelVector& bb, Tetrahedron *t) const;
		Vertex* findSegmentEncroachedReference(Vertex *end, const Tetrahedron& intersected) const;

		void detectCoplanarFaces(const Triangle& f, REAL facetRadianTol,
			std::vector<Vertex *>& coplanarVertices, std::vector<Segment>& boundaries, TriMeshDataStructure& surfRep) const;
		void propagateDetectCoplanarFaces(Vertex *ref, const Segment& s, REAL facetRadianTol,
			std::vector<Vertex *>& coplanarVertices, std::vector<Segment>& boundaries, TriMeshDataStructure& surfRep, int depth) const;
		void makeHole(Vertex *u, std::vector<Triangle>& boundaries, TetMeshDataStructure& meshRep) const;
		void propagateMakeHole(const Triangle& f, std::vector<Triangle>& boundaries, TetMeshDataStructure& meshRep, int depth) const;

		void propagateClean(const Triangle &f, int depth);

		static Predicator<REAL> predicator;

		std::deque<Tetrahedron> maySkinnyTets;
		std::deque<Triangle> mayMissingFacets;
		std::deque<TriangleWithGeometry> mayRefineFacets;
		std::deque<Segment> mayRefineSegs, mayMissingSegs;
		std::vector<Segment> markedSegments;
		std::vector<uintptr_t *> verticesPerPolygon;
		std::vector<DelVector> abovePoints;

		std::vector<TriangleWithIndex> tobeDeletedFaces;
		std::vector<SegmentWithIndex> newSegsOfFaces;

		std::vector<Tetrahedron> tobeDeletedTets;
		std::vector<Triangle> newFacesOfTets;

		std::vector<Vertex *> oriVertices;
		std::vector<Segment> oriSegments;
		MemoryPool<Vertex> vertPool;
		MemoryArena<uintptr_t> pointerArena;

		TetMeshDataStructure meshRep;
		TriMeshDataStructure surfaceRep;

		DelMesherCriteria criteria;
		AABB<REAL> boundBox;
	};

	inline Vertex* DelMesher::allocVertex(const DelVector &vert, REAL weight, VertexType type){
		Vertex *newVertex = meshRep.allocVertex(vert, weight, type);
		Assert(boundBox.Inside(newVertex->point));
		return newVertex;
	}

	inline Vertex* DelMesher::allocVertex(const Vertex &vertex){
		Vertex *newVertex = meshRep.allocVertex(vertex.point, vertex.weight, vertex.getVertexType());
		Assert(boundBox.Inside(newVertex->point));
		return newVertex;
	}

	inline void DelMesher::deallocVertex(Vertex *vert) {
		meshRep.deallocVertex(vert);
	}

	inline bool DelMesher::skinnyTetTest(Tetrahedron& t) const {
		t.setGeometricProperties();
		return t.getRelaxedRadiusEdgeRation() > criteria.maxTetRatio || t.getRadius() > criteria.maxTetRadius;
	}

	inline size_t DelMesher::getPolygonVertices(int facetIndex, Vertex ***verts) const {
		uintptr_t *p = verticesPerPolygon[facetIndex];
		size_t count = *(size_t *)p;
		*verts = (Vertex **)(p + 1);
		return count;
	}

	inline bool DelMesher::Encroached(const Segment &s, Vertex *v) const {
		return predicator.inDiametricBall(v->point, s.v[0]->point, s.v[1]->point) <= 0;
	}

	inline bool DelMesher::Encroached(const DelVector &circumcenter, REAL radiusSquare, Vertex *v) const {
		constexpr REAL epsilon = REAL(1e-8);
		return !v->isGhost() && ((v->point - circumcenter).length2() - radiusSquare) <= epsilon * epsilon * radiusSquare;
	}
}

#endif




