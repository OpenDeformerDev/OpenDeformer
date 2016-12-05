#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_MESHER_DELMESHER_H
#define ODER_MESHER_DELMESHER_H

#include "predicate.h"
#include "delPrimitive.h"
#include "memory.h"
#include "mesher.h"
#include "tetmesh.h"
#include "datastructure.h"
#include "aabb.h"
#include "meshRelabeler.h"
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <deque>

namespace ODER{
	class DelTriangulator{
	public:
		DelTriangulator(){
			ghost = allocAligned<Vertex>();
			ghost->setGhost();
		}
		~DelTriangulator(){
			freeAligned(ghost);
		}
		void generateSubPolygons(Vertex **vertices, int *segments, int vertexCount, int segmentCount, const Face& ref, bool boundaryOnly);
		void generateSubPolygons(Vertex **vertices, Segment *segments, int vertexCount, int segmentCount, const Face& ref, bool boundaryOnly);
		void insertSegments(const Face *triangles, const Segment *segments, int triangleCount, int segmentCount);
		void outPut(std::vector<Face> &meshVec) { meshRep.getTriangles(false, meshVec); }
	private:
		void calculateAbovePoint(int vertexCount, Vertex** vertices, const Face& ref);
		void insertSegment(const Segment& s);
		void triangulateHalfHole(const std::vector<Vertex *>& vertices);
		void findCavity(const Segment& s, std::vector<Vertex *>& positive, std::vector<Vertex *>& negtive);
		void insertVertexToCavity(Vertex *u, Vertex *v, Vertex *w,
			bool mark, bool oriTest, std::deque<Vertex *>& marked, int depth);
		void insertVertexToCavity(Vertex *u, Vertex *v, Vertex *w, int depth);
		void triangulateConvexPoly(Vertex *u, const std::vector<Vertex *>& convexPoly);
		Face findPosition(Vertex *u, const Face &f) const;
		void digCavity(Vertex *u, const Segment &s, Face *rf, int depth);
		void propagateClean(const Segment& s, int depth);
		inline int findGhost(const Face &f){
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

		std::unordered_set<Segment, segment_unordered_hash> segments;

		TriMeshDataStructure cavityRep;
	};

	class DelMesher :public Mesher{
	public:
		DelMesher(Vector *surfvs, int *triangls, int numv, int numtri, REAL maxRation, REAL maxRadius, REAL facetAngleTol = REAL(179));
		DelMesher(const DelMesher& mesher) = delete;
		DelMesher& operator=(const DelMesher& mesher) = delete;
		~DelMesher();
		Reference<Mesh> generateMesh(int *vertexLableMap);

	private:
		bool findSegment(const Segment& s) const{ return segments.find(s) != segments.end(); }
		Tetrahedron findPosition(Vertex *u, const Tetrahedron& t, const TetMeshDataStructure& meshRep) const;
		Face findPosition(Vertex *u, const DelVector& above, const Face& f) const;

		void constrainedTriangulation();
		void triangulation3D(std::vector<Vertex *>& vertices, TetMeshDataStructure& meshRep, bool insertToSkinny);
	
		void splitSubSegment(const Segment& s, Vertex* ref);

		void splitSubSegment(const Segment& s);
		void splitSubPolygon(const Face& f);
		void splitTetrahedron(const Tetrahedron& tet);

		//face recovery
		bool faceRecovery(Face& f, std::vector<Vertex *>& regionVertices,
			std::vector<Segment>& regionBoundaries, std::vector<Face>& regionFaces,
			std::vector<Vertex *>& positiveVertices, std::vector<Face>& positiveFaces,
			std::vector<Vertex *>& negativeVertices, std::vector<Face>& negativeFaces,
			std::vector<Tetrahedron>& deleted, std::vector<Tetrahedron>& inserted);
		void findMissingRegion(const Face& missed, std::vector<Vertex *> &regionVertices, 
			std::vector<Segment>& regionBoundary, std::vector<Face>& regionFaces);
		void propagateFindRegion(const Segment& edge, std::vector<Vertex *> &regionVertices,
			std::vector<Segment>& regionBoundary, std::vector<Face>& regionFaces, int depth);
		bool findCavity(const std::vector<Segment>& regionBoundary, std::vector<Face>& regionFaces,
			std::vector<Vertex *>& positiveVertices, std::vector<Face>& positiveFaces,
			std::vector<Vertex *>& negativeVertices, std::vector<Face>& negativeFaces,
			std::vector<Tetrahedron>& deleted);
		bool findCrossEdge(const Segment& boundary, const std::vector<Face>& regionFaces, Segment& cross) const;
		bool triangulateCavity(const std::vector<Face>& regionFaces, int encroachFaceTest, std::vector<Face>& boundaryFaces,
			std::vector<Vertex *>& cavityVertices, std::vector<Tetrahedron>& deleted, std::vector<Tetrahedron>& inserted, Face& encroached);
		void propagateCleanCavity(const Face& f, int depth);
		void refineRegion(const Face& regionFace);

		Vertex* allocVertex();
		Vertex* allocVertex(const DelVector &vert, REAL weight, VertexType type);
		Vertex* allocVertex(const Vertex &vert);
		void deallocVertex(Vertex *vert);

		struct VertexInsertionFlags {
			VertexInsertionFlags() {
				cdt = false;
				skinnyTetTest = false;
				trueInsertion = true;
				encroachSegTest = false;
				encroachFaceTest = 0;
			}
			bool cdt, skinnyTetTest, trueInsertion, encroachSegTest;
			int encroachFaceTest; //1 for missing test only; 2 for full encroachment test
		};
		void insertVertex(Vertex *u, const Tetrahedron& tet, TetMeshDataStructure& meshRep, 
			const VertexInsertionFlags& vifs = VertexInsertionFlags(), Tetrahedron *rt = NULL);
		void insertSurfaceVertex(Vertex *u, const Face &f, bool insertToQueue = true);
		void insertSurfaceSegmentVertex(Vertex *u, const Segment &s, Face *inFace = NULL, bool insertToQueue = true);
		bool digCavity(Vertex *u, const Face& f, TetMeshDataStructure& meshRep,
			const VertexInsertionFlags& vifs, Tetrahedron *rt = NULL);
		void digCavity(Vertex *u, const DelVector& above, const Segment &s, int index, bool insertToQueue = true, bool trulyDeleteOrAdd = true);

		////daling segments need to be fixed
		bool Encroached(const Segment &s) const;
		bool Encroached(const Face &f) const;
		bool Encroached(const Segment &s, const Vertex &v) const;
		bool Encroached(const Face &f, const Vertex &v) const;

		bool findIntersectedTetrahedron(Vertex *a, const DelVector& bb, Tetrahedron *t) const;
		Vertex* findSegmentEncroachedReference(Vertex *end, const Tetrahedron& intersected) const;

		void detectCoplanarFaces(const Face& f, REAL facetRadianTol,
			std::vector<Vertex *>& coplanarVertices, std::vector<Segment>& boundaries, TriMeshDataStructure& surfRep) const;
		void propagateDetectCoplanarFaces(Vertex *ref, const Segment& s, REAL facetRadianTol,
			std::vector<Vertex *>& coplanarVertices, std::vector<Segment>& boundaries, TriMeshDataStructure& surfRep, int depth) const;

		void propagateClean(const Face &f, int depth);

		REAL maxRatio;
		REAL maxRadius;
		static Predicator<REAL> predicator;

		std::unordered_set<Segment, segment_ordered_hash> segments;
		std::unordered_map<Vertex *, Segment, vertex_hash> vertSegMap;
		//unordered_set<Face, face_hash, std::equal_to<Face>, NONE_SYNC_CHUNK_ALLOC<Face>> polygons;
		//unordered_set<Tetrahedron, tet_hash, std::equal_to<Tetrahedron>, NONE_SYNC_CHUNK_ALLOC<Tetrahedron>> tets;

		//unordered_map<Face, Vertex *, face_hash, std::equal_to<Face>, NONE_SYNC_CHUNK_ALLOC<std::pair<Face, Vertex*>>> fvHash;
		//unordered_map<Segment, Vertex *, segment_hash, std::equal_to<Segment>, NONE_SYNC_CHUNK_ALLOC<std::pair<Segment, Vertex*>>> svHash;
		//unordered_map<Vertex *, Segment, vertex_hash> vfHash;
		//unordered_map<Vertex *, Vertex *, vertex_hash> vsHash;

		std::priority_queue<Tetrahedron, std::vector<Tetrahedron>, std::function<bool(const Tetrahedron&, const Tetrahedron&)>> skinnyTets;
		std::deque<Face> mayEncroachedFaces;
		std::deque<Segment> mayEncroachedSegs;
		std::unordered_set<Segment, segment_ordered_hash> mayEncroachedSegsSet;
		std::vector<uintptr_t *> polygonVertices;

		std::vector<Face> tobeDeletedFaces;
		std::vector<Segment> newSegsOfFaces;

		std::vector<Tetrahedron> tobeDeletedTets;
		std::vector<Tetrahedron> newTets;

		std::vector<Vertex *> oriVertices;
		MemoryPool<Vertex> vertPool;

		TetMeshDataStructure meshRep;
		TetMeshDataStructure cavityRep;
		TriMeshDataStructure surfaceRep;

		Vertex *ghost;
		AABB<REAL> boundBox;
		VertexLabeler labeler;
	};

	inline Vertex* DelMesher::allocVertex() { return vertPool.Alloc(); }

	inline Vertex* DelMesher::allocVertex(const DelVector &vert, REAL weight, VertexType type){
		Vertex *newVertex = vertPool.Alloc();
		newVertex->vert = vert;
		newVertex->weight = weight;
		newVertex->type = type;
		newVertex->setLabel(labeler.getLabel());
		Assert(boundBox.Inside(newVertex->vert));
		return newVertex;
	}

	inline Vertex* DelMesher::allocVertex(const Vertex &vertex){
		Vertex *newVertex = vertPool.Alloc();
		*newVertex = vertex;
		newVertex->setVertexPointer(NULL);
		Assert(boundBox.Inside(newVertex->vert));
		return newVertex;
	}

	inline void DelMesher::deallocVertex(Vertex *vert) {
		vertPool.Dealloc(vert);
	}
}

#endif




