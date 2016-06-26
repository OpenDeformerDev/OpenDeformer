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
	class DelMesher;

	class DelTriangulator{
	public:
		DelTriangulator(){
			ghost = allocAligned<Vertex>();
			ghost->setGhost();
			invertion = false;
		}
		~DelTriangulator(){
			freeAligned(ghost);
		}
		void generateSubPolygons(Vertex **vertices, int *segments, int vertexCount, int segmentCount, bool boundaryOnly);
		void outPut(DelMesher *mesher);
	private:
		void calculateAbovePoint(int vertexCount, Vertex** vertices);
		void insertSegment(const Segment& s);
		void triangulateHalfHole(const std::vector<Vertex *>& vertices);
		void findCavity(const Segment& s, std::vector<Vertex *>& positive, std::vector<Vertex *>& negtive);
		void insertVertexToCavity(Vertex *u, Vertex *v, Vertex *w,
			bool mark, bool oriTest, std::deque<Vertex *>& marked, int depth);
		void insertVertexToCavity(Vertex *u, Vertex *v, Vertex *w, int depth);
		void triangulateConvexPoly(Vertex *u, const std::vector<Vertex *>& convexPoly);
		Face findPosition(Vertex *u, const Face &f) const;
		void digCavity(Vertex *u, const Segment &s, Face *rf, int depth);
		bool detectInversion(Face &ghostFace) const;
		void cleanRigion(const Face& ghostFace);
		void propagateClean(const Segment& s, int depth);
		inline int findGhost(const Face &f){
			for (int i = 0; i < 3; i++){
				if (f.v[i] == ghost)
					return i;
			}
			return -1;
		}
		inline int findGhost(Vertex **v, int n){
			for (int i = 0; i < n; i++){
				if (v[i] == ghost)
					return i;
			}
			return -1;
		}

		static Predicator<REAL> predicator;
		TriMeshDataStructure meshRep;
		Vertex *ghost;
		bool invertion;

		std::unordered_set<Segment, segment_hash> segments;

		TriMeshDataStructure cavityRep;
	};

	class DelMesher :public Mesher{
	public:
		DelMesher(Vector *survs, int *segis, int *subpolygons, int numv, int numseg, int numpol, int *numsubpol, REAL minRation, REAL minRadius);
		DelMesher(const DelMesher& mesher) = delete;
		DelMesher& operator=(const DelMesher& mesher) = delete;
		~DelMesher() = default;
		Reference<Mesh> generateMesh(int *vertexLableMap);

	private:
		bool findSegment(const Segment& s) const{ return segments.find(s) != segments.end(); }

		Tetrahedron findPosition(Vertex *u, const Tetrahedron& t, bool boundary = false) const;
		Face findPosition(Vertex *u, const Face& f) const;

		void triangulation3D();
	
		void splitSubSegment(const Segment& s, Vertex* ref, bool onFace);

		void splitSubSegment(const Segment& s);
		void splitSubPolygon(const Face& f);
		void splitTetrahedron(const Tetrahedron& tet);

		//face recovery
		enum CavityTriangleType { TwoPositive, TwoNegative, PositiveNegativeCoplanar, NegativePositiveCoplanar };
		void findMissingRegion(const Segment& edge, std::vector<Vertex *> &regionVertices, std::vector<Face>& regionFaces, 
			Segment *boundary, int depth);
		void findCavity(const std::vector<Face>& regionFaces,
			const Face& intersect, const CavityTriangleType& faceType,
			std::vector<Vertex *>& positiveVertices, std::deque<Face>& positiveFaces, 
			std::vector<Vertex *>& negativeVertices, std::deque<Face>& negativeFaces, 
			std::vector<Tetrahedron>& deleted, int depth);
		bool findCrossEdge(const Segment& boundary, const std::vector<Face>& regionFaces, Segment& cross) const;
		bool faceRecovery(Face& f, std::vector<Face>& regionFaces, std::vector<Vertex *>& regionVertices, 
			std::vector<Vertex *>& positiveVertices, std::deque<Face>& positiveFaces,
			std::vector<Vertex *>& negativeVertices, std::deque<Face>& negativeFaces, 
			std::vector<Tetrahedron>& deleted, std::vector<Tetrahedron>& inserted);
		bool triangulateCavity(const std::vector<Vertex *>& regionVertices, const std::vector<Face>& regionFaces,
			const std::vector<Vertex *>& boundaryVertices, std::deque<Face>& boundaryFaces, 
			std::vector<Tetrahedron>& deleted, std::vector<Tetrahedron>& inserted, Face& encroached);
		void insertVertexToCavity(Vertex *u, const Face& f);
		Vertex *allocCavityVertex(const Vertex &vert);
		void deallocCavityVertex(Vertex *vert);

		Vertex* allocVertex(const DelVector &vert, REAL weight = 0.f);
		Vertex* allocVertex(const Vertex &vert);
		void insertVertex(Vertex *u, const Tetrahedron& tet, Tetrahedron *rt = NULL, bool insertToSkinny = false);
		void insertSurfaceVertex(Vertex *u, const Face &f, bool insertToQueue = true);

		void digCavity(Vertex *u, const Face& f, Tetrahedron *rt = NULL, bool insertToSkinny = false, bool trulyDeleteOrAdd = true);
		void digCavity(Vertex *u, const Vertex& aboveVert, const Segment &f, bool insertToQueue = true, bool trulyDeleteOrAdd = true);

		////daling segments need to be fixed
		bool Encroached(const Segment &s) const;
		bool Encroached(const Face &f) const;
		bool Encroached(const Segment &s, const Vertex &v) const;
		bool Encroached(const Face &f, const Vertex &v) const;

		void Protect();
		Vertex* Cover(const Segment &s);
		REAL estimateLocalGapSize2(const DelVector &c) const;

		Face findFaceAroundOnPlane(const Vertex& origin, Vertex *center) const;
		bool findIntersectedTetrahedron(Vertex *a, const DelVector& bb, Tetrahedron *t) const;
		Vertex* findSegmentEncroachedReference(Vertex *end, const Tetrahedron& intersected) const;

		REAL maxRatio;
		REAL maxRadius;
		static Predicator<REAL> predicator;

		std::unordered_set<Segment, segment_hash> segments;
		std::vector<Segment> oriSegments;
		std::vector<int> oriVertexIndices;
		std::unordered_map<Vertex *, int, vertex_hash> vorisHash;
		//unordered_set<Face, face_hash, std::equal_to<Face>, NONE_SYNC_CHUNK_ALLOC<Face>> polygons;
		//unordered_set<Tetrahedron, tet_hash, std::equal_to<Tetrahedron>, NONE_SYNC_CHUNK_ALLOC<Tetrahedron>> tets;

		//unordered_map<Face, Vertex *, face_hash, std::equal_to<Face>, NONE_SYNC_CHUNK_ALLOC<std::pair<Face, Vertex*>>> fvHash;
		//unordered_map<Segment, Vertex *, segment_hash, std::equal_to<Segment>, NONE_SYNC_CHUNK_ALLOC<std::pair<Segment, Vertex*>>> svHash;
		//unordered_map<Vertex *, Segment, vertex_hash> vfHash;
		//unordered_map<Vertex *, Vertex *, vertex_hash> vsHash;

		std::priority_queue<Tetrahedron, std::vector<Tetrahedron>> skinnyTets;
		std::deque<Face> mayEncroachedFaces;
		std::deque<Segment> mayEncroachedSegs;

		std::vector<Face> tobeDeletedFaces;
		std::vector<Segment> newSegsOfFaces;

		std::vector<Tetrahedron> tobeDeletedTets;
		std::vector<Face> newFacesOfTets;

		std::vector<Vertex *> vertices;
		MemoryArena<Vertex> vertArena;
		MemoryPool<Vertex> vertPool;

		TetMeshDataStructure meshRep;
		TetMeshDataStructure cavityRep;
		TriMeshDataStructure surfaceRep;

		Vertex *ghost;
		AABB<REAL> boundBox;
		friend class DelTriangulator;
	};

	inline Vertex* DelMesher::allocVertex(const DelVector &vert, REAL weight){
		Vertex *newVertex = vertArena.Alloc();
		newVertex->vert = vert;
		newVertex->weight = weight;
		newVertex->setLabel();
		Assert(boundBox.Inside(newVertex->vert));
		return newVertex;
	}

	inline Vertex* DelMesher::allocVertex(const Vertex &vertex){
		Vertex *newVertex = vertArena.Alloc();
		*newVertex = vertex;
		newVertex->setLabel();
		Assert(boundBox.Inside(newVertex->vert));
		return newVertex;
	}

	inline Vertex* DelMesher::allocCavityVertex(const Vertex &vert) {
		Vertex *newVertex = vertPool.Alloc();
		vert.Copy(*newVertex);
		return newVertex;
	}

	inline void DelMesher::deallocCavityVertex(Vertex *vert) {
		vertPool.Dealloc(vert);
	}


	inline REAL interiorAngle(const DelVector& o, const DelVector& a, const DelVector& b) {
		DelVector oa = a - o;
		DelVector ob = b - o;
		REAL len = oa.length() * ob.length();
		Assert(len != 0);

		return acos(Clamp((oa * ob) / len, REAL(-1), REAL(1)));
	}

}

#endif




