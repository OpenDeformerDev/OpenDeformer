#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_MESHER_DELPRIMITIVE_H
#define ODER_MESHER_DELPRIMITIVE_H

#include "oder.h"
#include "latool.h"
#include "memory.h"
#include "allocator.h"
#include <numeric>

#define NEXT_F(i) (((i)+1)%3)
#define NEXT_T(i) (((i)+1)%4)
#define SQRT_2 1.4142135623730950488016887242096980785696718753769480732

namespace ODER {
	using DelReal = double;
	using DelVector = VectorBase<DelReal>;

	class VertexLabeler {
	public:
		VertexLabeler() : label(1) {}
		int getLabel() {
			return label++;
		}
		static constexpr int getSpecilGhostLabel() {
			return 0;
		}
		void restartLable() {
			label = 1;
		}
	private:
		int label;
	};

	enum VertexType {
		Vertex_Undefined = 0,
		Vertex_Fixed = 1 << 0,
		Vertex_Free = 1 << 1,
		Vertex_Segment = 1 << 2,
		Vertex_Facet = 1 << 3,
		Vertex_Volume = 1 << 4,
		Vertex_Acute = 1 << 5,
		Vertex_Refined = 1 << 6,

		Vertex_FixedSegment = Vertex_Fixed | Vertex_Segment,
		Vertex_FixedFacet = Vertex_Fixed | Vertex_Facet,
		Vertex_FixedVolume = Vertex_Fixed | Vertex_Volume,

		Vertex_FreeSegment = Vertex_Free | Vertex_Segment,
		Vertex_FreeFacet = Vertex_Free | Vertex_Facet,
		Vertex_FreeVolume = Vertex_Free | Vertex_Facet,

		Vertex_RefinedFreeSegment = Vertex_FreeSegment | Vertex_Refined,
		Vertex_RefinedFreeFacet = Vertex_FreeFacet | Vertex_Refined,
		Vertex_RefinedFreeVolume = Vertex_FreeVolume | Vertex_Refined,

		Vertex_LowDimension = Vertex_Segment | Vertex_Facet,
	};

	class EdgeListNode;
	class TriVertexListNode;
	class TriMeshDataStructure;
	class TetMeshDataStructure;

	struct Vertex {
		Vertex() : weight(0), relaxedInsetionRadius(std::numeric_limits<DelReal>::max()), label(-1), 
			vertexPointer(NULL), pointers(0), type(VertexType::Vertex_Undefined){}
		template<class FT> explicit Vertex(const VectorBase<FT>& p, VertexType t = VertexType::Vertex_Undefined)
			: point{ p.x, p.y, p.z }, weight(0), relaxedInsetionRadius(std::numeric_limits<DelReal>::max()), label(-1), 
			vertexPointer(NULL), pointers(0), type(t) {}
		explicit Vertex(const DelVector& p, DelReal w = 0, VertexType t = VertexType::Vertex_Undefined)
			: point(p), weight(w), relaxedInsetionRadius(std::numeric_limits<DelReal>::max()), label(-1), pointers(0), vertexPointer(NULL), type(t) {}
		void setGhost() {
			constexpr DelReal inf = std::numeric_limits<DelReal>::infinity();
			point.x = inf; point.y = inf; point.z = inf;
			weight = -inf;
			label = VertexLabeler::getSpecilGhostLabel();
			pointers = 0;
			vertexPointer = NULL;
			relaxedInsetionRadius = inf;
		}
		bool isGhost() const {
			return label == VertexLabeler::getSpecilGhostLabel();
		}
		void setLabel(int label) {
			this->label = label;
		}
		int getLabel() const {
			return label;
		}
		void setMark() {
			pointers |= 0x1;
		}
		void unSetMark() {
			pointers &= (~0x1);
		}
		bool isMarked() const {
			return (pointers & 0x1) == 0x1;
		}
		void setVertexPointer(Vertex *v) {
			vertexPointer = v;
		}
		Vertex *getPointedVertex() const {
			return vertexPointer;
		}
		VertexType getVertexType() const { return type; }
		void setAcute() { type = VertexType(type | VertexType::Vertex_Acute); }
		void setOriSegmentIndex(int index);
		int getOriSegmentIndex() const;

		DelVector point;
		DelReal weight;
		DelReal relaxedInsetionRadius;
	private:
		void setEdgeList(EdgeListNode *n);
		EdgeListNode *getEdgeList();
		bool hasEdgeList() const;
		void setSupplyVertexPointer(Vertex *v);
		Vertex *getSupplyVertex() const;
		void setFaceLink(TriVertexListNode *link);
		TriVertexListNode *getFaceLink() const;
		void setEnforcedEdgeList(EdgeListNode *node);
		EdgeListNode *getEnforcedEdgeList() const;
		void setEnforcedEdgeMark();
		bool isEnforcedEdgeMarked() const;
		void setPointers(uintptr_t *p);
		uintptr_t *getPointers() const;

		int label;
		VertexType type;
		uintptr_t pointers;
		Vertex *vertexPointer;

		friend class TriMeshDataStructure;
		friend class TetMeshDataStructure;
	};

	class TriVertexListNode {
	public:
		TriVertexListNode() :currentVert(0x1), next(NULL), faceIndex(-1) {}
		void setVertex(Vertex *vert) {
			currentVert = uintptr_t(vert);
		}
		Vertex *getVertex() const {
			return (Vertex *)(currentVert & (~0x1));
		}
		TriVertexListNode *getNextNode() const {
			return (TriVertexListNode *)(next & (~0x3));
		}
		void setNextNode(TriVertexListNode *nextNode) {
			next = uintptr_t(nextNode) | (next & 0x3);
		}
		void setDeletedMark() {
			currentVert |= 0x1;
		}
		void unSetDeletedMark() {
			currentVert &= (~0x1);
		}
		bool isPreFaceDeleted() const {
			return (currentVert & 0x1) == 0x1;
		}
		void setMark() {
			next |= 0x1;
		}
		void unSetMark() {
			next &= (~0x1);
		}
		bool isMarked() {
			return (next & 0x1) == 0x1;
		}
		void setEdgeMark() {
			next |= 0x2;
		}
		void unSetEdgeMark() {
			next &= (~0x2);
		}
		bool isEdgeMarked() const {
			return (next & 0x2) == 0x2;
		}
		void setIndex(int index) {
			faceIndex = index;
		}
		int getIndex() const { 
			return faceIndex; 
		}
	private:
		uintptr_t currentVert;
		uintptr_t next;
		int faceIndex;
	};


	class TetVertexListNode {
	public:
		TetVertexListNode() :currentVert(0x1), next(NULL) {}
		void setVertex(Vertex *vert) {
			currentVert = uintptr_t(vert);
		}
		Vertex *getVertex() const {
			return (Vertex *)(currentVert & (~0x1));
		}
		TetVertexListNode *getNextNode() const {
			return (TetVertexListNode *)(next & (~0x3));
		}
		void setNextNode(TetVertexListNode *nextNode) {
			next = uintptr_t(nextNode) | (next & 0x3);
		}
		void setDeletedMark() {
			currentVert |= 0x1;
		}
		void unSetDeletedMark() {
			currentVert &= (~0x1);
		}
		bool isPreFaceDeleted() const {
			return (currentVert & 0x1) == 0x1;
		}
		void setForwardMark() {
			next |= 0x1;
		}
		void unSetForwardMark() {
			next &= (~0x1);
		}
		bool isForwardMarked() {
			return (next & 0x1) == 0x1;
		}
		void setBackwardMark() {
			next |= 0x2;
		}
		void unSetBackwardMark() {
			next &= (~0x2);
		}
		bool isBackwardMarked() {
			return (next & 0x2) == 0x2;
		}
	private:
		uintptr_t currentVert;
		uintptr_t next;
	};

	class EdgeListNode {
	public:
		EdgeListNode() :endVertex(NULL), link(NULL), next(NULL) {}
		void setEndVertex(Vertex *vert) {
			endVertex = vert;
		}
		Vertex *getEndVertex() const {
			return endVertex;
		}
		void setLink(TetVertexListNode *Link) {
			link = uintptr_t(Link) | (link & 0x1);
		}
		TetVertexListNode* getLink() const {
			return (TetVertexListNode*)(link & (~0x1));
		}
		EdgeListNode *getNextNode() const {
			return next;
		}
		void setNextNode(EdgeListNode *nextNode) {
			next = nextNode;
		}
		void setOtherVertex(Vertex *vert) {
			link = uintptr_t(vert) | (link & 0x1);
		}
		Vertex *getOtherVertex() {
			return (Vertex *)(link & (~0x1));
		}
		void setMark() {
			link |= 0x1;
		}
		void unSetMark() {
			link &= (~0x1);
		}
		bool isMarked() const {
			return (link & 0x1) == 0x1;
		}
	private:
		Vertex *endVertex;
		uintptr_t link;
		EdgeListNode *next;
	};

	struct Segment {
		Segment() {
			v[0] = v[1] = NULL;
		}
		Segment(Vertex *v0, Vertex *v1, bool ordered = false) {
			if (ordered && (v0->getLabel() > v1->getLabel())) std::swap(v0, v1);
			v[0] = v0;
			v[1] = v1;
		}
		bool operator==(const Segment &s) const {
			return v[0] == s.v[0] && v[1] == s.v[1];
		}
		DelReal getLength() const { return (v[0]->point - v[1]->point).length(); }
		Vertex *v[2];
	};

	struct SegmentWithIndex : public Segment {
		SegmentWithIndex() : Segment(), index(-1) {}
		SegmentWithIndex(Vertex *v0, Vertex *v1, int index, bool ordered = false)
			: Segment(v0, v1, ordered), index(index) {}

		int index;
	};

	struct Triangle {
		Triangle() {
			v[0] = v[1] = v[2] = NULL;
		}
		Triangle(Vertex *v0, Vertex *v1, Vertex *v2, bool ordered = false) {
			initVertices(v0, v1, v2, ordered);
		}
		void sortVertices();
		bool operator==(const Triangle& f) const {
			return v[0] == f.v[0] && v[1] == f.v[1] && v[2] == f.v[2];
		}
		bool operator!=(const Triangle &f) const {
			return v[0] != f.v[0] || v[1] != f.v[1] || v[2] != f.v[2];
		}
		Vertex *v[3];
	private:
		void initVertices(Vertex *v0, Vertex *v1, Vertex *v2, bool ordered);
	};

	struct TriangleWithIndex : public Triangle {
		TriangleWithIndex() : Triangle(), index(-1) {}
		TriangleWithIndex(Vertex *v0, Vertex *v1, Vertex *v2, int index, bool ordered = false)
			: Triangle(v0, v1, v2, ordered), index(index) {}

		int index;
	};

	struct TriangleWithGeometry : public Triangle {
	public:
		TriangleWithGeometry() : Triangle(), radius(-1) {}
		TriangleWithGeometry(Vertex *v0, Vertex *v1, Vertex *v2, bool ordered = false) : 
			Triangle(v0, v1, v2, ordered), radius(-1) {}
		TriangleWithGeometry(const Triangle& f): Triangle(f), radius(-1) {}
		void setGeometricProperties();
		DelReal getRadius() const { return radius; }
		DelVector getCircumcenter() const { return circumcenter; }
	private:
		DelReal radius;
		DelVector circumcenter;
	};

	struct Tetrahedron {
		Tetrahedron() {
			v[0] = v[1] = v[2] = v[3] = NULL;
			reRation = radius = DelReal(-1);
		}
		Tetrahedron(Vertex *v0, Vertex *v1, Vertex *v2, Vertex *v3, bool ordered = false) {
			reRation = radius = DelReal(-1);
			v[0] = v0; v[1] = v1; v[2] = v2; v[3] = v3;
			if (ordered) sortVertices();
		}
		void setGeometricProperties();
		void sortVertices();
		bool operator==(const Tetrahedron& t) const {
			return v[0] == t.v[0] &&
				v[1] == t.v[1] &&
				v[2] == t.v[2] &&
				v[3] == t.v[3];
		}
		bool operator!=(const Tetrahedron& t) const {
			return v[0] != t.v[0] ||
				v[1] != t.v[1] ||
				v[2] != t.v[2] ||
				v[3] != t.v[3];
		}
		DelReal getRelaxedRadiusEdgeRation() const { return reRation; }
		DelReal getRadius() const { return radius; }
		DelVector getCircumcenter() const { return circumcenter; }
		Vertex *v[4];
	private:
		DelReal reRation;
		DelReal radius;
		DelVector circumcenter;
	};

	struct face_compare {
		bool operator()(const Triangle &left, const Triangle &right) const {
			if (left.v[2] == right.v[2]) {
				if (left.v[1] == right.v[1])
					return left.v[0] < right.v[0];
				return left.v[1] < right.v[1];
			}
			return left.v[2] < right.v[2];
		}
	};

	struct face_hash {
		size_t operator()(const Triangle& f) const {
			size_t hashVal = 0;
			hashCombine(hashVal, f.v[0]->getLabel(), f.v[1]->getLabel(), f.v[2]->getLabel());
			return hashVal;
		}
	};

	struct segment_ordered_hash {
		size_t operator()(const Segment &s) const {
			constexpr int shift = 8 * sizeof(int);
			size_t smallerLable = s.v[0]->getLabel();
			uint64_t biggerLable = s.v[1]->getLabel();
			uint64_t hashVal = ((biggerLable * (biggerLable + 1)) >> 1) + smallerLable;
			return size_t((hashVal >> shift) ^ hashVal);
		}
	};

	struct segment_unordered_hash {
		size_t operator()(const Segment &s) const {
			size_t hashVal = 0;
			hashCombine(hashVal, s.v[0]->getLabel(), s.v[1]->getLabel());
			return hashVal;
		}
	};

	struct vertex_hash {
		size_t operator()(const Vertex *v) const {
			return v->getLabel();
		}
	};

	struct tet_hash {
		size_t operator()(const Tetrahedron &t) const {
			size_t hashVal = 0;
			hashCombine(hashVal, t.v[0]->getLabel(), t.v[1]->getLabel(), t.v[2]->getLabel(), t.v[3]->getLabel());
			return hashVal;
		}
	};

	struct tet_compare {
		bool operator()(const Tetrahedron &left, const Tetrahedron &right) const {
			if (left.v[3] == right.v[3]) {
				if (left.v[2] == right.v[2]) {
					if (left.v[1] == right.v[1])
						return left.v[0] < right.v[0];
					return left.v[1] < right.v[1];
				}
				return left.v[2] < right.v[2];
			}
			return left.v[3] < right.v[3];
		}
	};

	class TriMeshDataStructure {
	public:
		TriMeshDataStructure();
		Vertex *getGhostVertex();
		Vertex* allocVertex(const DelVector &point, DelReal weight);
		void deallocVertex(Vertex *vert);
		void bindVolumeVertex(Vertex *vert);
		void unbindVolumeVertex(Vertex *vert);
		void addTriangle(Vertex *a, Vertex *b, Vertex *c, int index = -1);
		void deleteTriangle(Vertex *a, Vertex *b, Vertex *c);
		void setDeletedMark(Vertex *u, Vertex *v);
		void unSetDeletedMark(Vertex *u, Vertex *v);
		void setDeletedMark(Vertex *a, Vertex *b, Vertex *c);
		void unSetDeletedMark(Vertex *a, Vertex *b, Vertex *c);
		void setMark(Vertex *u, Vertex *v);
		void unSetMark(Vertex *u, Vertex *v);
		bool isMarked(Vertex *u, Vertex *v) const;
		bool Adjacent(const Segment &s, Vertex **w, int *index = NULL) const;
		bool adjacent2Vertex(Vertex *w, Triangle *f, int *index = NULL) const;
		int getTriangleIndex(Vertex *a, Vertex *b, Vertex *c) const;
		bool Contain(Vertex *v) const;
		bool Contain(const Segment &s) const;
		bool Contain(const Triangle &f) const;
		void Clear();
		void Reserve(size_t n) { vertices.reserve(n); }
		std::vector<Triangle> getTriangles(bool ghost) const;
		void getTriangles(bool ghost, std::vector<Triangle>& triangles) const;

		//Ensure: edge ab or ba must in mesh
		void addSegment(Vertex *a, Vertex *b);
		//Ensure: edge ab or ba must in mesh
		void deleteSegment(Vertex *a, Vertex *b);
		bool isSegment(const Segment& s) const;

		~TriMeshDataStructure();

		class TriMeshConstIterator {
		public:
			using iterator_category = std::input_iterator_tag;
			using value_type = Triangle;
			using difference_type = ptrdiff_t;
			using size_type = size_t;
			using reference = const Triangle&;
			using pointer = const Triangle *;

			TriMeshConstIterator(const std::vector<Vertex *>::const_iterator& iter,
				const std::vector<Vertex *>::const_iterator& end);
			TriMeshConstIterator(const TriMeshConstIterator&) = default;
			TriMeshConstIterator(TriMeshConstIterator&&) = default;
			TriMeshConstIterator& operator=(const TriMeshConstIterator&) = default;
			TriMeshConstIterator& operator=(TriMeshConstIterator&&) = default;

			bool operator==(const TriMeshConstIterator& right) const;
			bool operator!=(const TriMeshConstIterator& right) const;

			TriMeshConstIterator& operator++();
			TriMeshConstIterator operator++(int);

			reference operator*() const { return current; }
			pointer operator->() const { return &current; }

		private:
			void findNext();

			Triangle current;
			TriVertexListNode *parent, *child;
			std::vector<Vertex *>::const_iterator vertIter;
			std::vector<Vertex *>::const_iterator vertEnd;
		};

		class TriMeshConstCirculator {
		public:
			using iterator_category = std::input_iterator_tag;
			using value_type = Triangle;
			using difference_type = ptrdiff_t;
			using size_type = size_t;
			using reference = const Triangle&;
			using pointer = const Triangle *;

			TriMeshConstCirculator(Vertex *v);
			TriMeshConstCirculator(const TriMeshConstCirculator&) = default;
			TriMeshConstCirculator(TriMeshConstCirculator&&) = default;
			TriMeshConstCirculator& operator=(const TriMeshConstCirculator&) = default;
			TriMeshConstCirculator& operator=(TriMeshConstCirculator&&) = default;

			bool operator==(const TriMeshConstCirculator& right) const;
			bool operator!=(const TriMeshConstCirculator& right) const;
			bool operator==(std::nullptr_t p) const;
			bool operator!=(std::nullptr_t p) const;

			reference operator*() const { return current; }
			pointer operator->() const { return &current; }

			TriMeshConstCirculator& operator++();
			TriMeshConstCirculator operator++(int);
		private:
			void findNext();

			Vertex *vert;
			TriVertexListNode *link;
			Triangle current;
		};

		using const_iterator = TriMeshConstIterator;
		const_iterator begin() const { return const_iterator(vertices.begin(), vertices.end()); }
		const_iterator end() const { return const_iterator(vertices.end(), vertices.end()); }
		TriMeshConstCirculator getIncidentTriangles(Vertex *v) const { return TriMeshConstCirculator(v); }
	private:
		bool getAdjacentListNode(Vertex* u, Vertex* v, TriVertexListNode **w) const;
		bool getMarkedNode(Vertex* u, Vertex* v, TriVertexListNode **w) const;
		void insertToTopology(Vertex *a, Vertex *b, Vertex *c, int index);
		void removeFromTopology(Vertex *a, Vertex *b, Vertex *c);
		static bool verticesOrderCheck(const Vertex *a, const Vertex *b, const Vertex *c);

		std::vector<Vertex *> vertices;
		MemoryPool<TriVertexListNode> *nodePool;
		MemoryArena<Vertex> *vertPool;
		VertexLabeler labeler;
		Vertex *ghost;
		Vertex *deadVerticesStack;
	};

	class TetMeshDataStructure {
	public:
		TetMeshDataStructure();
		Vertex* getGhostVertex();
		Vertex* allocVertex(const DelVector &point, DelReal weight, VertexType extraType = VertexType::Vertex_Undefined);
		void deallocVertex(Vertex *vert);
		void addTetrahedron(Vertex *a, Vertex *b, Vertex *c, Vertex *d);
		void deleteTetrahedron(Vertex *a, Vertex *b, Vertex *c, Vertex *d);
		void setMark(Vertex *u, Vertex *v, Vertex *w);
		void unSetMark(Vertex *u, Vertex *v, Vertex *w);
		bool isMarked(Vertex *u, Vertex *v, Vertex *w) const;
		bool Adjacent(const Triangle &f, Vertex **w) const;
		bool adjacent2Vertex(Vertex *w, Tetrahedron *t) const;
		std::vector<Tetrahedron> getTetrahedrons(bool ghost) const;
		void getTetrahedrons(bool ghost, std::vector<Tetrahedron>& tets) const;
		bool Contain(Vertex *v) const;
		bool Contain(const Triangle &f) const;
		bool Contain(const Tetrahedron& t) const;
		void Clear();
		void Reserve(size_t n) { vertices.reserve(n); }

		//Ensure: the segment ab must not in mesh
		void addSegment(Vertex *a, Vertex *b);
		void deleteSegment(Vertex *a, Vertex *b);
		bool adjacent2SegmentFast(const Segment &s, Tetrahedron *t) const;
		bool adjacent2Segment(const Segment &s, Tetrahedron *t) const;
		//Ensure: the segment s must be added
		bool Contain(const Segment &s) const;
		bool isSegment(const Segment &s) const;
		//Mark segment
		void setMark(Vertex *a, Vertex *b);
		void unSetMark(Vertex *a, Vertex *b);
		bool isMarked(Vertex *a, Vertex *b) const;
		bool testAndMark(Vertex *a, Vertex *b);

		//fast adjacency query check
		bool fastVertexQueryCheck(Vertex *u) const;
		bool fastSegmentQueryCheck(const Segment &s) const;

		~TetMeshDataStructure();

		class TetMeshConstIterator {
		public:
			using iterator_category = std::input_iterator_tag;
			using value_type = Tetrahedron;
			using difference_type = ptrdiff_t;
			using size_type = size_t;
			using reference = const Tetrahedron&;
			using pointer = const Tetrahedron *;
		
			TetMeshConstIterator() : linkHead(NULL), parent(NULL), child(NULL) {}
			TetMeshConstIterator(const std::vector<Vertex *>::const_iterator& iter,
				const std::vector<Vertex *>::const_iterator& end);
			TetMeshConstIterator(const TetMeshConstIterator&) = default;
			TetMeshConstIterator(TetMeshConstIterator&&) = default;
			TetMeshConstIterator& operator=(const TetMeshConstIterator&) = default;
			TetMeshConstIterator& operator=(TetMeshConstIterator&&) = default;

			bool operator==(const TetMeshConstIterator& right) const;
			bool operator!=(const TetMeshConstIterator& right) const;

			TetMeshConstIterator& operator++();
			TetMeshConstIterator operator++(int);

			reference operator*() const { return current; }
			pointer operator->() const { return &current; }

		private:
			void findNext();

			Tetrahedron current;
			EdgeListNode *linkHead;
			TetVertexListNode *parent, *child;

			std::vector<Vertex *>::const_iterator vertIter;
			std::vector<Vertex *>::const_iterator vertEnd;
		};

		class SegmentConstIterator {
		public:
			using iterator_category = std::input_iterator_tag;
			using value_type = Segment;
			using difference_type = ptrdiff_t;
			using size_type = size_t;
			using reference = const Segment&;
			using pointer = const Segment *;

			SegmentConstIterator(const std::vector<Vertex *>::const_iterator &iter, const std::vector<Vertex *>::const_iterator &end);
			SegmentConstIterator& operator++();
			SegmentConstIterator operator++(int);

			bool operator==(const SegmentConstIterator& right) const;
			bool operator!=(const SegmentConstIterator& right) const;

			reference operator*() const { return segment; }
			pointer operator->() const { return &segment; }
		private:
			void findNext();
			std::vector<Vertex *>::const_iterator edgeVertIter;
			std::vector<Vertex *>::const_iterator edgeVertEnd;
			EdgeListNode *edgeNode;
			Segment segment;
		};

		class TetMeshFacetConstCirculator {
		public:
			using value_type = Triangle;
			using difference_type = ptrdiff_t;
			using size_type = size_t;
			using reference = const Triangle&;
			using pointer = const Triangle *;

			TetMeshFacetConstCirculator(const Segment& s, const TetMeshDataStructure *ds);
			TetMeshFacetConstCirculator(const TetMeshFacetConstCirculator&) = default;
			TetMeshFacetConstCirculator(TetMeshFacetConstCirculator&&) = default;
			TetMeshFacetConstCirculator& operator=(const TetMeshFacetConstCirculator&) = default;
			TetMeshFacetConstCirculator& operator=(TetMeshFacetConstCirculator&&) = default;

			bool operator==(const TetMeshFacetConstCirculator& right) const;
			bool operator!=(const TetMeshFacetConstCirculator& right) const;
			bool operator==(std::nullptr_t p) const;
			bool operator!=(std::nullptr_t p) const;

			reference operator*() const { return current; }
			pointer operator->() const { return &current; }

			TetMeshFacetConstCirculator& operator++();
			TetMeshFacetConstCirculator operator++(int);
		private:
			void findNext();

			Vertex *ori, *dest;
			EdgeListNode *linkListHead;
			TetVertexListNode *listNode;
			Vertex *initAdjacentVert;
			Vertex *movingVert;
			bool reverse;

			const TetMeshDataStructure *ds;

			Triangle current;
		};

		using const_iterator = TetMeshConstIterator;
		const_iterator begin() const { return const_iterator(vertices.begin(), vertices.end()); }
		const_iterator end() const { return const_iterator(vertices.end(), vertices.end()); }
		SegmentConstIterator segmentBegin() const { return SegmentConstIterator(enforcedEdgeVertices.begin(), enforcedEdgeVertices.end()); }
		SegmentConstIterator segmentEnd() const { return SegmentConstIterator(enforcedEdgeVertices.end(), enforcedEdgeVertices.end()); }
		TetMeshFacetConstCirculator getIncidentFacets(const Segment& s) const { return TetMeshFacetConstCirculator(s, this); }
	private:
		bool getAdjacentListNode(const Triangle& f, TetVertexListNode **z) const;
		bool getAdjacentListNode(Vertex* w, Vertex *x, Vertex *y, TetVertexListNode **z) const;
		int getMarkedNode(const Triangle& f, TetVertexListNode **marked) const;
		int getMarkedNode(Vertex *w, Vertex *x, Vertex *y, TetVertexListNode **marked) const;
		void insertToTopology(const Segment& s, Vertex *mayC, Vertex *mayD);
		void removeFromTopology(const Segment &s, Vertex *mayC, Vertex *mayD);
		void addSupplyVerts(Vertex *a, Vertex *b, Vertex *c, Vertex *d, int mode);		
		void addEnforcedSegmentsAdjacency(Vertex *a, Vertex *b, Vertex *c, Vertex *d, int mode);
		void addEnforcedSegmentAdjacency(Vertex *a, Vertex *b, Vertex *c, Vertex *d);

		bool adjacent2VertexSlow(Vertex *w, Tetrahedron *t) const;
		bool adjacent2SegmentSlow(const Segment &s, Tetrahedron *t) const;

		static bool parityCheck(const Vertex *x, const Vertex *y);
		static bool edgeOrderCheck(const Vertex *a, const Vertex *b);
		static bool verticesOrderCheck(const Vertex *ori, const Vertex *end, const Vertex *c, const Vertex *d);

		std::vector<Vertex *> vertices;
		std::vector<Vertex *> enforcedEdgeVertices;
		MemoryPool<TetVertexListNode, std::alignment_of<TetVertexListNode>::value> *nodePool;
		MemoryPool<EdgeListNode, std::alignment_of<EdgeListNode>::value> *edgeNodePool;
		MemoryArena<Vertex> *vertPool;
		ThreadUnsafeFreelist<sizeof(uintptr_t), 4 * sizeof(uintptr_t)> *pointerList;
		VertexLabeler labeler;
		Vertex *ghost;
		Vertex *deadVerticesStack;
	};

	inline bool matchVertexFlag(VertexType type, VertexType flag) {
		return (type & flag) == flag;
	}

	inline bool matchOneOfVertexFlags(VertexType type, VertexType flag) {
		return (type & flag) != 0;
	}

	inline void Vertex::setEdgeList(EdgeListNode *n) {
		if (!matchOneOfVertexFlags(type, Vertex_LowDimension))
			pointers = uintptr_t(n) | (pointers & 0x1);
		else {
			uintptr_t *p = (uintptr_t *)(pointers & (~0x1));
			p[0] = uintptr_t(n);
		}
	}

	inline EdgeListNode *Vertex::getEdgeList() {
		if (!matchOneOfVertexFlags(type, Vertex_LowDimension))
			return (EdgeListNode *)(pointers & (~0x1));
		else {
			uintptr_t *p = (uintptr_t *)(pointers & (~0x1));
			return (EdgeListNode *)p[0];
		}
	}

	inline void Vertex::setSupplyVertexPointer(Vertex *v) {
		if (!matchOneOfVertexFlags(type, Vertex_LowDimension))
			pointers = (uintptr_t(v) | 0x2) | (pointers & 0x1);
		else {
			uintptr_t *p = (uintptr_t *)(pointers & (~0x1));
			p[0] = uintptr_t(v) | 0x2;
		}
	}

	inline Vertex *Vertex::getSupplyVertex() const {
		if (!matchOneOfVertexFlags(type, Vertex_LowDimension))
			return (Vertex *)(pointers & (~0x3));
		else {
			uintptr_t *p = (uintptr_t *)(pointers & (~0x1));
			return (Vertex *)(p[0] & (~0x2));
		}
	}

	inline bool Vertex::hasEdgeList() const {
		if (!matchOneOfVertexFlags(type, Vertex_LowDimension))
			return (pointers & (~0x1)) != 0 && (pointers & 0x2) == 0;
		else {
			uintptr_t *p = (uintptr_t *)(pointers & (~0x1));
			return p[0] != 0 && (p[0] & 0x2) == 0;
		}
	}

	inline void Vertex::setFaceLink(TriVertexListNode *link) {
		Assert(matchVertexFlag(type, VertexType::Vertex_Facet));
		if (matchVertexFlag(type, VertexType::Vertex_Volume)) {
			uintptr_t *p = (uintptr_t *)(pointers & (~0x1));
			p[1] = uintptr_t(link);
		}
		else
			pointers = uintptr_t(link) | (pointers & 0x1);
	}

	inline TriVertexListNode *Vertex::getFaceLink() const {
		Assert(matchVertexFlag(type, VertexType::Vertex_Facet));
		if (matchVertexFlag(type, VertexType::Vertex_Volume)) {
			uintptr_t *p = (uintptr_t *)(pointers & (~0x1));
			return (TriVertexListNode *)p[1];
		}
		else
			return (TriVertexListNode *)(pointers & (~0x1));
	}

	inline void Vertex::setEnforcedEdgeList(EdgeListNode *node) {
		Assert(matchVertexFlag(type, VertexType::Vertex_Segment));
		uintptr_t *p = (uintptr_t *)(pointers & (~0x1));
		int pos = matchVertexFlag(type, VertexType::Vertex_Facet) ? 2 : 1;
		p[pos] = uintptr_t(node) | (p[pos] & 0x1);
	}

	inline EdgeListNode *Vertex::getEnforcedEdgeList() const {
		Assert(matchVertexFlag(type, VertexType::Vertex_Segment));
		uintptr_t *p = (uintptr_t *)(pointers & (~0x1));
		int pos = matchVertexFlag(type, VertexType::Vertex_Facet) ? 2 : 1;
		return (EdgeListNode *)(p[pos] & (~0x1));
	}

	inline void Vertex::setEnforcedEdgeMark() {
		Assert(matchVertexFlag(type, VertexType::Vertex_Segment));
		uintptr_t *p = (uintptr_t *)(pointers & (~0x1));
		int pos = matchVertexFlag(type, VertexType::Vertex_Facet) ? 2 : 1;
		p[pos] |= 0x1;
	}

	inline bool Vertex::isEnforcedEdgeMarked() const {
		Assert(matchVertexFlag(type, VertexType::Vertex_Segment));
		uintptr_t *p = (uintptr_t *)(pointers & (~0x1));
		int pos = matchVertexFlag(type, VertexType::Vertex_Facet) ? 2 : 1;
		return (p[pos] & 0x1) == 0x1;
	}

	inline void Vertex::setOriSegmentIndex(int index) {
		Assert(matchVertexFlag(type, VertexType::Vertex_FreeSegment));
		int pos = matchVertexFlag(type, VertexType::Vertex_Facet) ? 3 : 2;
		uintptr_t *p = (uintptr_t *)(pointers & (~0x1));
		p[pos] = index;
	}

	inline int Vertex::getOriSegmentIndex() const {
		Assert(matchVertexFlag(type, VertexType::Vertex_FreeSegment));
		int pos = matchVertexFlag(type, VertexType::Vertex_Facet) ? 3 : 2;
		uintptr_t *p = (uintptr_t *)(pointers & (~0x1));
		return p[pos];
	}

	inline void Vertex::setPointers(uintptr_t *p) {
		pointers = uintptr_t(p);
	}

	inline uintptr_t *Vertex::getPointers() const {
		return (uintptr_t *)pointers;
	}

	inline Vertex *TriMeshDataStructure::getGhostVertex() {
		if (ghost == NULL) {
			if (vertPool == NULL) vertPool = new MemoryArena<Vertex>(128);
			ghost = vertPool->Alloc();
			ghost->type = VertexType::Vertex_Facet;
			ghost->setGhost();
			vertices.push_back(ghost);
		}
		return ghost;
	}

	inline Vertex* TetMeshDataStructure::getGhostVertex() {
		if (ghost == NULL) {
			ghost = vertPool->Alloc();
			ghost->type = VertexType::Vertex_Volume;
			ghost->setGhost();
		}
		return ghost;
	}

	inline bool TetMeshDataStructure::parityCheck(const Vertex *x, const Vertex *y) {
		constexpr int odd_mark = 1;
		return (x->getLabel() & odd_mark) == (y->getLabel() & odd_mark);
	}

	inline bool TetMeshDataStructure::edgeOrderCheck(const Vertex *a, const Vertex *b) {
		constexpr int mask = 2;
		int aLabel = a->getLabel(), bLabel = b->getLabel();
		return b->isGhost() || (!a->isGhost() && ((aLabel < bLabel) ^ ((aLabel & mask) == (bLabel & mask))));
	}

	inline bool TriMeshDataStructure::Contain(Vertex *v) const {
		Triangle f;
		return adjacent2Vertex(v, &f);
	}

	inline bool TriMeshDataStructure::Contain(const Segment &s) const {
		Vertex *w = NULL;
		return Adjacent(s, &w);
	}

	inline bool TriMeshDataStructure::Contain(const Triangle &f) const {
		if (matchVertexFlag(f.v[2]->getVertexType(), VertexType::Vertex_Facet)) {
			Vertex *w = NULL;
			if (Adjacent(Segment(f.v[0], f.v[1]), &w)) 
				return w == f.v[2];
		}

		return false;
	}

	inline bool TriMeshDataStructure::verticesOrderCheck(const Vertex *a, const Vertex *b, const Vertex *c) {
		return a->getLabel() < b->getLabel() && a->getLabel() < c->getLabel();
	}

	inline bool TetMeshDataStructure::Contain(Vertex *v) const {
		Tetrahedron t;
		return adjacent2Vertex(v, &t);
	}

	inline bool TetMeshDataStructure::Contain(const Triangle &f) const {
		Vertex *x = NULL;
		return Adjacent(f, &x);
	}

	inline bool TetMeshDataStructure::Contain(const Tetrahedron& t) const {
		Vertex *x = NULL;
		Adjacent(Triangle(t.v[1], t.v[2], t.v[3]), &x);
		return x == t.v[0];
	}

	inline bool TetMeshDataStructure::Contain(const Segment &s) const {
		Tetrahedron t;
		if (matchVertexFlag(s.v[0]->getVertexType(), VertexType::Vertex_Segment) && 
			matchVertexFlag(s.v[1]->getVertexType(), VertexType::Vertex_Segment)) 
			return adjacent2SegmentFast(s, &t);

		return false;
	}

	inline bool TetMeshDataStructure::fastVertexQueryCheck(Vertex *u) const {
		return u->hasEdgeList() || u->getSupplyVertex();
	}

	inline bool TetMeshDataStructure::fastSegmentQueryCheck(const Segment &s) const {
		return parityCheck(s.v[0], s.v[1]);
	}

	inline bool TetMeshDataStructure::isSegment(const Segment &s) const {
		Vertex *a = s.v[0], *b = s.v[1];
		if (!matchVertexFlag(a->type, VertexType::Vertex_Segment) || !matchVertexFlag(b->type, VertexType::Vertex_Segment))
			return false;
		if (!edgeOrderCheck(a, b)) std::swap(a, b);
		EdgeListNode *node = a->getEnforcedEdgeList();
		while (node != NULL && node->getEndVertex() != b)
			node = node->getNextNode();

		return node != NULL;
	}

	inline bool TetMeshDataStructure::adjacent2Segment(const Segment &s, Tetrahedron *t) const {
		if (!adjacent2SegmentFast(s, t)) return adjacent2SegmentSlow(s, t);
		return true;
	}

	inline bool TetMeshDataStructure::verticesOrderCheck(const Vertex *ori, const Vertex *end, const Vertex *c, const Vertex *d) {
		return (!parityCheck(c, ori) || (ori->getLabel() < c->getLabel() && end->getLabel() < c->getLabel())) &&
			(!parityCheck(d, ori) || (ori->getLabel() < d->getLabel() && end->getLabel() < d->getLabel())) &&
			(!parityCheck(c, d) || (ori->getLabel() < c->getLabel() && ori->getLabel() < d->getLabel()) ||
			(end->getLabel() < c->getLabel() && end->getLabel() < d->getLabel()));
	}

	inline TriMeshDataStructure::TriMeshConstIterator& TriMeshDataStructure::TriMeshConstIterator::operator++() {
		findNext();
		return *this;
	}

	inline TriMeshDataStructure::TriMeshConstIterator TriMeshDataStructure::TriMeshConstIterator::operator++(int) {
		TriMeshConstIterator temp = *this;
		findNext();
		return temp;
	}

	inline bool TriMeshDataStructure::TriMeshConstIterator::operator==(const TriMeshConstIterator& right) const {
		return vertIter == right.vertIter && parent == right.parent && child == right.child;
	}

	inline bool TriMeshDataStructure::TriMeshConstIterator::operator!=(const TriMeshConstIterator& right) const {
		return !(*this == right);
	}

	inline TriMeshDataStructure::TriMeshConstCirculator& TriMeshDataStructure::TriMeshConstCirculator::operator++() {
		findNext();
		return *this;
	}

	inline TriMeshDataStructure::TriMeshConstCirculator TriMeshDataStructure::TriMeshConstCirculator::operator++(int) {
		TriMeshConstCirculator tmp = *this;
		findNext();
		return tmp;
	}

	inline bool TriMeshDataStructure::TriMeshConstCirculator::operator==(const TriMeshConstCirculator& right) const {
		return vert == right.vert && link == right.link;
	}

	inline bool TriMeshDataStructure::TriMeshConstCirculator::operator!=(const TriMeshConstCirculator& right) const {
		return !(*this == right);
	}

	inline bool TriMeshDataStructure::TriMeshConstCirculator::operator==(std::nullptr_t p) const {
		return link == NULL;
	}

	inline bool TriMeshDataStructure::TriMeshConstCirculator::operator!=(std::nullptr_t p) const {
		return !(*this == p);
	}

	inline TetMeshDataStructure::TetMeshConstIterator& TetMeshDataStructure::TetMeshConstIterator::operator++() {
		findNext();
		return *this;
	}

	inline TetMeshDataStructure::TetMeshConstIterator TetMeshDataStructure::TetMeshConstIterator::operator++(int) {
		TetMeshConstIterator temp = *this;
		findNext();
		return temp;
	}

	inline bool TetMeshDataStructure::TetMeshConstIterator::operator==(const TetMeshConstIterator& right) const {
		return vertIter == right.vertIter && linkHead == right.linkHead && parent == right.parent && child == right.child;
	}

	inline bool TetMeshDataStructure::TetMeshConstIterator::operator!=(const TetMeshConstIterator& right) const {
		return !(*this == right);
	}

	inline TetMeshDataStructure::SegmentConstIterator& TetMeshDataStructure::SegmentConstIterator::operator++() {
		findNext();
		return *this;
	}

	inline TetMeshDataStructure::SegmentConstIterator TetMeshDataStructure::SegmentConstIterator::operator++(int) {
		SegmentConstIterator temp = *this;
		findNext();
		return temp;
	}

	inline bool TetMeshDataStructure::SegmentConstIterator::operator==(const SegmentConstIterator& right) const {
		return edgeVertIter == right.edgeVertIter && edgeNode == right.edgeNode;
	}

	inline bool TetMeshDataStructure::SegmentConstIterator::operator!=(const SegmentConstIterator& right) const {
		return !(*this == right);
	}

	inline TetMeshDataStructure::TetMeshFacetConstCirculator& TetMeshDataStructure::TetMeshFacetConstCirculator::operator++() {
		findNext();
		return *this;
	}

	inline TetMeshDataStructure::TetMeshFacetConstCirculator TetMeshDataStructure::TetMeshFacetConstCirculator::operator++(int) {
		TetMeshFacetConstCirculator tmp = *this;
		findNext();
		return tmp;
	}

	inline bool TetMeshDataStructure::TetMeshFacetConstCirculator::operator==(const TetMeshFacetConstCirculator& right) const {
		if (ori == right.ori && dest == right.dest) {
			if (linkListHead && right.linkListHead)
				return listNode == right.listNode;
			else if (initAdjacentVert && right.initAdjacentVert)
				return movingVert == right.movingVert;
		}
		return false;
	}

	inline bool TetMeshDataStructure::TetMeshFacetConstCirculator::operator!=(const TetMeshFacetConstCirculator& right) const {
		return !(*this == right);
	}

	inline bool TetMeshDataStructure::TetMeshFacetConstCirculator::operator==(std::nullptr_t p) const {
		return linkListHead == NULL && initAdjacentVert == NULL;
	}

	inline bool TetMeshDataStructure::TetMeshFacetConstCirculator::operator!=(std::nullptr_t p) const {
		return !(*this == p);
	}
}

#endif