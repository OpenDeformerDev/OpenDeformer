#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_MESHER_DELPRIMITIVE_H
#define ODER_MESHER_DELPRIMITIVE_H

#include "oder.h"
#include "latool.h"
#include <unordered_map>
#include <set>
#include <numeric>

#define REAL double

#define NEXT_F(i) (((i)+1)%3)
#define NEXT_T(i) (((i)+1)%4)
#define SQRTF_2 REAL(1.41422)

namespace ODER {
	using DelVector = VectorBase<REAL>;

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

		Vertex_FixedSegment = Vertex_Fixed | Vertex_Segment,
		Vertex_FixedFacet = Vertex_Fixed | Vertex_Facet,
		Vertex_FixedVolume = Vertex_Fixed | Vertex_Volume,
		
		Vertex_FreeSegment = Vertex_Free | Vertex_Segment,
		Vertex_FreeFacet = Vertex_Free | Vertex_Facet,
		Vertex_FreeVolume = Vertex_Free | Vertex_Facet,
	};

	class EdgeListNode;

	struct Vertex {
		Vertex() : weight(0), relaxedInsetionRadius(std::numeric_limits<REAL>::max()), label(-1), vertexPointer(0), type(VertexType::Vertex_Undefined){}
		template<class FT> explicit Vertex(const VectorBase<FT>& vv, VertexType t = VertexType::Vertex_Undefined)
			: vert{ vv.x, vv.y, vv.z }, weight(0), relaxedInsetionRadius(std::numeric_limits<REAL>::max()), label(-1), vertexPointer(0), type(t) {}
		explicit Vertex(const DelVector& vv, REAL w = 0, VertexType t = VertexType::Vertex_Undefined)
			: vert(vv), weight(w), label(-1), vertexPointer(0), type(t) {}
		void setGhost() {
			constexpr REAL max = std::numeric_limits<REAL>::max();
			vert.x = max; vert.y = max; vert.z = max;
			weight = -max;
			label = VertexLabeler::getSpecilGhostLabel();
			vertexPointer = 0;
			relaxedInsetionRadius = max;
			type = VertexType::Vertex_Undefined;
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
		void setVertexPointer(Vertex *v) {
			vertexPointer = uintptr_t(v) | (vertexPointer & 0x1);
		}
		Vertex *getPointedVertex() const {
			return (Vertex *)(vertexPointer & (~0x1));
		}
		void setMark() {
			vertexPointer |= 0x1;
		}
		void unSetMark() {
			vertexPointer &= (~0x1);
		}
		bool isMarked() const {
			return (vertexPointer & 0x1) == 0x1;
		}
		DelVector vert;
		REAL weight;
		REAL relaxedInsetionRadius;
		VertexType type;
	private:
		int label;
		uintptr_t vertexPointer;//point to a vertex which is one of the vertices of latest inserted vertex
	};

	class TriVertexListNode {
	public:
		TriVertexListNode() :currentVert(0x1), next(NULL), faceIndex(-1) {}
		void setVertex(Vertex *vert) {
			currentVert = uintptr_t(vert);
		}
		Vertex *getVertex() const {
			return (Vertex *)(currentVert & (~0x3));
		}
		TriVertexListNode *getNextNode() const {
			return next;
		}
		void setNextNode(TriVertexListNode *nextNode) {
			next = nextNode;
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
			currentVert |= 0x2;
		}
		void unSetMark() {
			currentVert &= (~0x2);
		}
		void setIndex(int index) { 
			faceIndex = index;
		}
		int getIndex() const { 
			return faceIndex; 
		}
		bool isMarked() {
			return (currentVert & 0x2) == 0x2;
		}
	private:
		uintptr_t currentVert;
		TriVertexListNode *next;
		int faceIndex;
	};


	class TetVertexListNode {
	public:
		TetVertexListNode() :currentVert(0x1), next(NULL) {}
		void setVertex(Vertex *vert) {
			currentVert = uintptr_t(vert);
		}
		Vertex *getVertex() const {
			return (Vertex *)(currentVert & (~0x3));
		}
		TetVertexListNode *getNextNode() const {
			return next;
		}
		void setNextNode(TetVertexListNode *nextNode) {
			next = nextNode;
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
			currentVert |= 0x2;
		}
		void unSetMark() {
			currentVert &= (~0x2);
		}
		bool isMarked() {
			return (currentVert & 0x2) == 0x2;
		}
	private:
		uintptr_t currentVert;
		TetVertexListNode *next;
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
			link = Link;
		}
		TetVertexListNode* getLink() const {
			return link;
		}
		EdgeListNode *getNextNode() const {
			return next;
		}
		void setNextNode(EdgeListNode *nextNode) {
			next = nextNode;
		}
	private:
		Vertex *endVertex;
		TetVertexListNode *link;
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
		Vertex *v[2];
	};
	struct Face {
		Face() {
			v[0] = v[1] = v[2] = NULL;
			index = -1;
		}
		Face(Vertex *v0, Vertex *v1, Vertex *v2, bool ordered = false) {
			initVertices(v0, v1, v2, ordered);
			index = -1;
		}
		Face(Vertex *v0, Vertex *v1, Vertex *v2, int index, bool ordered = false) {
			initVertices(v0, v1, v2, ordered);
			this->index = index;
		}
		void sortVertices();
		bool operator==(const Face& f) const {
			return v[0] == f.v[0] && v[1] == f.v[1] && v[2] == f.v[2];
		}
		bool operator!=(const Face &f) const {
			return v[0] != f.v[0] || v[1] != f.v[1] || v[2] != f.v[2];
		}
		Vertex *v[3];
		int index;
	private:
		void initVertices(Vertex *v0, Vertex *v1, Vertex *v2, bool ordered);
	};

	struct Tetrahedron {
		Tetrahedron() {
			v[0] = v[1] = v[2] = v[3] = NULL;
			reRation = r = 0.f;
		}
		Tetrahedron(Vertex *v0, Vertex *v1, Vertex *v2, Vertex *v3, bool ordered = false);
		void setRationAndRadius();
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
		inline REAL getREration() const { return reRation; }
		inline REAL getRadius() const { return r; }
		Vertex *v[4];
	private:
		REAL reRation;
		REAL r;
	};

	struct face_compare {
		bool operator()(const Face &left, const Face &right) const {
			if (left.v[2] == right.v[2]) {
				if (left.v[1] == right.v[1])
					return left.v[0] < right.v[0];
				return left.v[1] < right.v[1];
			}
			return left.v[2] < right.v[2];
		}
	};

	struct face_hash {
		size_t operator()(const Face& f) const {
			size_t hashVal = 0;
			hashCombine(hashVal, f.v[0]->getLabel(), f.v[1]->getLabel(), f.v[2]->getLabel());
			return hashVal;
		}
	};

	struct segment_ordered_hash {
		size_t operator()(const Segment &s) const {
			int64_t smallerLable = s.v[0]->getLabel();
			int64_t biggerLable = s.v[1]->getLabel();
			return std::hash<int64_t>()(((biggerLable * (biggerLable + 1)) >> 1) + smallerLable);
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
			return std::hash<int>()(v->getLabel());
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
		bool adjacent2Vertex(Vertex *w, Face *f) const;
		bool Contain(Vertex *v) const;
		bool Contain(const Segment &s) const;
		bool Contain(const Face &f) const;
		bool findIntersectedFace(Vertex *a, const DelVector& bb, Face *f) const;
		void Clear();
		void Reserve(size_t n) { topology.reserve(n); }
		std::vector<Face> getTriangles(bool ghost) const;
		void getTriangles(bool ghost, std::vector<Face>& triangles) const;
		~TriMeshDataStructure();

		class TriMeshConstIterator {
		public:
			using iterator_category = std::input_iterator_tag;
			using value_type = Face;
			using difference_type = ptrdiff_t;
			using size_type = size_t;
			using reference = const Face&;
			using pointer = const Face *;

			TriMeshConstIterator() : parent(NULL), child(NULL) {}
			TriMeshConstIterator(const std::unordered_map<Vertex *, TriVertexListNode *, vertex_hash>::const_iterator& iter,
				const std::unordered_map<Vertex *, TriVertexListNode *, vertex_hash>::const_iterator& end);
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

			Face current;
			TriVertexListNode *parent, *child;
			std::unordered_map<Vertex *, TriVertexListNode *, vertex_hash>::const_iterator topologyIter;
			std::unordered_map<Vertex *, TriVertexListNode *, vertex_hash>::const_iterator topologyEnd;
		};

		using const_iterator = TriMeshConstIterator;
		const_iterator begin() const { return const_iterator(topology.begin(), topology.end()); }
		const_iterator end() const { return const_iterator(topology.end(), topology.end()); }

	private:
		bool getAdjacentListNode(Vertex* u, Vertex* v, TriVertexListNode **w) const;
		void insertToTopology(Vertex *a, Vertex *b, Vertex *c, int index);
		void removeFromTopology(Vertex *a, Vertex *b, Vertex *c);
		static bool verticesOrderCheck(const Vertex *a, const Vertex *b, const Vertex *c);

		std::unordered_map<Vertex *, TriVertexListNode *, vertex_hash> topology;
		MemoryPool<TriVertexListNode> *nodePool;
	};

	class TetMeshDataStructure {
	public:
		TetMeshDataStructure();
		void addTetrahedron(Vertex *a, Vertex *b, Vertex *c, Vertex *d);
		void deleteTetrahedron(Vertex *a, Vertex *b, Vertex *c, Vertex *d);
		void setMark(Vertex *u, Vertex *v, Vertex *w);
		void unSetMark(Vertex *u, Vertex *v, Vertex *w);
		bool isMarked(Vertex *u, Vertex *v, Vertex *w) const;
		bool Adjacent(const Face &f, Vertex **w) const;
		bool adjacent2Vertex(Vertex *w, Tetrahedron *t) const;
		std::vector<Tetrahedron> getTetrahedrons(bool ghost) const;
		void getTetrahedrons(bool ghost, std::vector<Tetrahedron>& tets) const;
		bool Contain(Vertex *v) const;
		bool Contain(const Face &f) const;
		bool Contain(const Tetrahedron& t) const;
		void Clear();
		void Reserve(size_t n) { topology.reserve(n); }
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
			TetMeshConstIterator(const std::unordered_map<Vertex *, EdgeListNode *, vertex_hash>::const_iterator& iter,
				const std::unordered_map<Vertex *, EdgeListNode *, vertex_hash>::const_iterator& end);
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

			std::unordered_map<Vertex *, EdgeListNode *, vertex_hash>::const_iterator topologyIter;
			std::unordered_map<Vertex *, EdgeListNode *, vertex_hash>::const_iterator topologyEnd;
		};

		using const_iterator = TetMeshConstIterator;
		const_iterator begin() const { return const_iterator(topology.begin(), topology.end()); }
		const_iterator end() const { return const_iterator(topology.end(), topology.end()); }

	private:
		bool getAdjacentListNode(const Face& f, TetVertexListNode **z) const;
		bool getAdjacentListNode(Vertex* w, Vertex *x, Vertex *y, TetVertexListNode **z) const;
		void insertToTopology(const Segment& s, Vertex *mayC, Vertex *mayD);
		void removeFromTopology(const Segment &s, Vertex *mayC, Vertex *mayD);
		void addSupplyVerts(Vertex *a, Vertex *b, Vertex *c, Vertex *d, int mode);
		static bool parityCheck(const Vertex *x, const Vertex *y);
		static bool edgeOrderCheck(const Vertex *a, const Vertex *b);
		static bool verticesOrderCheck(const Vertex *ori, const Vertex *end, const Vertex *c, const Vertex *d);

		std::unordered_map<Vertex *, EdgeListNode *, vertex_hash> topology;
		MemoryPool<TetVertexListNode> *nodePool;
		MemoryPool<EdgeListNode> *edgeNodePool;
	};


	inline bool TetMeshDataStructure::parityCheck(const Vertex *x, const Vertex *y) {
		constexpr int odd_mark = 1;
		return (std::hash<int>()(x->getLabel()) & odd_mark) == (std::hash<int>()(y->getLabel()) & odd_mark);
	}

	inline bool TetMeshDataStructure::edgeOrderCheck(const Vertex *a, const Vertex *b) {
		constexpr int mask = 2;
		int aLabel = a->getLabel(), bLabel = b->getLabel();
		return b->isGhost() || (!a->isGhost() && ((aLabel < bLabel) ^ ((std::hash<int>()(aLabel) & mask) == (std::hash<int>()(bLabel) & mask))));
	}

	inline bool TriMeshDataStructure::Contain(Vertex *v) const {
		Face f;
		return adjacent2Vertex(v, &f);
	}

	inline bool TriMeshDataStructure::Contain(const Segment &s) const {
		Vertex *w = NULL;
		return Adjacent(s, &w);
	}

	inline bool TriMeshDataStructure::Contain(const Face &f) const {
		Vertex *w = NULL;
		bool find = Adjacent(Segment(f.v[0], f.v[1]), &w);
		if (find) return w == f.v[2];

		return false;
	}

	inline bool TriMeshDataStructure::verticesOrderCheck(const Vertex *a, const Vertex *b, const Vertex *c) {
		return a->getLabel() < b->getLabel() && a->getLabel() < c->getLabel();
	}

	inline bool TetMeshDataStructure::Contain(Vertex *v) const {
		Tetrahedron t;
		return adjacent2Vertex(v, &t);
	}

	inline bool TetMeshDataStructure::Contain(const Face &f) const {
		Vertex *x = NULL;
		return Adjacent(f, &x);
	}

	inline bool TetMeshDataStructure::Contain(const Tetrahedron& t) const {
		Vertex *x = NULL;
		Adjacent(Face(t.v[1], t.v[2], t.v[3]), &x);
		return x == t.v[0];
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
		return topologyIter == right.topologyIter && parent == right.parent && child == right.child;
	}

	inline bool TriMeshDataStructure::TriMeshConstIterator::operator!=(const TriMeshConstIterator& right) const {
		return !(*this == right);
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
		return topologyIter == right.topologyIter && linkHead == right.linkHead && parent == right.parent && child == right.child;
	}

	inline bool TetMeshDataStructure::TetMeshConstIterator::operator!=(const TetMeshConstIterator& right) const {
		return !(*this == right);
	}

	inline bool matchVertexFlag(VertexType type, VertexType flag) {
		return (type & flag) == flag;
	}
}

#endif
