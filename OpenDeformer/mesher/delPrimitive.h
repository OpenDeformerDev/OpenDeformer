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

	class Labeler {
	public:
		Labeler() : label(1) {}
		int getLabel() {
			return label++;
		}
		int getSpecilGhostLabel() {
			return 0;
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

		Vertex_FixedSegment = Vertex_Fixed | Vertex_Segment,
		Vertex_FixedFacet = Vertex_Fixed | Vertex_Facet,
		Vertex_FixedVolume = Vertex_Fixed | Vertex_Volume,
		
		Vertex_FreeSegment = Vertex_Free | Vertex_Segment,
		Vertex_FreeFacet = Vertex_Free | Vertex_Facet,
		Vertex_FreeVolume = Vertex_Free | Vertex_Facet,

		Vertex_Low_Dimen = Vertex_Segment | Vertex_Facet
	};

	class EdgeListNode;

	struct Vertex {
		Vertex() : weight(0.0), label(-1), pointer(0), type(VertexType::Vertex_Undefined){}
		template<class FT> explicit Vertex(const VectorBase<FT>& vv, VertexType t = VertexType::Vertex_Undefined)
			: vert{ vv.x, vv.y, vv.z }, weight(0), label(-1), pointer(0), type(t) {}
		explicit Vertex(const DelVector& vv, REAL w = 0, VertexType t = VertexType::Vertex_Undefined)
			: vert(vv), weight(w), label(-1), pointer(0), type(t) {}
		void setGhost() {
			constexpr REAL max = std::numeric_limits<REAL>::max();
			vert.x = max; vert.y = max; vert.z = max;
			weight = -max;
			label = labeler.getSpecilGhostLabel();
			type = VertexType::Vertex_Undefined;
		}
		bool isGhost() const {
			return label == labeler.getSpecilGhostLabel();
		}
		void setLabel() {
			label = labeler.getLabel();
		}
		int getLabel() const {
			return label;
		}
		void setListPointer(EdgeListNode *n) {
			pointer = uintptr_t(n) | (pointer & 0x2);
		}
		void setVertexPointer(Vertex *v) {
			pointer = (uintptr_t(v) | 0x1) | (pointer & 0x2);
		}
		EdgeListNode *getListHead() {
			return (EdgeListNode *)(pointer & (~0x2));
		}
		Vertex *getPointedVertex() const {
			return (Vertex *)(pointer & (~0x3));
		}
		void setMark() {
			pointer |= 0x2;
		}
		void unSetMark() {
			pointer &= (~0x2);
		}
		uintptr_t getRawPointer() const {
			return pointer;
		}
		void setRawPointer(uintptr_t p) {
			pointer = p;
		}
		bool isMarked() const {
			return (pointer & 0x2) == 0x2;
		}
		bool hasList() const {
			return (pointer & (~0x2)) != 0 && (pointer & 0x1) == 0;
		}
		bool operator<(const Vertex& v) const {
			return label < v.label;
		}
		bool operator>(const Vertex& v) const {
			return label > v.label;
		}
		DelVector vert;
		REAL weight;
		VertexType type;
	private:
		int label;
		uintptr_t pointer;//point to the adjaceny list or another vertex
		static Labeler labeler;
	};


	class VertexListNode {
	public:
		VertexListNode() :currentVert(0x1), next(NULL) {}
		void setVertex(Vertex *vert) {
			currentVert = uintptr_t(vert);
		}
		Vertex *getVertex() const {
			return (Vertex *)(currentVert & (~0x3));
		}
		VertexListNode *getNextNode() const {
			return next;
		}
		void setNextNode(VertexListNode *nextNode) {
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
		VertexListNode *next;
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
		void setLink(VertexListNode *Link) {
			link = Link;
		}
		VertexListNode* getLink() const {
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
		VertexListNode *link;
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
		}
		Face(Vertex *v0, Vertex *v1, Vertex *v2, bool ordered = false) {
			if (ordered) {
				int min = v0->getLabel() < v1->getLabel() ? (v0->getLabel() < v2->getLabel() ? 0 : 2) : (v1->getLabel() < v2->getLabel() ? 1 : 2);
				switch (min) {
				case 0:
					v[0] = v0; v[1] = v1; v[2] = v2;
					break;
				case 1:
					v[0] = v1; v[1] = v2; v[2] = v0;
					break;
				case 2:
					v[0] = v2; v[1] = v0; v[2] = v1;
					break;
				}
			}
			else {
				v[0] = v0; v[1] = v1; v[2] = v2;
			}
		}
		void sortVertices() {
			int min = v[0]->getLabel() < v[1]->getLabel() ? (v[0]->getLabel() < v[2]->getLabel() ? 0 : 2) : (v[1]->getLabel() < v[2]->getLabel() ? 1 : 2);
			if (min == 1) {
				std::swap(v[0], v[1]);
				std::swap(v[1], v[2]);
			}
			else if (min == 2) {
				std::swap(v[0], v[1]);
				std::swap(v[0], v[2]);
			}
		}
		bool operator==(const Face& f) const {
			return v[0] == f.v[0] && v[1] == f.v[1] && v[2] == f.v[2];
		}
		bool operator!=(const Face &f) const {
			return v[0] != f.v[0] || v[1] != f.v[1] || v[2] != f.v[2];
		}
		Vertex *v[3];
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
		bool operator<(const Tetrahedron& t) const {
			if (reRation <= maxREration && t.reRation <= maxREration)
				return r < t.r;
			return reRation < t.reRation;
		}
		inline REAL getREration() const { return reRation; }
		inline REAL getRadius() const { return r; }
		Vertex *v[4];
		static REAL maxREration;
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
			hashCombine(hashVal, f.v[0]->getLabel());
			hashCombine(hashVal, f.v[1]->getLabel());
			hashCombine(hashVal, f.v[2]->getLabel());
			return hashVal;
		}
	};

	struct segment_hash {
		size_t operator()(const Segment &s) const {
			int64_t smallerLable = s.v[0]->getLabel();
			int64_t biggerLable = s.v[1]->getLabel();
			return std::hash<int64_t>()(((biggerLable * (biggerLable + 1)) >> 1) + smallerLable);
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
			hashCombine(hashVal, t.v[0]->getLabel());
			hashCombine(hashVal, t.v[1]->getLabel());
			hashCombine(hashVal, t.v[2]->getLabel());
			hashCombine(hashVal, t.v[3]->getLabel());
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
		void addTriangle(Vertex *a, Vertex *b, Vertex *c);
		void deleteTriangle(Vertex *a, Vertex *b, Vertex *c);
		void setDeletedMark(Vertex *u, Vertex *v);
		void unSetDeletedMark(Vertex *u, Vertex *v);
		void setDeletedMark(Vertex *a, Vertex *b, Vertex *c);
		void unSetDeletedMark(Vertex *a, Vertex *b, Vertex *c);
		void setMark(Vertex *u, Vertex *v);
		void unSetMark(Vertex *u, Vertex *v);
		bool isMarked(Vertex *u, Vertex *v) const;
		bool Adjacent(const Segment &s, Vertex **w) const;
		bool adjacent2Vertex(Vertex *w, Face *f) const;
		bool Contain(Vertex *v) const;
		bool Contain(const Segment &s) const;
		bool Contain(const Face &f) const;
		bool findIntersectedFace(Vertex *a, const DelVector& bb, Face *f) const;
		void Clear();
		std::vector<Face> getTriangles(bool ghost) const;
		void getTriangles(bool ghost, std::vector<Face>& triangles) const;
		~TriMeshDataStructure();
	private:
		bool getAdjacentListNode(Vertex* u, Vertex* v, VertexListNode **w) const;
		void insertToTopology(Vertex *a, Vertex *b, Vertex *c);
		void removeFromTopology(Vertex *a, Vertex *b, Vertex *c);

		std::unordered_map<Vertex *, VertexListNode *, vertex_hash> topology;
		MemoryPool<VertexListNode> *nodePool;
	};

	class TetMeshDataStructure {
	public:
		TetMeshDataStructure();
		void addTetrahedron(Vertex *a, Vertex *b, Vertex *c, Vertex *d);
		void deleteTetrahedron(Vertex *a, Vertex *b, Vertex *c, Vertex *d);
		void setMark(Vertex *u, Vertex *v, Vertex *w);
		void unSetMark(Vertex *u, Vertex *v, Vertex *w);
		bool testAndMark(Vertex *u, Vertex *v, Vertex *w);
		bool isMarked(Vertex *u, Vertex *v, Vertex *w) const;
		bool Adjacent(const Face &f, Vertex **w) const;
		bool adjacent2Vertex(Vertex *w, Tetrahedron *t) const;
		std::vector<Tetrahedron> getTetrahedrons(bool ghost) const;
		void getTetrahedrons(bool ghost, std::vector<Tetrahedron>& tets) const;
		bool Contain(Vertex *v) const;
		bool Contain(const Face &f) const;
		bool Contain(const Tetrahedron& t) const;
		void Clear();
		~TetMeshDataStructure();
	private:
		bool getAdjacentListNode(const Face& f, VertexListNode **z) const;
		bool getAdjacentListNode(Vertex* w, Vertex *x, Vertex *y, VertexListNode **z) const;
		void insertToTopology(const Segment& s, Vertex *mayC, Vertex *mayD);
		void removeFromTopology(const Segment &s, Vertex *mayC, Vertex *mayD);
		void addSupplyVerts(Vertex *a, Vertex *b, Vertex *c, Vertex *d, int mode);
		bool parityCheck(const Vertex *x, const Vertex *y) const;
		bool edgeOrderCheck(const Vertex *a, const Vertex *b) const;

		std::vector<Vertex *> vertices;
		MemoryPool<VertexListNode> *nodePool;
		MemoryPool<EdgeListNode> *edgeNodePool;
	};


	inline bool TetMeshDataStructure::parityCheck(const Vertex *x, const Vertex *y) const {
		constexpr int odd_mark = 1;
		return (std::hash<int>()(x->getLabel()) & odd_mark) == (std::hash<int>()(y->getLabel()) & odd_mark);
	}

	inline bool TetMeshDataStructure::edgeOrderCheck(const Vertex *a, const Vertex *b) const {
		constexpr int mask = 2;
		int aLabel = a->getLabel(), bLabel = b->getLabel();
		return b->isGhost() || ((aLabel < bLabel) ^ ((std::hash<int>()(aLabel) & mask) == (std::hash<int>()(bLabel) & mask)));
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

	inline bool matchVertexFlag(VertexType type, VertexType flag) {
		return (type & flag) == flag;
	}

}

#endif
