#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_MESHER_DELPRIMITIVE_H
#define ODER_MESHER_DELPRIMITIVE_H

#include "oder.h"
#include "latool.h"
#include <unordered_map>
#include <set>

namespace ODER {
#define REAL double

#define NEXT_F(i) (((i)+1)%3)
#define NEXT_T(i) (((i)+1)%4)
#define SQRTF_2 REAL(1.41422)
#define ODD_MAKE 1

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

	class EdgeListNode;

	struct Vertex {
		Vertex() : weight(0.0), label(-1), pointer(0) {}
		template<class FT> Vertex(const VectorBase<FT>& vv) : vert{ vv.x, vv.y, vv.z }, weight(0), label(-1), pointer(0){}
		Vertex(const DelVector& vv, REAL w = 0) :vert(vv), weight(w), label(-1), pointer(0) {}
		void setGhost() {
			vert.x = FLT_MAX; vert.y = FLT_MAX; vert.z = FLT_MAX;
			weight = -FLT_MAX;
			label = labeler.getSpecilGhostLabel();
		}
		bool isGhost() {
			return label == labeler.getSpecilGhostLabel();
		}
		void setLabel() {
			label = labeler.getLabel();
		}
		int getLabel() const {
			return label;
		}
		void setListPointer(EdgeListNode *n) {
			pointer = uintptr_t(n);
		}
		void setEndVertexPoint(Vertex *v) {
			pointer = uintptr_t(v) | 0x1;
		}
		EdgeListNode *getListHead() {
			return (EdgeListNode *)pointer;
		}
		Vertex *getEndVertex() const {
			return (Vertex *)(pointer & (~0x1));
		}
		bool hasList() const {
			return pointer != 0 && (pointer & 0x1) == 0;
		}
		bool operator<(const Vertex& v) const {
			return label < v.label;
		}
		bool operator>(const Vertex& v) const {
			return label > v.label;
		}
		DelVector vert;
		REAL weight;
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
			return (Vertex *)(currentVert & (~0x1));
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
		void setRationAndRadius(const Vertex *ghost);
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
			int smallerLable = s.v[0]->getLabel();
			int biggerLable = s.v[1]->getLabel();
			return std::hash<int>()(((biggerLable * (biggerLable + 1)) >> 1) + smallerLable);
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
		bool Adjacent(const Segment &s, Vertex **w) const;
		bool adjacent2Vertex(Vertex *w, Face *f) const;
		bool Contain(const Segment &s) const;
		bool Contain(const Face &f) const;
		bool findIntersectedFace(Vertex *a, const DelVector& bb, const DelVector& above, Face *f) const;
		void Clear();
		std::set<Face, face_compare> getTriangleSet(bool ghost) const;
		std::vector<Face> getTriangleVector(bool ghost) const;
		~TriMeshDataStructure();
	private:
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
		bool Adjacent(const Face &s, Vertex **w) const;
		bool adjacent2Vertex(Vertex *w, Tetrahedron *t) const;
		std::set<Tetrahedron, tet_compare> getTetraheronSet(bool ghost) const;
		std::vector<Tetrahedron> getTetraheronVector(bool ghost) const;
		~TetMeshDataStructure();
	private:
		bool Adjacent(Vertex* w, Vertex *x, Vertex *y, Vertex **z) const;
		void insertToTopology(const Segment& s, Vertex *mayC, Vertex *mayD);
		void removeFromTopology(const Segment &s, Vertex *mayC, Vertex *mayD);
		void addSupplyVerts(Vertex *a, Vertex *b, Vertex *c, Vertex *d, int mode);

		std::vector<Vertex *> vertices;
		MemoryPool<VertexListNode> *nodePool;
		MemoryPool<EdgeListNode> *edgeNodePool;
	};


	inline bool parityCheck(const Vertex *x, const Vertex *y) {
		return (std::hash<int>()(x->getLabel()) & ODD_MAKE) == (std::hash<int>()(y->getLabel()) & ODD_MAKE);
	}

	inline bool TriMeshDataStructure::Contain(const Segment &s) const {
		Vertex *w;
		return Adjacent(s, &w);
	}

	inline bool TriMeshDataStructure::Contain(const Face &f) const {
		Vertex *w;
		bool find = Adjacent(Segment(f.v[0], f.v[1]), &w);
		if (find) return w == f.v[2];

		return false;
	}
}

#endif