#include "stdafx.h"
#include "predicate.h"
#include "delPrimitive.h"
#include "memory.h"
#include <numeric>

namespace ODER {
	Labeler Vertex::labeler;

	REAL Tetrahedron::maxREration = std::numeric_limits<REAL>::max();

	Tetrahedron::Tetrahedron(Vertex *v0, Vertex *v1, Vertex *v2, Vertex *v3, bool ordered) {
		reRation = r = 0.0;
		v[0] = v0; v[1] = v1; v[2] = v2; v[3] = v3;
		if (ordered)
			sortVertices();
	}

	void Tetrahedron::setRationAndRadius() {
		if (v[0]->isGhost() || v[1]->isGhost() || v[2]->isGhost() || v[3]->isGhost()) {
			reRation = 0.0;
			r = 0.0;
			return;
		}

		DelVector da = v[0]->vert - v[3]->vert;
		DelVector db = v[1]->vert - v[3]->vert;
		DelVector dc = v[2]->vert - v[3]->vert;
		DelVector ca = v[0]->vert - v[2]->vert;
		DelVector cb = v[1]->vert - v[2]->vert;
		DelVector ba = v[0]->vert - v[1]->vert;

		const Predicator<REAL> predicator;
		predicator.Orthosphere(v[0]->vert, v[0]->weight, v[1]->vert, v[1]->weight, v[2]->vert, v[2]->weight, 
			v[3]->vert, v[3]->weight, NULL, &r);

		REAL minEdgeLength = sqrt(std::min(da.length2(), std::min(db.length2(), 
			std::min(dc.length2(), std::min(ca.length2(), std::min(cb.length2(), ba.length2()))))));
		reRation = r / minEdgeLength;
	}

	void Tetrahedron::sortVertices() {
		int min = 0, max = 0;
		for (int i = 1; i < 4; i++) {
			if (v[i]->getLabel() < v[min]->getLabel())
				min = i;
			else if (v[i]->getLabel() > v[max]->getLabel())
				max = i;
		}
		if (min == 0) {
			if (max == 1) {
				std::swap(v[1], v[2]);
				std::swap(v[2], v[3]);
			}
			else if (max == 2) {
				std::swap(v[2], v[3]);
				std::swap(v[1], v[2]);
			}
		}
		else if (min == 3) {
			std::swap(v[0], v[3]);
			if (max == 0) {
				std::swap(v[1], v[2]);
			}
			else {
				std::swap(v[max], v[3]);
			}
		}
		else {
			std::swap(v[0], v[min]);
			if (max == 3)
				std::swap(v[1], v[2]);
			else if (max == 0)
				std::swap(v[min], v[3]);
			else
				std::swap(v[max], v[3]);
		}
	}

	TriMeshDataStructure::TriMeshDataStructure() {
		topology.reserve(32);
		nodePool = new MemoryPool<VertexListNode>(64);
	}

	void TriMeshDataStructure::addTriangle(Vertex *a, Vertex *b, Vertex *c) {
		insertToTopology(a, b, c);
		insertToTopology(b, c, a);
		insertToTopology(c, a, b);
	}

	void TriMeshDataStructure::deleteTriangle(Vertex *a, Vertex *b, Vertex *c) {
		removeFromTopology(a, b, c);
		removeFromTopology(b, c, a);
		removeFromTopology(c, a, b);
	}

	void TriMeshDataStructure::setDeletedMark(Vertex *u, Vertex *v) {
		VertexListNode *node = NULL;
		if (getAdjacentListNode(u, v, &node))
			node->setDeletedMark();
	}

	void TriMeshDataStructure::unSetDeletedMark(Vertex *u, Vertex *v) {
		VertexListNode *node = NULL;
		if (getAdjacentListNode(u, v, &node))
			node->unSetDeletedMark();
	}

	void TriMeshDataStructure::setDeletedMark(Vertex *a, Vertex *b, Vertex *c) {
		setDeletedMark(a, b);
		setDeletedMark(b, c);
		setDeletedMark(c, a);
	}

	void TriMeshDataStructure::unSetDeletedMark(Vertex *a, Vertex *b, Vertex *c) {
		unSetDeletedMark(a, b);
		unSetDeletedMark(b, c);
		unSetDeletedMark(c, a);
	}

	void TriMeshDataStructure::setMark(Vertex *u, Vertex *v) {
		VertexListNode *node = NULL;
		if (getAdjacentListNode(u, v, &node))
			node->setMark();
	}

	void TriMeshDataStructure::unSetMark(Vertex *u, Vertex *v) {
		VertexListNode *node = NULL;
		if (getAdjacentListNode(u, v, &node))
			node->unSetMark();
	}

	bool TriMeshDataStructure::isMarked(Vertex *u, Vertex *v) const {
		VertexListNode *node = NULL;
		return getAdjacentListNode(u, v, &node) && node->isMarked();
	}

	bool TriMeshDataStructure::Adjacent(const Segment &s, Vertex **w) const {
		VertexListNode *node = NULL;
		if (getAdjacentListNode(s.v[0], s.v[1], &node) && !node->isPreFaceDeleted()) {
			*w = node->getVertex();
			return true;
		}
		return false;
	}

	bool TriMeshDataStructure::adjacent2Vertex(Vertex *w, Face *f) const {
		auto entry = topology.find(w);
		if (entry == topology.end())
			return false;
		VertexListNode *node = entry->second;
		VertexListNode *nextNode = NULL;
		bool found = false;
		while (node != NULL && !found) {
			nextNode = node->getNextNode();
			if (nextNode && !nextNode->isPreFaceDeleted()) {
				*f = Face(w, node->getVertex(), nextNode->getVertex());
				found = true;
			}
			node = nextNode;
		}
		return found;
	}

	bool TriMeshDataStructure::findIntersectedFace(Vertex *a, const DelVector& bb, Face *f) const {
		auto entry = topology.find(a);
		Assert(entry != topology.end());

		DelVector aa = a->vert;
		const Predicator<REAL> predicator;

		VertexListNode *parent = entry->second;
		if (parent == NULL) return false;
		Vertex *c = parent->getVertex();

		VertexListNode *child = parent->getNextNode();
		if (child == NULL) return false;
		Vertex *d = child->getVertex();

		if(c->isGhost() || d->isGhost()) {
			parent = child->getNextNode();
			if (parent == NULL) return false;

			child = parent->getNextNode();
			if (child == NULL) return false;

			c = parent->getVertex();
			d = child->getVertex();
		}

		REAL ori0 = predicator.orientCoplane(aa, bb, c->vert);
		REAL ori1 = predicator.orientCoplane(aa, bb, d->vert);
		Assert(ori0 != 0);
		Assert(ori1 != 0);
		while (child != NULL && (ori0 > 0) == (ori1 > 0)) {
			parent = child;
			child = parent->getNextNode();
			c = d;
			d = child->getVertex();
			if (!d->isGhost()) {
				ori0 = ori1;
				ori1 = predicator.orientCoplane(aa, bb, d->vert);
			}
			else {
				parent = child->getNextNode();
				if (parent == NULL) return false;

				child = parent->getNextNode();
				if (child == NULL) return false;

				c = parent->getVertex();
				d = child->getVertex();
				ori0 = predicator.orientCoplane(aa, bb, c->vert);
				ori1 = predicator.orientCoplane(aa, bb, d->vert);
				Assert(ori0 != 0);
			}
			Assert(ori1 != 0);
		}

		if (child == NULL) return false;

		*f = Face(a, c, d);

		return true;
	}

	std::vector<Face> TriMeshDataStructure::getTriangles(bool ghost) const {
		std::vector<Face> output;
		getTriangles(ghost, output);
		return output;
	}

	void TriMeshDataStructure::getTriangles(bool ghost, std::vector<Face>& triangles) const {
		VertexListNode *parent = NULL;
		VertexListNode *child = NULL;

		for (auto entry : topology) {
			Vertex *center = entry.first;
			if ((ghost || !center->isGhost()) && entry.second) {
				parent = entry.second;
				child = parent->getNextNode();
				while (child != NULL) {
					Vertex *b = parent->getVertex();
					Vertex *c = child->getVertex();
					if (!child->isPreFaceDeleted() &&
						*center < *b && *center < *c)
						triangles.push_back(Face(center, b, c));

					parent = child;
					child = child->getNextNode();
				}
				if (!entry.second->isPreFaceDeleted()) {
					Vertex *b = parent->getVertex();
					Vertex *c = entry.second->getVertex();
					if (*center < *b && *center < *c)
						triangles.push_back(Face(center, b, c));
				}
			}
		}
	}

	bool TriMeshDataStructure::getAdjacentListNode(Vertex* u, Vertex* v, VertexListNode **w) const {
		auto entry = topology.find(u);
		if (entry == topology.end()) return false;

		VertexListNode *node = entry->second;
		bool found = false;
		while (node != NULL) {
			if (node->getVertex() == v) {
				VertexListNode *nextNode = node->getNextNode();
				if (nextNode == NULL)
					nextNode = entry->second;
				*w = nextNode;
				found = true;
				break;
			}
			node = node->getNextNode();
		}

		return found;
	}

	void TriMeshDataStructure::insertToTopology(Vertex *a, Vertex *b, Vertex *c) {
		auto entry = topology.find(a);
		if (entry == topology.end()) {
			VertexListNode *node = nodePool->Alloc();
			node->setVertex(b);
			node->setDeletedMark();
			VertexListNode *nextNode = nodePool->Alloc();
			nextNode->setVertex(c);
			node->setNextNode(nextNode);
			topology[a] = node;
		}
		else {
			VertexListNode *head = entry->second;
			VertexListNode *parent = head;
			VertexListNode *child = parent->getNextNode();
			VertexListNode *foundNode[2];
			VertexListNode *parentFoundNode;
			int found = 0;
			if (head->getVertex() == b || head->getVertex() == c) {
				foundNode[found++] = head;
			}
			while (child != NULL && found < 2) {
				Vertex *foundVert = child->getVertex();
				if (foundVert == b)
					foundNode[found++] = child;
				else if (foundVert == c) {
					parentFoundNode = parent;
					foundNode[found++] = child;
				}
				parent = child;
				child = child->getNextNode();
			}

			switch (found) {
			case 0:
			{
				VertexListNode *newNode = nodePool->Alloc();
				newNode->setVertex(b);
				newNode->setDeletedMark();
				VertexListNode *nextNode = nodePool->Alloc();
				nextNode->setVertex(c);
				newNode->setNextNode(nextNode);

				if (head->isPreFaceDeleted()) {
					nextNode->setNextNode(head);
					head = newNode;
					entry->second = head;
				}
				else {
					parent = head;
					child = parent->getNextNode();
					while (child && !child->isPreFaceDeleted()) {
						parent = child;
						child = child->getNextNode();
					}
					Assert(child != NULL);
					nextNode->setNextNode(child);
					parent->setNextNode(newNode);
				}
				break;
			}
			case 1:
			{
				VertexListNode *newNode = nodePool->Alloc();
				if (foundNode[0]->isPreFaceDeleted()) {
					newNode->setVertex(b);
					if (foundNode[0] == head) {
						head->unSetDeletedMark();
						newNode->setDeletedMark();
						newNode->setNextNode(head);
						head = newNode;
						entry->second = head;
					}
					else {
						foundNode[0]->unSetDeletedMark();
						newNode->setNextNode(foundNode[0]);
						newNode->setDeletedMark();
						parentFoundNode->setNextNode(newNode);
					}
				}
				else {
					newNode->setVertex(c);
					newNode->setNextNode(foundNode[0]->getNextNode());
					foundNode[0]->setNextNode(newNode);
				}
				break;
			}
			case 2:
			{
				VertexListNode *bNode = NULL, *cNode = NULL;
				bool unoredered = (foundNode[0]->getVertex() == c);
				if (unoredered) {
					bNode = foundNode[1];  cNode = foundNode[0];
				}
				else {
					bNode = foundNode[0];  cNode = foundNode[1];
				}
				if (bNode->getNextNode() == cNode || (cNode == head && bNode->getNextNode() == NULL))
					cNode->unSetDeletedMark();
				else {
					//find the lost fan
					VertexListNode *start = bNode->getNextNode();
					VertexListNode *end = parentFoundNode;

					if (unoredered) {
						if (head == cNode) {
							VertexListNode *loop = start;
							while (loop != NULL) {
								end = loop;
								loop = loop->getNextNode();
							}
							bNode->setNextNode(NULL);
						}
						else {
							//the lost fan incule tail->head
							VertexListNode *loop = start;
							if (loop == NULL) {
								start = head;
								loop = start;
							}
							else {
								while (loop != NULL) {
									end = loop;
									loop = loop->getNextNode();
								}
								end->setNextNode(head);
								loop = head;
							}

							while (loop != cNode) {
								end = loop;
								loop = loop->getNextNode();
							}
							bNode->setNextNode(NULL);
							end->setNextNode(NULL);
							entry->second = cNode;
							head = cNode;
						}
					}

					else
						bNode->setNextNode(cNode);

					cNode->unSetDeletedMark();

					//find the hole to insert the lost fan
					VertexListNode *gapEnd = head;
					VertexListNode *gapStart = NULL;
					while (gapEnd && !gapEnd->isPreFaceDeleted()) {
						gapStart = gapEnd;
						gapEnd = gapEnd->getNextNode();
					}
					//insert it
					if (gapEnd != NULL) {
						if (gapStart != NULL) {
							gapStart->setNextNode(start);
							end->setNextNode(gapEnd);
						}
						else {
							end->setNextNode(head);
							head = start;
							entry->second = head;
						}
					}
					//the lost fan lost
					else {
						VertexListNode *loop = start;
						while (loop && loop != end) {
							VertexListNode *deleted = loop;
							loop = loop->getNextNode();
							nodePool->Dealloc(deleted);
						}
						nodePool->Dealloc(end);
					}
				}
				break;
			}
			default:
				Severe("Unexpted case in TriMeshDataStructure::insertToTopology");
				break;
			}
		}
	}

	void TriMeshDataStructure::removeFromTopology(Vertex *a, Vertex *b, Vertex *c) {
		auto entry = topology.find(a);
		Assert(entry != topology.end());
		VertexListNode *head = entry->second;
		if (head) {
			VertexListNode *grandparent = head;
			VertexListNode *parent = grandparent->getNextNode();
			VertexListNode *child = parent->getNextNode();

			if (grandparent->getVertex() == b && parent->getVertex() == c) {
				//tail case
				if (grandparent->isPreFaceDeleted()) {
					nodePool->Dealloc(grandparent);
					//check the node behide parent
					if (child && child->isPreFaceDeleted()) {
						entry->second = child;
						nodePool->Dealloc(parent);
					}
					//normal case
					else {
						parent->setDeletedMark();
						entry->second = parent;
					}
				}
				//head case
				else {
					parent->setDeletedMark();
					//check the node behide parent
					if (child && child->isPreFaceDeleted()) {
						grandparent->setNextNode(child);
						nodePool->Dealloc(parent);
					}
				}
			}
			else {
				bool found = false;
				while (child != NULL) {
					if (parent->getVertex() == b && child->getVertex() == c) {
						found = true;
						child->setDeletedMark();
						VertexListNode *grandchild = child->getNextNode();
						//head and tail case
						if (grandchild == NULL && head->isPreFaceDeleted()) {
							if (parent->isPreFaceDeleted()) {
								grandparent->setNextNode(NULL);
								nodePool->Dealloc(parent);
								parent = NULL;
							}
							else
								parent->setNextNode(NULL);

							nodePool->Dealloc(child);
							child = NULL;
						}
						//parent and child are all deleted
						else if (parent->isPreFaceDeleted()) {
							if (grandchild && grandchild->isPreFaceDeleted()) {
								grandparent->setNextNode(grandchild);
								nodePool->Dealloc(child);
								nodePool->Dealloc(parent);
								child = NULL;
								parent = NULL;
							}
							else {
								grandparent->setNextNode(child);
								nodePool->Dealloc(parent);
								parent = NULL;
							}
						}
						//child and grandchild are all deletd
						else if (grandchild && grandchild->isPreFaceDeleted()) {
							parent->setNextNode(grandchild);
							nodePool->Dealloc(child);
							child = NULL;
						}
						break;
					}
					grandparent = parent;
					parent = child;
					child = child->getNextNode();
				}

				if (!found) {
					if (parent->getVertex() == b && head->getVertex() == c) {
						VertexListNode *next = head->getNextNode();
						//the node behide the head deleted
						if (next && next->isPreFaceDeleted()) {
							next->setDeletedMark();
							entry->second = next;
							nodePool->Dealloc(head);
							head = next;
						}
						//normal case
						else {
							head->setDeletedMark();
							//check the tail
							if (parent->isPreFaceDeleted()) {
								grandparent->setNextNode(NULL);
								nodePool->Dealloc(parent);
							}
						}
					}
				}
			}
			//clean the deleted node
			if (head->getNextNode() == NULL) {
				nodePool->Dealloc(head);
				entry->second = NULL;
			}
		}
	}

	void TriMeshDataStructure::Clear() {
		topology.clear();
		nodePool->freeAll();
	}

	TriMeshDataStructure::~TriMeshDataStructure() {
		delete nodePool;
	}

	TetMeshDataStructure::TetMeshDataStructure() {
		nodePool = new MemoryPool<VertexListNode>(256);
		edgeNodePool = new MemoryPool<EdgeListNode>();
	}

	void TetMeshDataStructure::addTetrahedron(Vertex *a, Vertex *b, Vertex *c, Vertex *d) {
		bool ab = parityCheck(a, b);
		bool ac = parityCheck(a, c);
		bool ad = parityCheck(a, d);

		int condition = ab + (ac << 1) + (ad << 2);

		switch (condition) {
		case 0:
			insertToTopology(Segment(b, c), a, d);
			insertToTopology(Segment(b, d), c, a);
			insertToTopology(Segment(c, d), a, b);

			addSupplyVerts(b, c, d, a, 1);
			break;
		case 1:
			insertToTopology(Segment(a, b), c, d);
			insertToTopology(Segment(c, d), a, b);

			addSupplyVerts(a, b, c, d, 0);
			break;
		case 2:
			insertToTopology(Segment(a, c), d, b);
			insertToTopology(Segment(b, d), c, a);

			addSupplyVerts(a, c, b, d, 0);
			break;
		case 3:
			insertToTopology(Segment(a, b), c, d);
			insertToTopology(Segment(a, c), d, b);
			insertToTopology(Segment(b, c), a, d);

			addSupplyVerts(a, b, c, d, 1);
			break;
		case 4:
			insertToTopology(Segment(a, d), b, c);
			insertToTopology(Segment(b, c), a, d);

			addSupplyVerts(a, d, b, c, 0);
			break;
		case 5:
			insertToTopology(Segment(a, b), c, d);
			insertToTopology(Segment(a, d), b, c);
			insertToTopology(Segment(b, d), c, a);

			addSupplyVerts(a, b, d, c, 1);
			break;
		case 6:
			insertToTopology(Segment(a, c), d, b);
			insertToTopology(Segment(a, d), b, c);
			insertToTopology(Segment(c, d), a, b);

			addSupplyVerts(a, c, d, b, 1);
			break;
		case 7:
			insertToTopology(Segment(a, b), c, d);
			insertToTopology(Segment(a, c), d, b);
			insertToTopology(Segment(a, d), b, c);
			insertToTopology(Segment(b, c), a, d);
			insertToTopology(Segment(b, d), c, a);
			insertToTopology(Segment(c, d), a, b);

			addSupplyVerts(a, b, c, d, 2);
			break;
		default:
			Severe("Unexpected Case in TetMeshDataStructure::addTetrahedron");
			break;
		}
	}

	void TetMeshDataStructure::deleteTetrahedron(Vertex *a, Vertex *b, Vertex *c, Vertex *d) {
		bool ab = parityCheck(a, b);
		bool ac = parityCheck(a, c);
		bool ad = parityCheck(a, d);

		int condition = ab + (ac << 1) + (ad << 2);

		switch (condition) {
		case 0:
			removeFromTopology(Segment(b, c), a, d);
			removeFromTopology(Segment(b, d), c, a);
			removeFromTopology(Segment(c, d), a, b);
			break;
		case 1:
			removeFromTopology(Segment(a, b), c, d);
			removeFromTopology(Segment(c, d), a, b);
			break;
		case 2:
			removeFromTopology(Segment(a, c), d, b);
			removeFromTopology(Segment(b, d), c, a);
			break;
		case 3:
			removeFromTopology(Segment(a, b), c, d);
			removeFromTopology(Segment(a, c), d, b);
			removeFromTopology(Segment(b, c), a, d);
			break;
		case 4:
			removeFromTopology(Segment(a, d), b, c);
			removeFromTopology(Segment(b, c), a, d);
			break;
		case 5:
			removeFromTopology(Segment(a, b), c, d);
			removeFromTopology(Segment(a, d), b, c);
			removeFromTopology(Segment(b, d), c, a);
			break;
		case 6:
			removeFromTopology(Segment(a, c), d, b);
			removeFromTopology(Segment(a, d), b, c);
			removeFromTopology(Segment(c, d), a, b);
			break;
		case 7:
			removeFromTopology(Segment(a, b), c, d);
			removeFromTopology(Segment(a, c), d, b);
			removeFromTopology(Segment(a, d), b, c);
			removeFromTopology(Segment(b, c), a, d);
			removeFromTopology(Segment(b, d), c, a);
			removeFromTopology(Segment(c, d), a, b);
			break;
		default:
			Severe("Unexpected Case in TetMeshDataStructure::deleteTetrahedron");
			break;
		}
	}

	void TetMeshDataStructure::setMark(Vertex *u, Vertex *v, Vertex *w) {
		Face f(u, v, w, true);

		VertexListNode *node = NULL;
		if (getAdjacentListNode(f, &node))
			node->setMark();
	}

	void TetMeshDataStructure::unSetMark(Vertex *u, Vertex *v, Vertex *w) {
		Face f(u, v, w, true);

		VertexListNode *node = NULL;
		if (getAdjacentListNode(f, &node))
			node->unSetMark();
	}

	bool TetMeshDataStructure::testAndMark(Vertex *u, Vertex *v, Vertex *w) {
		Face f(u, v, w, true);

		VertexListNode *node = NULL;
		if (getAdjacentListNode(f, &node)) {
			if (node->isMarked()) return true;
			node->setMark();
		}
		return false;
	}

	bool TetMeshDataStructure::isMarked(Vertex *u, Vertex *v, Vertex *w) const {
		Face f(u, v, w, true);
		VertexListNode *node = NULL;

		return getAdjacentListNode(f, &node) && node->isMarked();
	}


	bool TetMeshDataStructure::Adjacent(const Face &f, Vertex **z) const {
		VertexListNode *node = NULL;
		if (getAdjacentListNode(f, &node)) {
			*z = node->getVertex();
			return true;
		}
		return false;
	}

	bool TetMeshDataStructure::getAdjacentListNode(const Face& f, VertexListNode **z) const {
		bool found = false;
		bool ab = parityCheck(f.v[0], f.v[1]);
		bool ac = parityCheck(f.v[0], f.v[2]);

		VertexListNode *node = NULL;
		if (ab)
			found = getAdjacentListNode(f.v[0], f.v[1], f.v[2], &node);
		else if (ac)
			found = getAdjacentListNode(f.v[2], f.v[0], f.v[1], &node);
		else
			found = getAdjacentListNode(f.v[1], f.v[2], f.v[0], &node);

		return found;
	}

	bool TetMeshDataStructure::getAdjacentListNode(Vertex *w, Vertex *x, Vertex *y, VertexListNode **z) const {
		bool found = false;
		if (edgeOrderCheck(x, w)) {
			if (x->hasList()) {
				EdgeListNode *linkHead = x->getListHead();
				while (linkHead != NULL && linkHead->getEndVertex() != w) {
					linkHead = linkHead->getNextNode();
				}
				if (linkHead != NULL) {
					VertexListNode *node = linkHead->getLink();
					while (node != NULL && !found) {
						if (node->getVertex() == y) {
							VertexListNode *foundNode = node->getNextNode();
							if (foundNode == NULL)
								foundNode = linkHead->getLink();
							if (!foundNode->isPreFaceDeleted()) {
								*z = foundNode;
								found = true;
							}
						}
						node = node->getNextNode();
					}
				}
			}
		}
		else {
			if (w->hasList()) {
				EdgeListNode *linkHead = w->getListHead();
				while (linkHead != NULL && linkHead->getEndVertex() != x) {
					linkHead = linkHead->getNextNode();
				}
				if (linkHead != NULL) {
					VertexListNode *parentNode = NULL;
					VertexListNode *node = linkHead->getLink();
					while (node != NULL && !found) {
						if (node->getVertex() == y) {
							if (!node->isPreFaceDeleted()) {
								VertexListNode *foundNode = parentNode;
								if (foundNode == NULL) {
									VertexListNode *loop = node;
									while (loop != NULL) {
										foundNode = loop;
										loop = loop->getNextNode();
									}
								}
								*z = foundNode;
								found = true;
							}
						}
						parentNode = node;
						node = node->getNextNode();
					}
				}
			}
		}
		return found;
	}

	bool TetMeshDataStructure::adjacent2Vertex(Vertex *w, Tetrahedron *t) const {
		bool found = false;
		EdgeListNode *linkHead = NULL;
		bool hasList = w->hasList();
		if (!hasList) {
			Vertex *end = w->getPointedVertex();
			linkHead = end->getListHead();
			bool foundHead = false;
			while (linkHead != NULL && !foundHead) {
				if (linkHead->getEndVertex() == w)
					foundHead = true;
				else
					linkHead = linkHead->getNextNode();
			}
		}
		else
			linkHead = w->getListHead();

		if (linkHead != NULL) {
			VertexListNode *head = linkHead->getLink();
			VertexListNode *parentLoop = linkHead->getLink();
			VertexListNode *loop = parentLoop->getNextNode();
			while (loop && loop->isPreFaceDeleted()) {
				parentLoop = loop;
				loop = loop->getNextNode();
			}
			if (loop == NULL && !head->isPreFaceDeleted()) {
				parentLoop = loop;
				loop = head;
			}
			if (loop) {
				if(hasList)
				    *t = Tetrahedron(w, linkHead->getEndVertex(), parentLoop->getVertex(), loop->getVertex());
				else
					*t = Tetrahedron(w, w->getPointedVertex(), loop->getVertex(), parentLoop->getVertex());
				found = true;
			}
		}
		else {
			if (!hasList) {
				Vertex *oppo = w->getPointedVertex();
				linkHead = oppo->getListHead();
				while (linkHead && !found) {
					VertexListNode *head = linkHead->getLink();
					VertexListNode *parentLoop = head;
					VertexListNode *loop = parentLoop->getNextNode();
					while (loop->getVertex() != w) {
						parentLoop = loop;
						loop = loop->getNextNode();
					}
					if (loop == NULL && head->getVertex() == w) {
						parentLoop = loop;
						loop = head;
					}
					if (loop) {
						VertexListNode *childLoop = loop->getNextNode();
						if (childLoop == NULL) childLoop = head;
						if (!loop->isPreFaceDeleted()) {
							*t = Tetrahedron(w, parentLoop->getVertex(), oppo, linkHead->getEndVertex());
							found = true;
						}
						else if (!childLoop->isPreFaceDeleted()) {
							*t = Tetrahedron(w, childLoop->getVertex(), linkHead->getEndVertex(), oppo);
							found = true;
						}
					}
					linkHead = linkHead->getNextNode();
				}
			}
		}

		return found;
	}

	void TetMeshDataStructure::getTetraherons(bool ghost, std::vector<Tetrahedron>& tets) const {
		for (auto vert : vertices) {
			EdgeListNode *linkHead = vert->getListHead();
			while (linkHead) {
				Vertex *endVert = linkHead->getEndVertex();
				VertexListNode *head = linkHead->getLink();
				VertexListNode *parentNode = head;
				VertexListNode *node = parentNode->getNextNode();
				while (node) {
					if (!node->isPreFaceDeleted()) {
						Vertex *c = parentNode->getVertex();
						Vertex *d = node->getVertex();
						if ((!parityCheck(c, vert) || (*vert < *c && *endVert < *c)) &&
							(!parityCheck(d, vert) || (*vert < *d && *endVert < *d)) &&
							(!parityCheck(c, d) || (*vert < *c && *vert < *d) || (*endVert < *c && *endVert < *d))) {
							if (ghost || (!endVert->isGhost() && !c->isGhost() && !d->isGhost()))
								tets.push_back(Tetrahedron(vert, endVert, c, d));
						}
					}
					parentNode = node;
					node = node->getNextNode();
				}
				if (!head->isPreFaceDeleted()) {
					Vertex *c = parentNode->getVertex();
					Vertex *d = node->getVertex();
					if ((!parityCheck(c, vert) || (*vert < *c && *endVert < *c)) &&
						(!parityCheck(d, vert) || (*vert < *d && *endVert < *d)) &&
						(!parityCheck(c, d) || (*vert < *c && *vert < *d) || (*endVert < *c && *endVert < *d))) {
						if (ghost || (!endVert->isGhost() && !c->isGhost() && !d->isGhost()))
							tets.push_back(Tetrahedron(vert, endVert, c, d));
					}
				}
				linkHead = linkHead->getNextNode();
			}
		}
	}

	std::vector<Tetrahedron> TetMeshDataStructure::getTetraherons(bool ghost) const {
		std::vector<Tetrahedron> tets;
		getTetraherons(ghost, tets);
		return tets;
	}

	void TetMeshDataStructure::Clear() {
		vertices.clear();
		nodePool->freeAll();
		edgeNodePool->freeAll();
	}

	void TetMeshDataStructure::insertToTopology(const Segment& s, Vertex *mayC, Vertex *mayD) {
		Vertex *a = s.v[0], *b = s.v[1], *c = mayC, *d = mayD;
		if (!edgeOrderCheck(a, b)) {
			std::swap(a, b);
			std::swap(c, d);
		}

		if (!a->hasList()) {
			vertices.push_back(a);
			//init link
			VertexListNode *node = nodePool->Alloc();
			node->setVertex(c);
			node->setDeletedMark();
			VertexListNode *nextNode = nodePool->Alloc();
			nextNode->setVertex(d);
			node->setNextNode(nextNode);
			//init edge list
			EdgeListNode *newNode = edgeNodePool->Alloc();
			newNode->setEndVertex(b);
			newNode->setLink(node);
			a->setListPointer(newNode);
		}

		else {
			EdgeListNode *linkHead = a->getListHead();
			bool foundHead = false;
			while (linkHead != NULL && !foundHead) {
				if (linkHead->getEndVertex() == b)
					foundHead = true;
				else
					linkHead = linkHead->getNextNode();
			}
			if (!foundHead) {
				//init link
				VertexListNode *node = nodePool->Alloc();
				node->setVertex(c);
				node->setDeletedMark();
				VertexListNode *nextNode = nodePool->Alloc();
				nextNode->setVertex(d);
				node->setNextNode(nextNode);
				//init edge list
				EdgeListNode *newNode = edgeNodePool->Alloc();
				newNode->setEndVertex(b);
				newNode->setLink(node);
				newNode->setNextNode(a->getListHead());
				a->setListPointer(newNode);
			}
			else {
				VertexListNode *head = linkHead->getLink();
				VertexListNode *parent = head;
				VertexListNode *child = parent->getNextNode();
				VertexListNode *foundNode[2];
				VertexListNode *parentFoundNode = NULL;
				int found = 0;
				if (head->getVertex() == c || head->getVertex() == d) {
					foundNode[found++] = head;
				}
				while (child != NULL && found < 2) {
					Vertex *foundVert = child->getVertex();
					if (foundVert == c)
						foundNode[found++] = child;
					else if (foundVert == d) {
						parentFoundNode = parent;
						foundNode[found++] = child;
					}
					parent = child;
					child = child->getNextNode();
				}

				switch (found) {
				case 0:
				{
					VertexListNode *newNode = nodePool->Alloc();
					newNode->setVertex(c);
					newNode->setDeletedMark();
					VertexListNode *nextNode = nodePool->Alloc();
					nextNode->setVertex(d);
					newNode->setNextNode(nextNode);

					if (head->isPreFaceDeleted()) {
						nextNode->setNextNode(head);
						head = newNode;
						linkHead->setLink(head);
					}
					else {
						parent = head;
						child = parent->getNextNode();
						while (child && !child->isPreFaceDeleted()) {
							parent = child;
							child = child->getNextNode();
						}
						Assert(child != NULL);
						nextNode->setNextNode(child);
						parent->setNextNode(newNode);
					}
					break;
				}
				case 1:
				{
					VertexListNode *newNode = nodePool->Alloc();
					if (foundNode[0]->isPreFaceDeleted()) {
						Assert(foundNode[0]->getVertex() == d);
						newNode->setVertex(c);
						if (foundNode[0] == head) {
							head->unSetDeletedMark();
							newNode->setDeletedMark();
							newNode->setNextNode(head);
							head = newNode;
							linkHead->setLink(head);
						}
						else {
							foundNode[0]->unSetDeletedMark();
							newNode->setNextNode(foundNode[0]);
							newNode->setDeletedMark();
							parentFoundNode->setNextNode(newNode);
						}
					}
					else {
						Assert(foundNode[0]->getVertex() == c);
						newNode->setVertex(d);
						newNode->setNextNode(foundNode[0]->getNextNode());
						foundNode[0]->setNextNode(newNode);
					}
					break;
				}
				case 2:
				{
					VertexListNode *cNode = NULL, *dNode = NULL;
					bool unoredered = (foundNode[0]->getVertex() == d);
					if (unoredered) {
						cNode = foundNode[1];  dNode = foundNode[0];
					}
					else {
						cNode = foundNode[0];  dNode = foundNode[1];
					}
					if (cNode->getNextNode() == dNode || (dNode == head && cNode->getNextNode() == NULL))
						dNode->unSetDeletedMark();
					else {
						//find the lost fan
						VertexListNode *start = cNode->getNextNode();
						VertexListNode *end = parentFoundNode;

						if (unoredered) {
							if (head == dNode) {
								VertexListNode *loop = start;
								while (loop != NULL) {
									end = loop;
									loop = loop->getNextNode();
								}
								cNode->setNextNode(NULL);
							}
							else {
								//the lost fan incule tail->head
								VertexListNode *loop = start;
								if (loop == NULL) {
									start = head;
									loop = start;
								}
								else {
									while (loop != NULL) {
										end = loop;
										loop = loop->getNextNode();
									}
									end->setNextNode(head);
									loop = head;
								}

								while (loop != dNode) {
									end = loop;
									loop = loop->getNextNode();
								}
								cNode->setNextNode(NULL);
								end->setNextNode(NULL);
								linkHead->setLink(dNode);
								head = dNode;
							}
						}

						else
							cNode->setNextNode(dNode);

						dNode->unSetDeletedMark();

						//find the hole to insert the lost fan
						VertexListNode *gapEnd = head;
						VertexListNode *gapStart = NULL;
						while (gapEnd && !gapEnd->isPreFaceDeleted()) {
							gapStart = gapEnd;
							gapEnd = gapEnd->getNextNode();
						}
						//insert it
						if (gapEnd != NULL) {
							if (gapStart != NULL) {
								gapStart->setNextNode(start);
								end->setNextNode(gapEnd);
							}
							else {
								end->setNextNode(head);
								head = start;
								linkHead->setLink(head);
							}
						}
						//the lost fan lost forever
						else {
							VertexListNode *loop = start;
							while (loop && loop != end) {
								VertexListNode *deleted = loop;
								loop = loop->getNextNode();
								nodePool->Dealloc(deleted);
							}
							nodePool->Dealloc(end);
						}
					}
					break;
				}
				default:
					Severe("Unexpted case in DelTriangulator::insertToTopology");
					break;
				}
			}
		}
	}

	void TetMeshDataStructure::removeFromTopology(const Segment &s, Vertex *mayC, Vertex *mayD) {
		Vertex *a = s.v[0], *b = s.v[1], *c = mayC, *d = mayD;
		if (!edgeOrderCheck(a, b)) {
			std::swap(a, b);
			std::swap(c, d);
		}

		Assert(a->hasList());
		EdgeListNode *parentLinkHead = NULL;
		EdgeListNode *linkHead = a->getListHead();
		bool foundHead = false;
		while (linkHead != NULL && !foundHead) {
			if (linkHead->getEndVertex() == b)
				foundHead = true;
			else {
				parentLinkHead = linkHead;
				linkHead = linkHead->getNextNode();
			}
		}
		Assert(foundHead);
		if (foundHead) {
			VertexListNode *head = linkHead->getLink();
			VertexListNode *grandparent = head;
			VertexListNode *parent = grandparent->getNextNode();
			VertexListNode *child = parent->getNextNode();

			if (grandparent->getVertex() == c && parent->getVertex() == d) {
				//tail case
				if (grandparent->isPreFaceDeleted()) {
					nodePool->Dealloc(grandparent);
					//check the node behide parent
					if (child && child->isPreFaceDeleted()) {
						head = child;
						linkHead->setLink(head);
						nodePool->Dealloc(parent);
					}
					else {
						parent->setDeletedMark();
						head = parent;
						linkHead->setLink(head);
					}
				}
				//head case
				else {
					parent->setDeletedMark();
					//check the node behide parent
					if (child && child->isPreFaceDeleted()) {
						grandparent->setNextNode(child);
						nodePool->Dealloc(parent);
					}
				}
			}
			else {
				bool found = false;
				while (child != NULL) {
					if (parent->getVertex() == c && child->getVertex() == d) {
						found = true;
						child->setDeletedMark();
						VertexListNode *grandchild = child->getNextNode();
						//head and tail case
						if (grandchild == NULL && head->isPreFaceDeleted()) {
							if (parent->isPreFaceDeleted()) {
								grandparent->setNextNode(NULL);
								nodePool->Dealloc(parent);
								parent = NULL;
							}
							else
								parent->setNextNode(NULL);

							nodePool->Dealloc(child);
							child = NULL;
						}
						//parent and child are all deleted
						else if (parent->isPreFaceDeleted()) {
							if (grandchild && grandchild->isPreFaceDeleted()) {
								grandparent->setNextNode(grandchild);
								nodePool->Dealloc(child);
								nodePool->Dealloc(parent);
								child = NULL;
								parent = NULL;
							}
							else {
								grandparent->setNextNode(child);
								nodePool->Dealloc(parent);
								parent = NULL;
							}
						}
						//child and grandchild are all deletd
						else if (grandchild && grandchild->isPreFaceDeleted()) {
							parent->setNextNode(grandchild);
							nodePool->Dealloc(child);
							child = NULL;
						}
						break;
					}
					grandparent = parent;
					parent = child;
					child = child->getNextNode();
				}

				if (!found) {
					if (parent->getVertex() == c && head->getVertex() == d) {
						//check the tail
						if (parent->isPreFaceDeleted()) {
							grandparent->setNextNode(NULL);
							nodePool->Dealloc(parent);
						}
						//the node behide the head deleted
						VertexListNode *next = head->getNextNode();
						if (next && next->isPreFaceDeleted()) {
							linkHead->setLink(next);
							nodePool->Dealloc(head);
							head = next;
						}
						//normal case
						else
							head->setDeletedMark();
					}
				}
			}
			//clean the deleted node
			if (head->getNextNode() == NULL) {
				nodePool->Dealloc(head);
				linkHead->setLink(NULL);
				if (parentLinkHead != NULL)
					parentLinkHead->setNextNode(linkHead->getNextNode());
				else
					a->setListPointer(linkHead->getNextNode());
				edgeNodePool->Dealloc(linkHead);
			}
		}
	}

	void TetMeshDataStructure::addSupplyVerts(Vertex *a, Vertex *b, Vertex *c, Vertex *d, int mode) {
		switch (mode) {
		case 0:
			if (!a->hasList() && !a->isGhost()) a->setVertexPointer(b);
			if (!b->hasList() && !b->isGhost()) b->setVertexPointer(a);
			if (!c->hasList() && !c->isGhost()) c->setVertexPointer(d);
			if (!d->hasList() && !d->isGhost()) d->setVertexPointer(c);
			break;
		case 1:
			if (!a->hasList() && !a->isGhost()) a->setVertexPointer(b);
			if (!b->hasList() && !b->isGhost()) b->setVertexPointer(c);
			if (!c->hasList() && !c->isGhost()) c->setVertexPointer(a);
			if (!d->hasList() && !d->isGhost()) d->setVertexPointer(a);
			break;
		case 2:
			if (!a->hasList() && !a->isGhost()) a->setVertexPointer(b);
			if (!b->hasList() && !b->isGhost()) b->setVertexPointer(c);
			if (!c->hasList() && !c->isGhost()) c->setVertexPointer(d);
			if (!d->hasList() && !d->isGhost()) d->setVertexPointer(a);
			break;
		default:
			Severe("Unexpected mode in DelMesher::addSupplyVerts");
			break;
		}
	}

	TetMeshDataStructure::~TetMeshDataStructure() {
		delete nodePool;
		delete edgeNodePool;
	}
}
