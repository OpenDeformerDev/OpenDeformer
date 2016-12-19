#include "stdafx.h"
#include "predicate.h"
#include "delPrimitive.h"
#include "geometer.h"
#include "memory.h"
#include <numeric>

namespace ODER {
	void Face::sortVertices() {
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

	void Face::initVertices(Vertex *v0, Vertex *v1, Vertex *v2, bool ordered) {
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

	Tetrahedron::Tetrahedron(Vertex *v0, Vertex *v1, Vertex *v2, Vertex *v3, bool ordered) {
		reRation = r = REAL(0);
		v[0] = v0; v[1] = v1; v[2] = v2; v[3] = v3;
		if (ordered)
			sortVertices();
	}

	void Tetrahedron::setRationAndRadius() {
		if (v[0]->isGhost() || v[1]->isGhost() || v[2]->isGhost() || v[3]->isGhost()) {
			reRation = REAL(0);
			r = REAL(0);
			return;
		}

		DelVector da = v[0]->vert - v[3]->vert;
		DelVector db = v[1]->vert - v[3]->vert;
		DelVector dc = v[2]->vert - v[3]->vert;
		DelVector ca = v[0]->vert - v[2]->vert;
		DelVector cb = v[1]->vert - v[2]->vert;
		DelVector ba = v[0]->vert - v[1]->vert;

		Geometer::Orthosphere(v[0]->vert, v[0]->weight, v[1]->vert, v[1]->weight, v[2]->vert, v[2]->weight, 
			v[3]->vert, v[3]->weight, (DelVector *)NULL, &r);

		REAL minEdgeLength = sqrt(std::min({da.length2(), db.length2(), dc.length2(), ca.length2(), cb.length2(), ba.length2()}));
		REAL relaxedLength = std::max(minEdgeLength, std::min({ v[0]->relaxedInsetionRadius, v[1]->relaxedInsetionRadius,
			v[2]->relaxedInsetionRadius, v[3]->relaxedInsetionRadius }));

		reRation = r / relaxedLength;
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
		nodePool = new MemoryPool<TriVertexListNode>(64);
	}

	void TriMeshDataStructure::addTriangle(Vertex *a, Vertex *b, Vertex *c, int index) {
		insertToTopology(a, b, c, index);
		insertToTopology(b, c, a, index);
		insertToTopology(c, a, b, index);
	}

	void TriMeshDataStructure::deleteTriangle(Vertex *a, Vertex *b, Vertex *c) {
		removeFromTopology(a, b, c);
		removeFromTopology(b, c, a);
		removeFromTopology(c, a, b);
	}

	void TriMeshDataStructure::setDeletedMark(Vertex *u, Vertex *v) {
		TriVertexListNode *node = NULL;
		if (getAdjacentListNode(u, v, &node))
			node->setDeletedMark();
	}

	void TriMeshDataStructure::unSetDeletedMark(Vertex *u, Vertex *v) {
		TriVertexListNode *node = NULL;
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
		TriVertexListNode *node = NULL;
		if (getAdjacentListNode(u, v, &node))
			node->setMark();
	}

	void TriMeshDataStructure::unSetMark(Vertex *u, Vertex *v) {
		TriVertexListNode *node = NULL;
		if (getAdjacentListNode(u, v, &node))
			node->unSetMark();
	}

	bool TriMeshDataStructure::isMarked(Vertex *u, Vertex *v) const {
		TriVertexListNode *node = NULL;
		return getAdjacentListNode(u, v, &node) && node->isMarked();
	}

	bool TriMeshDataStructure::Adjacent(const Segment &s, Vertex **w, int *index) const {
		TriVertexListNode *node = NULL;
		if (getAdjacentListNode(s.v[0], s.v[1], &node) && !node->isPreFaceDeleted()) {
			*w = node->getVertex();
			if (index) *index = node->getIndex();
			return true;
		}
		return false;
	}

	bool TriMeshDataStructure::adjacent2Vertex(Vertex *w, Face *f) const {
		auto entry = topology.find(w);
		if (entry == topology.end())
			return false;
		TriVertexListNode *node = entry->second;
		TriVertexListNode *nextNode = NULL;
		bool found = false;
		while (node != NULL && !found) {
			nextNode = node->getNextNode();
			if (nextNode && !nextNode->isPreFaceDeleted()) {
				*f = Face(w, node->getVertex(), nextNode->getVertex(), nextNode->getIndex());
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
		constexpr Predicator<REAL> predicator;

		TriVertexListNode *parent = entry->second;
		TriVertexListNode *child = parent->getNextNode();
		Vertex *c = NULL, *d = NULL;
		
		while (child != NULL) {
			if (child && !child->isPreFaceDeleted()) {
				c = parent->getVertex();
				d = child->getVertex();
				if(!c->isGhost() && !d->isGhost()) break;
			}
			parent = child;
			child = child->getNextNode();
		}
		if (child == NULL) {
			child = entry->second;
			if (!child->isPreFaceDeleted()) {
				c = parent->getVertex();
				d = child->getVertex();
				if (c->isGhost() || d->isGhost()) return false;
			}
			else 
				return false;
		}

		REAL ori0 = predicator.orientCoplane(aa, bb, c->vert);
		REAL ori1 = predicator.orientCoplane(aa, bb, d->vert);
		while (true) {
			if ((ori0 > 0) != (ori1 > 0) || ori0 == 0 || ori1 == 0){
				REAL ori2 = predicator.orientCoplane(aa, c->vert, d->vert);
				REAL ori3 = predicator.orientCoplane(bb, c->vert, d->vert);
				if ((ori2 > 0) != (ori3 > 0) || ori2 == 0 || ori3 == 0)
					break;
			}

			bool nextSuccess = false;
			parent = child;
			child = parent->getNextNode();
			c = d;

			if (child == NULL) break;
			if (!child->isPreFaceDeleted()) {
				d = child->getVertex();
				if (!d->isGhost()) {
					nextSuccess = true;
					ori0 = ori1;
					ori1 = predicator.orientCoplane(aa, bb, d->vert);
				}
			}

			if(!nextSuccess) {
				child = child->getNextNode();
				while (child != NULL) {
					parent = child;
					child = parent->getNextNode();
					if (child && !child->isPreFaceDeleted()) {
						c = parent->getVertex();
						d = child->getVertex();
						if (!c->isGhost() && !d->isGhost()) break;
					}
				}

				if (child != NULL) {
					ori0 = predicator.orientCoplane(aa, bb, c->vert);
					ori1 = predicator.orientCoplane(aa, bb, d->vert);
				}
				else
					break;
			}
		}

		if (child == NULL) {
			child = entry->second;
			if (!child->isPreFaceDeleted()) {
				c = parent->getVertex();
				d = child->getVertex();
				if (c->isGhost() || d->isGhost()) return false;

				ori0 = predicator.orientCoplane(aa, bb, c->vert);
				ori1 = predicator.orientCoplane(aa, bb, d->vert);

				if ((ori0 > 0) == (ori1 > 0) && ori0 != 0 && ori1 != 0)
					return false;
				else {
					REAL ori2 = predicator.orientCoplane(aa, c->vert, d->vert);
					REAL ori3 = predicator.orientCoplane(bb, c->vert, d->vert);
					if ((ori2 > 0) == (ori3 > 0) && ori2 != 0 && ori3 != 0)
						return false;
				}
			}
			else
				return false;
		}

		*f = Face(a, c, d, child->getIndex());

		return true;
	}

	TriMeshDataStructure::TriMeshConstIterator::TriMeshConstIterator(const std::unordered_map<Vertex *, TriVertexListNode *, vertex_hash>::const_iterator& iter,
		const std::unordered_map<Vertex *, TriVertexListNode *, vertex_hash>::const_iterator& end) : topologyIter(iter), topologyEnd(end){
		while (topologyIter != topologyEnd) {
			Vertex *center = topologyIter->first;
			parent = topologyIter->second;
			child = parent->getNextNode();
			while (child != NULL) {
				if (!child->isPreFaceDeleted()) {
					Vertex *b = parent->getVertex();
					Vertex *c = child->getVertex();
					if (verticesOrderCheck(center, b, c)) {
						current = Face(center, b, c, child->getIndex());
						return;
					}
				}

				parent = child;
				child = parent->getNextNode();
			}

			child = topologyIter->second;
			if (!child->isPreFaceDeleted()) {
				Vertex *b = parent->getVertex();
				Vertex *c = child->getVertex();
				if (verticesOrderCheck(center, b, c)) {
					current = Face(center, b, c, child->getIndex());
					return;
				}
			}

			++topologyIter;
		}
		parent = NULL; child = NULL;
	}

	void TriMeshDataStructure::TriMeshConstIterator::findNext() {
		do {
			if (child != topologyIter->second) {
				parent = child;
				child = parent->getNextNode();
				if (child == NULL) child = topologyIter->second;
			}
			else {
				++topologyIter;

				if (topologyIter != topologyEnd) {
					parent = topologyIter->second;
					child = parent->getNextNode();
				}
				else {
					current = Face(NULL, NULL, NULL);
					parent = NULL; child = NULL;
					return;
				}
			}

			if (!child->isPreFaceDeleted()) {
				Vertex *center = topologyIter->first;
				Vertex *b = parent->getVertex();
				Vertex *c = child->getVertex();
				if (verticesOrderCheck(center, b, c)) {
					current =  Face(center, b, c, child->getIndex());
					return;
				}
			}
		} while (true);
	}

	std::vector<Face> TriMeshDataStructure::getTriangles(bool ghost) const {
		std::vector<Face> output;
		getTriangles(ghost, output);
		return output;
	}

	void TriMeshDataStructure::getTriangles(bool ghost, std::vector<Face>& triangles) const {
		TriVertexListNode *parent = NULL;
		TriVertexListNode *child = NULL;
		triangles.reserve(topology.size() * 2);

		for (auto entry : topology) {
			Vertex *center = entry.first;
			if ((ghost || !center->isGhost())) {
				parent = entry.second;
				child = parent->getNextNode();
				while (child != NULL) {
					if (!child->isPreFaceDeleted()) {
						Vertex *b = parent->getVertex();
						Vertex *c = child->getVertex();
						if (verticesOrderCheck(center, b, c))
							triangles.push_back(Face(center, b, c, child->getIndex()));
					}

					parent = child;
					child = child->getNextNode();
				}
				if (!entry.second->isPreFaceDeleted()) {
					child = entry.second;
					Vertex *b = parent->getVertex();
					Vertex *c = child->getVertex();
					if (verticesOrderCheck(center, b, c))
						triangles.push_back(Face(center, b, c, child->getIndex()));
				}
			}
		}
	}

	bool TriMeshDataStructure::getAdjacentListNode(Vertex* u, Vertex* v, TriVertexListNode **w) const {
		auto entry = topology.find(u);
		if (entry == topology.end()) return false;

		TriVertexListNode *node = entry->second;
		bool found = false;
		while (node != NULL) {
			if (node->getVertex() == v) {
				TriVertexListNode *nextNode = node->getNextNode();
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

	void TriMeshDataStructure::insertToTopology(Vertex *a, Vertex *b, Vertex *c, int index) {
		auto entry = topology.find(a);
		if (entry == topology.end()) {
			TriVertexListNode *node = nodePool->Alloc();
			node->setVertex(b);
			node->setDeletedMark();
			TriVertexListNode *nextNode = nodePool->Alloc();
			nextNode->setVertex(c);
			nextNode->setIndex(index);
			node->setNextNode(nextNode);
			topology[a] = node;
		}
		else {
			TriVertexListNode *head = entry->second;
			TriVertexListNode *parent = head;
			TriVertexListNode *child = parent->getNextNode();
			TriVertexListNode *foundNode[2];
			TriVertexListNode *parentFoundNode;
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
				TriVertexListNode *newNode = nodePool->Alloc();
				newNode->setVertex(b);
				newNode->setDeletedMark();
				TriVertexListNode *nextNode = nodePool->Alloc();
				nextNode->setVertex(c);
				nextNode->setIndex(index);
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
				TriVertexListNode *newNode = nodePool->Alloc();
				if (foundNode[0]->isPreFaceDeleted()) {
					newNode->setVertex(b);
					foundNode[0]->setIndex(index);
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
					newNode->setIndex(index);
					newNode->setNextNode(foundNode[0]->getNextNode());
					foundNode[0]->setNextNode(newNode);
				}
				break;
			}
			case 2:
			{
				TriVertexListNode *bNode = NULL, *cNode = NULL;
				bool unoredered = (foundNode[0]->getVertex() == c);
				if (unoredered) {
					bNode = foundNode[1];  cNode = foundNode[0];
					Assert(cNode->getNextNode() != bNode);
				}
				else {
					bNode = foundNode[0];  cNode = foundNode[1];
				}
				cNode->setIndex(index);
				if (bNode->getNextNode() == cNode || (cNode == head && bNode->getNextNode() == NULL))
					cNode->unSetDeletedMark();
				else {
					//find the lost fan
					TriVertexListNode *start = bNode->getNextNode();
					TriVertexListNode *end = parentFoundNode;

					if (unoredered) {
						if (head == cNode) {
							TriVertexListNode *loop = start;
							while (loop != NULL) {
								end = loop;
								loop = loop->getNextNode();
							}
							bNode->setNextNode(NULL);
						}
						else {
							//the lost fan incule tail->head
							TriVertexListNode *loop = start;
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
					TriVertexListNode *gapEnd = head;
					TriVertexListNode *gapStart = NULL;
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
						TriVertexListNode *loop = start;
						while (loop && loop != end) {
							TriVertexListNode *deleted = loop;
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
		TriVertexListNode *head = entry->second;
		if (head) {
			TriVertexListNode *grandparent = head;
			TriVertexListNode *parent = grandparent->getNextNode();
			TriVertexListNode *child = parent->getNextNode();

			if (grandparent->getVertex() == b && parent->getVertex() == c) {
				//tail case
				if (grandparent->isPreFaceDeleted()) {
					nodePool->Dealloc(grandparent);
					//check the node behide parent
					if (child && child->isPreFaceDeleted()) {
						nodePool->Dealloc(parent);
						entry->second = child;
						head = child;
					}
					//normal case
					else {
						parent->setDeletedMark();
						entry->second = parent;
						head = parent;
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
						TriVertexListNode *grandchild = child->getNextNode();
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
						TriVertexListNode *next = head->getNextNode();
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
				topology.erase(entry);
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
		nodePool = new MemoryPool<TetVertexListNode>(256);
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
			addEnforcedSegmentsAdjacency(b, d, c, a, 1);
			break;
		case 1:
			insertToTopology(Segment(a, b), c, d);
			insertToTopology(Segment(c, d), a, b);

			addSupplyVerts(a, b, c, d, 0);
			addEnforcedSegmentsAdjacency(a, b, c, d, 0);
			break;
		case 2:
			insertToTopology(Segment(a, c), d, b);
			insertToTopology(Segment(b, d), c, a);

			addSupplyVerts(a, c, b, d, 0);
			addEnforcedSegmentsAdjacency(a, c, d, b, 0);
			break;
		case 3:
			insertToTopology(Segment(a, b), c, d);
			insertToTopology(Segment(a, c), d, b);
			insertToTopology(Segment(b, c), a, d);

			addSupplyVerts(a, b, c, d, 1);
			addEnforcedSegmentsAdjacency(a, b, c, d, 1);
			break;
		case 4:
			insertToTopology(Segment(a, d), b, c);
			insertToTopology(Segment(b, c), a, d);

			addSupplyVerts(a, d, b, c, 0);
			addEnforcedSegmentsAdjacency(a, d, b, c, 0);
			break;
		case 5:
			insertToTopology(Segment(a, b), c, d);
			insertToTopology(Segment(a, d), b, c);
			insertToTopology(Segment(b, d), c, a);

			addSupplyVerts(a, b, d, c, 1);
			addEnforcedSegmentsAdjacency(a, d, b, c, 1);
			break;
		case 6:
			insertToTopology(Segment(a, c), d, b);
			insertToTopology(Segment(a, d), b, c);
			insertToTopology(Segment(c, d), a, b);

			addSupplyVerts(a, c, d, b, 1);
			addEnforcedSegmentsAdjacency(a, c, d, b, 1);
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

		TetVertexListNode *node = NULL;
		if (getAdjacentListNode(f, &node))
			node->setMark();
	}

	void TetMeshDataStructure::unSetMark(Vertex *u, Vertex *v, Vertex *w) {
		Face f(u, v, w, true);

		TetVertexListNode *node = NULL;
		if (getAdjacentListNode(f, &node))
			node->unSetMark();
	}

	bool TetMeshDataStructure::isMarked(Vertex *u, Vertex *v, Vertex *w) const {
		Face f(u, v, w, true);
		TetVertexListNode *node = NULL;

		return getAdjacentListNode(f, &node) && node->isMarked();
	}


	bool TetMeshDataStructure::Adjacent(const Face &f, Vertex **z) const {
		TetVertexListNode *node = NULL;
		if (getAdjacentListNode(f, &node)) {
			*z = node->getVertex();
			return true;
		}
		return false;
	}

	bool TetMeshDataStructure::getAdjacentListNode(const Face& f, TetVertexListNode **z) const {
		bool found = false;
		bool ab = parityCheck(f.v[0], f.v[1]);
		bool ac = parityCheck(f.v[0], f.v[2]);

		if (ab)
			found = getAdjacentListNode(f.v[0], f.v[1], f.v[2], z);
		else if (ac)
			found = getAdjacentListNode(f.v[2], f.v[0], f.v[1], z);
		else
			found = getAdjacentListNode(f.v[1], f.v[2], f.v[0], z);

		return found;
	}

	bool TetMeshDataStructure::getAdjacentListNode(Vertex *w, Vertex *x, Vertex *y, TetVertexListNode **z) const {
		bool found = false;
		if (edgeOrderCheck(x, w)) {
			auto pair = topology.find(x);
			if (pair != topology.end()) {
				EdgeListNode *linkHead = pair->second;
				while (linkHead != NULL && linkHead->getEndVertex() != w) {
					linkHead = linkHead->getNextNode();
				}
				if (linkHead != NULL) {
					TetVertexListNode *node = linkHead->getLink();
					while (node != NULL && !found) {
						if (node->getVertex() == y) {
							TetVertexListNode *foundNode = node->getNextNode();
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
			auto pair = topology.find(w);
			if (pair != topology.end()) {
				EdgeListNode *linkHead = pair->second;
				while (linkHead != NULL && linkHead->getEndVertex() != x) {
					linkHead = linkHead->getNextNode();
				}
				if (linkHead != NULL) {
					TetVertexListNode *parentNode = NULL;
					TetVertexListNode *node = linkHead->getLink();
					while (node != NULL && !found) {
						if (node->getVertex() == y) {
							if (!node->isPreFaceDeleted()) {
								TetVertexListNode *foundNode = parentNode;
								if (foundNode == NULL) {
									TetVertexListNode *loop = node;
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
		auto pair = topology.find(w);
		bool hasList = (pair != topology.end() && pair->second != NULL);
		if (!hasList) {
			Vertex *end = w->getPointedVertex();
			linkHead = topology.find(end)->second;
			bool foundHead = false;
			while (linkHead != NULL && !foundHead) {
				if (linkHead->getEndVertex() == w)
					foundHead = true;
				else
					linkHead = linkHead->getNextNode();
			}
		}
		else
			linkHead = pair->second;

		if (linkHead != NULL) {
			TetVertexListNode *head = linkHead->getLink();
			TetVertexListNode *parentLoop = linkHead->getLink();
			TetVertexListNode *loop = parentLoop->getNextNode();
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
				linkHead = topology.find(oppo)->second;
				while (linkHead && !found) {
					TetVertexListNode *head = linkHead->getLink();
					TetVertexListNode *parentLoop = head;
					TetVertexListNode *loop = parentLoop->getNextNode();
					while (loop && loop->getVertex() != w) {
						parentLoop = loop;
						loop = loop->getNextNode();
					}
					if (loop == NULL && head->getVertex() == w)
						loop = head;
					if (loop) {
						TetVertexListNode *childLoop = loop->getNextNode();
						if (childLoop == NULL) childLoop = head;
						if (!loop->isPreFaceDeleted()) {
							*t = Tetrahedron(w, parentLoop->getVertex(), linkHead->getEndVertex(), oppo);
							found = true;
						}
						else if (!childLoop->isPreFaceDeleted()) {
							*t = Tetrahedron(w, childLoop->getVertex(), oppo, linkHead->getEndVertex());
							found = true;
						}
					}
					linkHead = linkHead->getNextNode();
				}
			}
		}

		return found;
	}

	bool TetMeshDataStructure::adjacent2Segment(const Segment &s, Tetrahedron *t) const {
		Vertex *a = s.v[0], *b = s.v[1];
		if (parityCheck(a, b)) {
			bool ordered = edgeOrderCheck(a, b);
			if (!ordered) std::swap(a, b);
			auto pair = topology.find(a);
			if (pair != topology.end()) {
				EdgeListNode *linkHead = pair->second;
				while (linkHead != NULL && linkHead->getEndVertex() != b)
					linkHead = linkHead->getNextNode();

				if (linkHead != NULL) {
					TetVertexListNode *parentNode = linkHead->getLink();
					TetVertexListNode *node = parentNode != NULL ? parentNode->getNextNode() : NULL;
					while (node != NULL) {
						if (!node->isPreFaceDeleted()) {
							if (ordered) *t = Tetrahedron(a, b, parentNode->getVertex(), node->getVertex());
							else *t = Tetrahedron(b, a, node->getVertex(), parentNode->getVertex());
							return true;
						}
						parentNode = node;
						node = node->getNextNode();
					}
					node = linkHead->getLink();
					if (parentNode != node && !node->isPreFaceDeleted()) {
						if (ordered) *t = Tetrahedron(a, b, parentNode->getVertex(), node->getVertex());
						else *t = Tetrahedron(b, a, node->getVertex(), parentNode->getVertex());
						return true;
					}
				}
			}
		}
		else {
			auto found = segmentAdjacency.find(Segment(a, b, true));
			Vertex *c = found->second;
			if (found != segmentAdjacency.end() && c != NULL) {
				Vertex *d = NULL;
				if (Adjacent(Face(b, a, c), &d)) {
					*t = Tetrahedron(a, b, c, d);
					return true;
				}
			}
		}

		return false;
	}


	TetMeshDataStructure::TetMeshConstIterator::TetMeshConstIterator(const std::unordered_map<Vertex *, EdgeListNode *, vertex_hash>::const_iterator& iter,
		const std::unordered_map<Vertex *, EdgeListNode *, vertex_hash>::const_iterator& end) : topologyIter(iter), topologyEnd(end){
		while (topologyIter != topologyEnd) {
			Vertex *ori = topologyIter->first;
			linkHead = topologyIter->second;
			while (linkHead) {
				Vertex *end = linkHead->getEndVertex();
				parent = linkHead->getLink();
				child = parent->getNextNode();
				while (child) {
					if (!child->isPreFaceDeleted()) {
						Vertex *c = parent->getVertex();
						Vertex *d = child->getVertex();
						if (verticesOrderCheck(ori, end, c, d)) {
							current = Tetrahedron(ori, end, c, d);
							return;
						}
					}

					parent = child;
					child = parent->getNextNode();
				}
				child = linkHead->getLink();
				if (!child->isPreFaceDeleted()) {
					Vertex *c = parent->getVertex();
					Vertex *d = child->getVertex();
					if (verticesOrderCheck(ori, end, c, d)) {
						current = Tetrahedron(ori, end, c, d);
						return;
					}
				}
				linkHead = linkHead->getNextNode();
			}
			++topologyIter;
		}
		linkHead = NULL;
		parent = NULL; child = NULL;
	}

	void TetMeshDataStructure::TetMeshConstIterator::findNext() {
		do {
			if (child != linkHead->getLink()) {
				parent = child;
				child = parent->getNextNode();
				if (!child) child = linkHead->getLink();
			}
			else {
				linkHead = linkHead->getNextNode();
				if (!linkHead) {
					++topologyIter;
					if (topologyIter != topologyEnd)
						linkHead = topologyIter->second;
					else {
						current = Tetrahedron(NULL, NULL, NULL, NULL);
						linkHead = NULL;
						parent = NULL; child = NULL;
						return;
					}
				}
				parent = linkHead->getLink();
				child = parent->getNextNode();
			}

			if (!child->isPreFaceDeleted()) {
				Vertex *ori = topologyIter->first;
				Vertex *end = linkHead->getEndVertex();
				Vertex *c = parent->getVertex();
				Vertex *d = child->getVertex();
				if (verticesOrderCheck(ori, end, c, d)) {
					current = Tetrahedron(ori, end, c, d);
					return;
				}
			}
		} while (true);
	}

	void TetMeshDataStructure::getTetrahedrons(bool ghost, std::vector<Tetrahedron>& tets) const {
		for (auto pair : topology) {
			Vertex *vert = pair.first;
			EdgeListNode *linkHead = pair.second;
			while (linkHead) {
				Vertex *endVert = linkHead->getEndVertex();
				if (ghost || !endVert->isGhost()) {
					TetVertexListNode *head = linkHead->getLink();
					TetVertexListNode *parentNode = head;
					TetVertexListNode *node = parentNode->getNextNode();
					while (node) {
						if (!node->isPreFaceDeleted()) {
							Vertex *c = parentNode->getVertex();
							Vertex *d = node->getVertex();
							if (verticesOrderCheck(vert, endVert, c, d)) {
								if (ghost || (!c->isGhost() && !d->isGhost()))
									tets.push_back(Tetrahedron(vert, endVert, c, d));
							}
						}
						parentNode = node;
						node = node->getNextNode();
					}
					if (!head->isPreFaceDeleted()) {
						Vertex *c = parentNode->getVertex();
						Vertex *d = head->getVertex();
						if (verticesOrderCheck(vert, endVert, c, d)) {
							if (ghost || (!c->isGhost() && !d->isGhost()))
								tets.push_back(Tetrahedron(vert, endVert, c, d));
						}
					}
				}
				linkHead = linkHead->getNextNode();
			}
		}
	}

	std::vector<Tetrahedron> TetMeshDataStructure::getTetrahedrons(bool ghost) const {
		std::vector<Tetrahedron> tets;
		getTetrahedrons(ghost, tets);
		return tets;
	}

	void TetMeshDataStructure::Clear() {
		topology.clear();
		nodePool->freeAll();
		edgeNodePool->freeAll();
	}

	void TetMeshDataStructure::insertToTopology(const Segment& s, Vertex *mayC, Vertex *mayD) {
		Vertex *a = s.v[0], *b = s.v[1], *c = mayC, *d = mayD;
		if (!edgeOrderCheck(a, b)) {
			std::swap(a, b);
			std::swap(c, d);
		}

		auto pair = topology.find(a);
		if (pair == topology.end()) {
			//init link
			TetVertexListNode *node = nodePool->Alloc();
			node->setVertex(c);
			node->setDeletedMark();
			TetVertexListNode *nextNode = nodePool->Alloc();
			nextNode->setVertex(d);
			node->setNextNode(nextNode);
			//init edge list
			EdgeListNode *newNode = edgeNodePool->Alloc();
			newNode->setEndVertex(b);
			newNode->setLink(node);
			topology.insert(std::make_pair(a, newNode));
		}

		else {
			EdgeListNode *linkHead = pair->second;
			bool foundHead = false;
			while (linkHead != NULL && !foundHead) {
				if (linkHead->getEndVertex() == b)
					foundHead = true;
				else
					linkHead = linkHead->getNextNode();
			}
			if (!foundHead) {
				//init link
				TetVertexListNode *node = nodePool->Alloc();
				node->setVertex(c);
				node->setDeletedMark();
				TetVertexListNode *nextNode = nodePool->Alloc();
				nextNode->setVertex(d);
				node->setNextNode(nextNode);
				//init edge list
				EdgeListNode *newNode = edgeNodePool->Alloc();
				newNode->setEndVertex(b);
				newNode->setLink(node);
				newNode->setNextNode(pair->second);
				pair->second = newNode;
			}
			else {
				TetVertexListNode *head = linkHead->getLink();
				TetVertexListNode *parent = head;
				TetVertexListNode *child = parent->getNextNode();
				TetVertexListNode *foundNode[2];
				TetVertexListNode *parentFoundNode = NULL;
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
					TetVertexListNode *newNode = nodePool->Alloc();
					newNode->setVertex(c);
					newNode->setDeletedMark();
					TetVertexListNode *nextNode = nodePool->Alloc();
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
					TetVertexListNode *newNode = nodePool->Alloc();
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
					TetVertexListNode *cNode = NULL, *dNode = NULL;
					bool unoredered = (foundNode[0]->getVertex() == d);
					if (unoredered) {
						cNode = foundNode[1];  dNode = foundNode[0];
						Assert(dNode->getNextNode() != cNode);
					}
					else {
						cNode = foundNode[0];  dNode = foundNode[1];
					}
					if (cNode->getNextNode() == dNode || (dNode == head && cNode->getNextNode() == NULL))
						dNode->unSetDeletedMark();
					else {
						//find the lost fan
						TetVertexListNode *start = cNode->getNextNode();
						TetVertexListNode *end = parentFoundNode;

						if (unoredered) {
							if (head == dNode) {
								TetVertexListNode *loop = start;
								while (loop != NULL) {
									end = loop;
									loop = loop->getNextNode();
								}
								cNode->setNextNode(NULL);
							}
							else {
								//the lost fan incule tail->head
								TetVertexListNode *loop = start;
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
						TetVertexListNode *gapEnd = head;
						TetVertexListNode *gapStart = NULL;
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
							TetVertexListNode *loop = start;
							while (loop && loop != end) {
								TetVertexListNode *deleted = loop;
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

		auto pair = topology.find(a);
		Assert(pair != topology.end());
		EdgeListNode *parentLinkHead = NULL;
		EdgeListNode *linkHead = pair->second;
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
			TetVertexListNode *head = linkHead->getLink();
			TetVertexListNode *grandparent = head;
			TetVertexListNode *parent = grandparent->getNextNode();
			TetVertexListNode *child = parent->getNextNode();

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
						TetVertexListNode *grandchild = child->getNextNode();
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
						//child and grandchild are all deleted
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
						TetVertexListNode *next = head->getNextNode();
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
				else {
					EdgeListNode *newLinkHead = linkHead->getNextNode();
					if (newLinkHead) pair->second = newLinkHead;
					else topology.erase(pair);
				}
				edgeNodePool->Dealloc(linkHead);
			}
		}
	}

	void TetMeshDataStructure::addSupplyVerts(Vertex *a, Vertex *b, Vertex *c, Vertex *d, int mode) {
		switch (mode) {
		case 0:
			a->setVertexPointer(b);
			b->setVertexPointer(a);
			c->setVertexPointer(d);
			d->setVertexPointer(c);
			break;
		case 1:
			a->setVertexPointer(b);
			b->setVertexPointer(c);
			c->setVertexPointer(a);
			d->setVertexPointer(edgeOrderCheck(a, b) ? a : b);
			break;
		case 2:
			a->setVertexPointer(b);
			b->setVertexPointer(c);
			c->setVertexPointer(d);
			d->setVertexPointer(a);
			break;
		default:
			Severe("Unexpected mode in DelMesher::addSupplyVerts");
			break;
		}
	}

	void TetMeshDataStructure::addEnforcedSegmentsAdjacency(Vertex *a, Vertex *b, Vertex *c, Vertex *d, int mode) {
		switch (mode) {
		case 0:
			addEnforcedSegmentAdjacency(a, c, d, b);
			addEnforcedSegmentAdjacency(a, d, b, c);
			addEnforcedSegmentAdjacency(b, c, a, d);
			addEnforcedSegmentAdjacency(b, d, c, a);
			break;
		case 1:
			addEnforcedSegmentAdjacency(a, d, b, c);
			addEnforcedSegmentAdjacency(b, d, c, a);
			addEnforcedSegmentAdjacency(c, d, a, b);
			break;
		default:
			Severe("Unexpected mode in DelMesher::addEnforcedSegmentsAdjacency");
			break;
		}
	}

	void TetMeshDataStructure::addEnforcedSegmentAdjacency(Vertex *a, Vertex *b, Vertex *c, Vertex *d) {
		if (matchVertexFlag(a->type, VertexType::Vertex_Segment) && 
			matchVertexFlag(b->type, VertexType::Vertex_Segment)) {
			auto found = segmentAdjacency.find(Segment(a, b, true));
			if (found != segmentAdjacency.end()) {
				if (!edgeOrderCheck(a, b)) {
					std::swap(a, b);
					std::swap(c, d);
				}

				found->second = c;
			}
		}
	}

	TetMeshDataStructure::~TetMeshDataStructure() {
		delete nodePool;
		delete edgeNodePool;
	}
}
