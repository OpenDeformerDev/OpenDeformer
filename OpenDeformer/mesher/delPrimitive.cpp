#include "stdafx.h"
#include "delPrimitive.h"
#include "memory.h"
#include "predicate.h"

namespace ODER {
	Labeler Vertex::labeler;

	REAL Tetrahedron::maxREration = FLT_MAX;

	Tetrahedron::Tetrahedron(Vertex *v0, Vertex *v1, Vertex *v2, Vertex *v3, bool ordered) {
		reRation = r = 0.0;
		v[0] = v0; v[1] = v1; v[2] = v2; v[3] = v3;
		if (ordered)
			sortVertices();
	}

	void Tetrahedron::setRationAndRadius(const Vertex *ghost) {
		if (v[0] == ghost || v[1] == ghost || v[2] == ghost || v[3] == ghost) {
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

	bool TriMeshDataStructure::Adjacent(const Segment &s, Vertex **w) const {
		auto entry = topology.find(s.v[0]);
		if (entry == topology.end())
			return false;
		VertexListNode *node = entry->second;
		bool found = false;
		while (node != NULL) {
			if (node->getVertex() == s.v[1]) {
				VertexListNode *nextNode = node->getNextNode();
				if (nextNode == NULL)
					nextNode = entry->second;
				if (!nextNode->isPreFaceDeleted()) {
					*w = nextNode->getVertex();
					found = true;
				}
				break;
			}
			node = node->getNextNode();
		}
		return found;
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

	std::set<Face, face_compare> TriMeshDataStructure::getTriangles(const Vertex *ghost) const {
		std::set<Face, face_compare> output;
		VertexListNode *parent = NULL;
		VertexListNode *child = NULL;

		for (auto entry : topology) {
			Vertex *center = entry.first;
			if ((ghost == NULL || center != ghost) && entry.second) {
				parent = entry.second;
				child = parent->getNextNode();
				while (child != NULL) {
					Vertex *b = parent->getVertex();
					Vertex *c = child->getVertex();
					if (!child->isPreFaceDeleted()) {
						Face f(center, b, c, true);
						if (ghost == NULL || f.v[0] != ghost)
							output.insert(f);
					}
					parent = child;
					child = child->getNextNode();
				}
				if (!entry.second->isPreFaceDeleted()) {
					Vertex *b = parent->getVertex();
					Vertex *c = entry.second->getVertex();
					Face f(center, b, c, true);
					if (ghost == NULL || f.v[0] != ghost)
						output.insert(f);
				}
			}
		}
		return output;
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

			addSupplyVerts(b, c, d, NULL, 1);
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

			addSupplyVerts(a, b, c, NULL, 1);
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

			addSupplyVerts(a, b, d, NULL, 1);
			break;
		case 6:
			insertToTopology(Segment(a, c), d, b);
			insertToTopology(Segment(a, d), b, c);
			insertToTopology(Segment(c, d), a, b);

			addSupplyVerts(a, c, d, NULL, 1);
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
			Severe("Unexpected Case in DelMesher::addTet");
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
			Severe("Unexpected Case in DelMesher::deleteTet");
			break;
		}
	}

	bool TetMeshDataStructure::Adjacent(const Face &f, Vertex **z) const {
		bool found = false;
		bool ab = parityCheck(f.v[0], f.v[1]);
		bool ac = parityCheck(f.v[0], f.v[2]);

		if (ab)
			found = Adjacent(f.v[0], f.v[1], f.v[2], z);
		else if (ac)
			found = Adjacent(f.v[2], f.v[0], f.v[1], z);
		else
			found = Adjacent(f.v[1], f.v[2], f.v[0], z);

		//if (oldFound != found)
		return found;
	}

	bool TetMeshDataStructure::Adjacent(Vertex *w, Vertex *x, Vertex *y, Vertex **z) const {
		bool found = false;
		if (w->getLabel() > x->getLabel()) {
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
								*z = foundNode->getVertex();
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
								*z = foundNode->getVertex();
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
		if (!w->hasList()) {
			Vertex *end = w->getEndVertex();
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
			VertexListNode *parentLoop = linkHead->getLink();
			VertexListNode *loop = parentLoop->getNextNode();
			while (loop && loop->isPreFaceDeleted()) {
				parentLoop = loop;
				loop = loop->getNextNode();
			}
			if (loop) {
				*t = Tetrahedron(w, linkHead->getEndVertex(), parentLoop->getVertex(), loop->getVertex(), true);
				found = true;
			}
		}
		return found;
	}

	std::set<Tetrahedron, tet_compare> TetMeshDataStructure::getTetraherons(const Vertex *ghost) const {
		std::set<Tetrahedron, tet_compare> tets;
		for (auto vert : vertices) {
			if (vert->hasList()) {
				EdgeListNode *linkHead = vert->getListHead();
				while (linkHead) {
					Vertex *endVert = linkHead->getEndVertex();
					VertexListNode *head = linkHead->getLink();
					VertexListNode *parentNode = head;
					VertexListNode *node = parentNode->getNextNode();
					while (node) {
						if (!node->isPreFaceDeleted()) {
							Tetrahedron t(vert, endVert, parentNode->getVertex(), node->getVertex(), true);
							if (ghost == NULL || t.v[0] != ghost)
								tets.insert(t);
						}
						parentNode = node;
						node = node->getNextNode();
					}
					if (!head->isPreFaceDeleted()) {
						Tetrahedron t(vert, endVert, parentNode->getVertex(), head->getVertex(), true);
						if (ghost == NULL || t.v[0] != ghost)
							tets.insert(t);
					}
					linkHead = linkHead->getNextNode();
				}
			}
		}

		return tets;
	}

	void TetMeshDataStructure::insertToTopology(const Segment& s, Vertex *mayC, Vertex *mayD) {
		Vertex *a = s.v[0], *b = s.v[1], *c = mayC, *d = mayD;
		if (*a > *b) {
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
		if (*a > *b) {
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
			VertexListNode *loop = linkHead->getLink();
			int startnum = 0;
			while (loop) {
				loop = loop->getNextNode();
				startnum++;
			}

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
		int aIndex = a->getLabel();
		int bIndex = b->getLabel();
		int cIndex = c->getLabel();

		int dIndex = -1;
		if (d)
			dIndex = d->getLabel();

		switch (mode) {
		case 0:
			if (aIndex > bIndex) {
				if (!a->hasList())
					a->setEndVertexPoint(b);
			}
			else {
				if (!b->hasList())
					b->setEndVertexPoint(a);
			}
			if (cIndex > dIndex) {
				if (!c->hasList())
					c->setEndVertexPoint(d);
			}
			else {
				if (!d->hasList())
					d->setEndVertexPoint(c);
			}
			break;
		case 1:
			if (aIndex > bIndex) {
				if (aIndex > cIndex) {
					if (!a->hasList())
						a->setEndVertexPoint(c);
				}
				else {
					if (!b->hasList())
						b->setEndVertexPoint(c);
				}
			}
			else if (bIndex > cIndex) {
				if (!b->hasList())
					b->setEndVertexPoint(a);
			}
			else {
				if (!c->hasList())
					c->setEndVertexPoint(b);
			}
			break;
		case 2:
		{
			int v[4] = { aIndex, bIndex, cIndex, dIndex };
			int max = 0;
			for (int i = 1; i < 4; i++) {
				if (v[i] > v[max])
					max = i;
			}
			if (max == 0) {
				if (!a->hasList())
					a->setEndVertexPoint(b);
			}
			else if (max == 1) {
				if (!b->hasList())
					b->setEndVertexPoint(c);
			}
			else if (max == 2) {
				if (!c->hasList())
					c->setEndVertexPoint(d);
			}
			else {
				if (!d->hasList())
					d->setEndVertexPoint(a);
			}
			break;
		}
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
