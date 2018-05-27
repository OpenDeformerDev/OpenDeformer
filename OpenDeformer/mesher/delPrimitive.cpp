#include "stdafx.h"
#include "predicate.h"
#include "delPrimitive.h"
#include "geometer.h"

namespace ODER {
	void Triangle::initVertices(Vertex *v0, Vertex *v1, Vertex *v2, bool ordered) {
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

	void Triangle::sortVertices() {
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

	void TriangleWithGeometry::setGeometricProperties() {
		if (radius >= 0) return;
		if (v[0]->isGhost() || v[1]->isGhost() || v[2]->isGhost()) {
			radius = DelReal(0);
			return;
		}
		Geometer::Circumcircle(v[0]->point, v[1]->point, v[2]->point, &circumcenter, &radius);
	}

	void Tetrahedron::setGeometricProperties() {
		if (radius >= 0) return;

		if (v[0]->isGhost() || v[1]->isGhost() || v[2]->isGhost() || v[3]->isGhost()) {
			reRation = DelReal(0);
			radius = DelReal(0);
			constexpr DelReal inf = std::numeric_limits<DelReal>::infinity();
			circumcenter = DelVector(inf, inf, inf);
			return;
		}
		DelVector da = v[0]->point - v[3]->point;
		DelVector db = v[1]->point - v[3]->point;
		DelVector dc = v[2]->point - v[3]->point;
		DelVector ca = v[0]->point - v[2]->point;
		DelVector cb = v[1]->point - v[2]->point;
		DelVector ba = v[0]->point - v[1]->point;

		Geometer::Circumsphere(v[0]->point, v[1]->point, v[2]->point, v[3]->point, &circumcenter, &radius);

		DelReal minEdgeLength = sqrt(std::min({ da.length2(), db.length2(), dc.length2(), ca.length2(), cb.length2(), ba.length2() }));
		DelReal relaxedLength = std::max(minEdgeLength, std::min({ v[0]->relaxedInsetionRadius, v[1]->relaxedInsetionRadius,
			v[2]->relaxedInsetionRadius, v[3]->relaxedInsetionRadius }));

		reRation = radius / relaxedLength;
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
		nodePool = new MemoryPool<TriVertexListNode>(64);
		vertPool = NULL;
		ghost = NULL;
		deadVerticesStack = NULL;
	}

	Vertex* TriMeshDataStructure::allocVertex(const DelVector &point, DelReal weight) {
		if (vertPool == NULL) vertPool = new MemoryArena<Vertex>(128);

		Vertex *newVert = NULL;
		if (deadVerticesStack == NULL) {
			newVert = vertPool->Alloc(1, point, weight, VertexType::Vertex_Facet);
			newVert->setLabel(labeler.getLabel());
			vertices.push_back(newVert);
		}
		else {
			newVert = deadVerticesStack;
			deadVerticesStack = newVert->getPointedVertex();
			int lable = newVert->getLabel();
			Construct(newVert, point, weight, VertexType::Vertex_Facet);
			newVert->setLabel(lable);
		}

		return newVert;
	}

	void TriMeshDataStructure::deallocVertex(Vertex *vert) {
		Assert(vertPool != NULL);
		TriVertexListNode *node = vert->getFaceLink();
		while (node) {
			TriVertexListNode *tmp = node->getNextNode();
			nodePool->Dealloc(node);
			node = tmp;
		}
		vert->type = Vertex_Undefined;
		vert->setVertexPointer(deadVerticesStack);
		deadVerticesStack = vert;
	}

	void TriMeshDataStructure::bindVolumeVertex(Vertex *vert) {
		Assert(matchVertexFlag(vert->type, VertexType(VertexType::Vertex_Volume | VertexType::Vertex_Facet)));
		Assert(vert->getFaceLink() == NULL);
		vertices.push_back(vert);
	}

	void TriMeshDataStructure::unbindVolumeVertex(Vertex *vert) {
		Assert(matchVertexFlag(vert->type, VertexType(VertexType::Vertex_Volume | VertexType::Vertex_Facet)));
		TriVertexListNode *node = vert->getFaceLink();
		while (node) {
			TriVertexListNode *tmp = node->getNextNode();
			nodePool->Dealloc(node);
			node = tmp;
		}
		vertices.erase(std::remove(vertices.begin(), vertices.end(), vert), vertices.end());
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
		if (getMarkedNode(u, v, &node))
			node->setMark();
	}

	void TriMeshDataStructure::unSetMark(Vertex *u, Vertex *v) {
		TriVertexListNode *node = NULL;
		if (getMarkedNode(u, v, &node)) 
			node->unSetMark();
	}

	bool TriMeshDataStructure::isMarked(Vertex *u, Vertex *v) const {
		TriVertexListNode *node = NULL;

		return getMarkedNode(u, v, &node) && node->isMarked();
	}

	bool TriMeshDataStructure::Adjacent(const Segment &s, Vertex **w, int *index) const {
		if (!matchVertexFlag(s.v[0]->getVertexType(), VertexType::Vertex_Facet) ||
			!matchVertexFlag(s.v[1]->getVertexType(), VertexType::Vertex_Facet))
			return false;

		TriVertexListNode *node = NULL;
		if (getAdjacentListNode(s.v[0], s.v[1], &node) && !node->isPreFaceDeleted()) {
			*w = node->getVertex();
			if (index) *index = node->getIndex();
			return true;
		}
		return false;
	}

	bool TriMeshDataStructure::adjacent2Vertex(Vertex *w, Triangle *f, int *index) const {
		if (!matchVertexFlag(w->type, VertexType::Vertex_Facet)) return false;

		TriVertexListNode *node = w->getFaceLink();
		TriVertexListNode *nextNode = node->getNextNode();

		while (nextNode != NULL) {
			if (!nextNode->isPreFaceDeleted()) {
				*f = Triangle(w, node->getVertex(), nextNode->getVertex());
				if (index) *index = nextNode->getIndex();
				return true;
			}
			node = nextNode;
			nextNode = node->getNextNode();
		}
		if (!node->isPreFaceDeleted()) {
			nextNode = w->getFaceLink();
			*f = Triangle(w, node->getVertex(), nextNode->getVertex());
			if (index) *index = nextNode->getIndex();
			return true;
		}

		return false;
	}

	int TriMeshDataStructure::getTriangleIndex(Vertex *a, Vertex *b, Vertex *c) const {
		TriVertexListNode *node = NULL;
		if (getAdjacentListNode(a, b, &node) && !node->isPreFaceDeleted() &&
			node->getVertex() == c)
			return node->getIndex();
		return -1;
	}

	TriMeshDataStructure::TriMeshConstIterator::TriMeshConstIterator(const std::vector<Vertex *>::const_iterator& iter,
		const std::vector<Vertex *>::const_iterator& end) : vertIter(iter), vertEnd(end){
		while (vertIter != vertEnd) {
			Vertex *center = *vertIter;
			if (matchVertexFlag(center->type, VertexType::Vertex_Facet) && 
				center->getFaceLink()) {
				parent = center->getFaceLink();
				child = parent->getNextNode();
				while (child != NULL) {
					if (!child->isPreFaceDeleted()) {
						Vertex *b = parent->getVertex();
						Vertex *c = child->getVertex();
						if (verticesOrderCheck(center, b, c)) {
							current = Triangle(center, b, c);
							return;
						}
					}

					parent = child;
					child = parent->getNextNode();
				}

				child = center->getFaceLink();
				if (!child->isPreFaceDeleted()) {
					Vertex *b = parent->getVertex();
					Vertex *c = child->getVertex();
					if (verticesOrderCheck(center, b, c)) {
						current = Triangle(center, b, c);
						return;
					}
				}
			}

			++vertIter;
		}
		parent = NULL; child = NULL;
	}

	void TriMeshDataStructure::TriMeshConstIterator::findNext() {
		if (vertIter == vertEnd) return;
		do {
			if (child != (*vertIter)->getFaceLink()) {
				parent = child;
				child = parent->getNextNode();
				if (child == NULL) child = (*vertIter)->getFaceLink();
			}
			else {
				do {
					++vertIter;
				} while (vertIter != vertEnd && 
					(!matchVertexFlag((*vertIter)->type, Vertex_Facet) || !((*vertIter)->getFaceLink())));

				if (vertIter != vertEnd) {
					parent = (*vertIter)->getFaceLink();
					child = parent->getNextNode();
				}
				else {
					current = Triangle(NULL, NULL, NULL);
					parent = NULL; child = NULL;
					return;
				}
			}

			if (!child->isPreFaceDeleted()) {
				Vertex *center = *vertIter;
				Vertex *b = parent->getVertex();
				Vertex *c = child->getVertex();
				if (verticesOrderCheck(center, b, c)) {
					current = Triangle(center, b, c);
					return;
				}
			}
		} while (true);
	}

	TriMeshDataStructure::TriMeshConstCirculator::TriMeshConstCirculator(Vertex *v) {
		vert = v;
		link = NULL;

		TriVertexListNode *node = vert->getFaceLink();
		if (node) {
			TriVertexListNode *nextNode = node->getNextNode();

			//find the first tiangle
			while (nextNode != NULL) {
				if (!nextNode->isPreFaceDeleted()) {
					link = nextNode;
					current = Triangle(vert, node->getVertex(), nextNode->getVertex());
					return;
				}
				node = nextNode;
				nextNode = node->getNextNode();
			}
			if (!node->isPreFaceDeleted()) {
				nextNode = vert->getFaceLink();
				link = nextNode;
				current = Triangle(vert, node->getVertex(), nextNode->getVertex());
			}
		}
	}

	void TriMeshDataStructure::TriMeshConstCirculator::findNext() {
		do {
			TriVertexListNode *parent = link;
			link = link->getNextNode();
			if (link == NULL) link = vert->getFaceLink();

			if (!link->isPreFaceDeleted()) {
				current = Triangle(vert, parent->getVertex(), link->getVertex());
				return;
			}
		} while (true);
	}

	std::vector<Triangle> TriMeshDataStructure::getTriangles(bool ghost) const {
		std::vector<Triangle> output;
		getTriangles(ghost, output);
		return output;
	}

	void TriMeshDataStructure::getTriangles(bool ghost, std::vector<Triangle>& triangles) const {
		TriVertexListNode *parent = NULL;
		TriVertexListNode *child = NULL;
		triangles.reserve(vertices.size() * 2);

		for (auto center : vertices) {
			if (matchVertexFlag(center->type, Vertex_Facet) 
				&& center->getFaceLink() && (ghost || !center->isGhost())) {
				parent = center->getFaceLink();
				child = parent->getNextNode();
				while (child != NULL) {
					if (!child->isPreFaceDeleted()) {
						Vertex *b = parent->getVertex();
						Vertex *c = child->getVertex();
						if (verticesOrderCheck(center, b, c))
							triangles.push_back(Triangle(center, b, c));
					}

					parent = child;
					child = child->getNextNode();
				}
				if (!center->getFaceLink()->isPreFaceDeleted()) {
					child = center->getFaceLink();
					Vertex *b = parent->getVertex();
					Vertex *c = child->getVertex();
					if (verticesOrderCheck(center, b, c))
						triangles.push_back(Triangle(center, b, c));
				}
			}
		}
	}

	void TriMeshDataStructure::addSegment(Vertex *a, Vertex *b) {
		if (a->getLabel() > b->getLabel()) std::swap(a, b);
		if (!matchVertexFlag(a->type, VertexType::Vertex_Facet)) return;

		TriVertexListNode *node = a->getFaceLink();
		while (node != NULL && node->getVertex() != b)
			node = node->getNextNode();
		if (node) node->setEdgeMark();
	}

	void TriMeshDataStructure::deleteSegment(Vertex *a, Vertex *b) {
		if (a->getLabel() > b->getLabel()) std::swap(a, b);
		if (!matchVertexFlag(a->type, VertexType::Vertex_Facet)) return;

		TriVertexListNode *node = a->getFaceLink();
		while (node != NULL && node->getVertex() != b)
			node = node->getNextNode();
		if (node) node->unSetEdgeMark();
	}

	bool TriMeshDataStructure::isSegment(const Segment& s) const {
		Vertex *a = s.v[0], *b = s.v[1];
		if (a->getLabel() > b->getLabel()) std::swap(a, b);
		if (!matchVertexFlag(a->type, VertexType::Vertex_Facet)) return false;

		TriVertexListNode *node = a->getFaceLink();
		while (node != NULL && node->getVertex() != b)
			node = node->getNextNode();
		if (node) return node->isEdgeMarked();

		return false;
	}

	bool TriMeshDataStructure::getAdjacentListNode(Vertex* u, Vertex* v, TriVertexListNode **w) const {
		TriVertexListNode *node = u->getFaceLink();
		bool found = false;
		while (node != NULL) {
			if (node->getVertex() == v) {
				TriVertexListNode *nextNode = node->getNextNode();
				if (nextNode == NULL)
					nextNode = u->getFaceLink();
				*w = nextNode;
				found = true;
				break;
			}
			node = node->getNextNode();
		}

		return found;
	}

	bool TriMeshDataStructure::getMarkedNode(Vertex* u, Vertex* v, TriVertexListNode **w) const {
		TriVertexListNode *node = u->getFaceLink();
		bool found = false;
		while (node != NULL) {
			if (node->getVertex() == v) {
				*w = node;
				found = true;
				break;
			}
			node = node->getNextNode();
		}

		return found;
	}

	void TriMeshDataStructure::insertToTopology(Vertex *a, Vertex *b, Vertex *c, int index) {
		if (!a->getFaceLink()) {
			TriVertexListNode *node = nodePool->Alloc();
			node->setVertex(b);
			node->setDeletedMark();
			TriVertexListNode *nextNode = nodePool->Alloc();
			nextNode->setVertex(c);
			nextNode->setIndex(index);
			node->setNextNode(nextNode);
			a->setFaceLink(node);
		}
		else {
			TriVertexListNode *head = a->getFaceLink();
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
					a->setFaceLink(head);
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
						a->setFaceLink(head);
					}
					else {
						foundNode[0]->unSetDeletedMark();
						newNode->setNextNode(foundNode[0]);
						newNode->setDeletedMark();
						parentFoundNode->setNextNode(newNode);
					}
				}
				else {
					Assert(foundNode[0]->getVertex() != c);
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
							a->setFaceLink(cNode);
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
							a->setFaceLink(head);
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
		TriVertexListNode *head = a->getFaceLink();
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
						a->setFaceLink(child);
						head = child;
					}
					//normal case
					else {
						parent->setDeletedMark();
						a->setFaceLink(parent);
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
					if (parent->getVertex() == b && head->getVertex() == c) {
						TriVertexListNode *next = head->getNextNode();
						//the node behide the head deleted
						if (next && next->isPreFaceDeleted()) {
							next->setDeletedMark();
							a->setFaceLink(next);
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
				a->setFaceLink(NULL);
			}
		}
	}

	void TriMeshDataStructure::Clear() {
		vertices.clear();
		nodePool->freeAll();
		if (vertPool) vertPool->freeAll();
		ghost = NULL;
		deadVerticesStack = NULL;
		labeler.restartLable();
	}

	TriMeshDataStructure::~TriMeshDataStructure() {
		delete nodePool;
		delete vertPool;
	}

	TetMeshDataStructure::TetMeshDataStructure() {
		nodePool = new MemoryPool<TetVertexListNode, std::alignment_of<TetVertexListNode>::value>(256);
		edgeNodePool = new MemoryPool<EdgeListNode, std::alignment_of<EdgeListNode>::value>();
		vertPool = new MemoryArena<Vertex>(256);
		pointerList = new ThreadUnsafeFreelist<sizeof(uintptr_t), 4 * sizeof(uintptr_t)>();
		ghost = NULL;
		deadVerticesStack = NULL;
	}

	Vertex* TetMeshDataStructure::allocVertex(const DelVector &point, DelReal weight, VertexType extraType) {
		constexpr size_t freeSegFacetPointerCount = 4, fixedSegFacetPointerCount = 3, facetPointerCount = 2;
		constexpr size_t freeSegPointerCount = 3, fixedSegPointerCount = 2;
		
		Vertex *newVert = NULL;
		if (deadVerticesStack == NULL) {
			newVert = vertPool->Alloc(1, point, weight, VertexType(Vertex_Volume | extraType));
			newVert->setLabel(labeler.getLabel());
			vertices.push_back(newVert);
		}
		else {
			newVert = deadVerticesStack;
			deadVerticesStack = newVert->getPointedVertex();
			int label = newVert->getLabel();
			Construct(newVert, point, weight, VertexType(Vertex_Volume | extraType));
			newVert->setLabel(label);
		}

		uintptr_t *pointers = NULL;
		if (matchVertexFlag(extraType, VertexType::Vertex_Facet)) {
			if (matchVertexFlag(extraType, VertexType::Vertex_Segment)) {
				if (matchVertexFlag(extraType, VertexType::Vertex_Free)) {
					pointers = pointerList->Alloc<uintptr_t>(freeSegFacetPointerCount);
					Initiation(pointers, freeSegFacetPointerCount);
				}
				else {
					pointers = pointerList->Alloc<uintptr_t>(fixedSegFacetPointerCount);
					Initiation(pointers, freeSegFacetPointerCount);
				}
			}
			else {
				pointers = pointerList->Alloc<uintptr_t>(facetPointerCount);
				Initiation(pointers, facetPointerCount);
			}
		}
		else if (matchVertexFlag(extraType, VertexType::Vertex_Segment)) {
			if (matchVertexFlag(extraType, VertexType::Vertex_Free)) {
				pointers = pointerList->Alloc<uintptr_t>(freeSegPointerCount);
				Initiation(pointers, freeSegPointerCount);
			}
			else {
				pointers = pointerList->Alloc<uintptr_t>(fixedSegPointerCount);
				Initiation(pointers, fixedSegPointerCount);
			}
		}

		newVert->setPointers(pointers);

		return newVert;
	}

	void TetMeshDataStructure::deallocVertex(Vertex *vert) {
		constexpr size_t freeSegFacetPointerCount = 4, fixedSegFacetPointerCount = 3, facetPointerCount = 2;
		constexpr size_t freeSegPointerCount = 3, fixedSegPointerCount = 2;
		VertexType type = vert->type;

		//clean edgelist
		EdgeListNode *linkHead = NULL;
		if (vert->hasEdgeList()) linkHead = vert->getEdgeList();
		while (linkHead) {
			TetVertexListNode *node = linkHead->getLink();
			while (node) {
				TetVertexListNode *tmp = node->getNextNode();
				nodePool->Dealloc(node);
				node = tmp;
			}
			EdgeListNode *tmp = linkHead->getNextNode();
			edgeNodePool->Dealloc(linkHead);
			linkHead = tmp;
		}

		//clean enforced edges
		if (matchVertexFlag(type, VertexType::Vertex_Segment)) {
			EdgeListNode *enfocedListHead = vert->getEnforcedEdgeList();
			while (enfocedListHead) {
				EdgeListNode *tmp = enfocedListHead->getNextNode();
				edgeNodePool->Dealloc(enfocedListHead);
				enfocedListHead = tmp;
			}
		}

		//clean pointers
		if (matchVertexFlag(type, VertexType::Vertex_Facet)) {
			if (matchVertexFlag(type, VertexType::Vertex_Segment)) {
				if (matchVertexFlag(type, VertexType::Vertex_Free))
					pointerList->Dealloc<uintptr_t>(vert->getPointers(), freeSegFacetPointerCount);
				else
					pointerList->Dealloc<uintptr_t>(vert->getPointers(), fixedSegFacetPointerCount);

				//erase it if vert in enforcedEdgeVertices
				if (vert->isEnforcedEdgeMarked())
					enforcedEdgeVertices.erase(std::remove(enforcedEdgeVertices.begin(), enforcedEdgeVertices.end(), vert),
						enforcedEdgeVertices.end());
			}
			else
				pointerList->Dealloc<uintptr_t>(vert->getPointers(), facetPointerCount);
		}
		else if (matchVertexFlag(type, VertexType::Vertex_Segment)) {
			if (matchVertexFlag(type, VertexType::Vertex_Free))
				pointerList->Dealloc<uintptr_t>(vert->getPointers(), freeSegPointerCount);
			else
				pointerList->Dealloc<uintptr_t>(vert->getPointers(), fixedSegPointerCount);

			//erase it if vert in enforcedEdgeVertices
			if (vert->isEnforcedEdgeMarked()) 
				enforcedEdgeVertices.erase(std::remove(enforcedEdgeVertices.begin(), enforcedEdgeVertices.end(), vert),
					enforcedEdgeVertices.end());
		}

		vert->type = VertexType::Vertex_Undefined;
		vert->setVertexPointer(deadVerticesStack);
		deadVerticesStack = vert;
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
		Triangle f(u, v, w, true);
		TetVertexListNode *node = NULL;

		int condition = getMarkedNode(f, &node);
		if (condition > 0) {
			if (condition == 1) node->setForwardMark();
			else node->setBackwardMark();
		}
	}

	void TetMeshDataStructure::unSetMark(Vertex *u, Vertex *v, Vertex *w) {
		Triangle f(u, v, w, true);
		TetVertexListNode *node = NULL;

		int condition = getMarkedNode(f, &node);
		if (condition > 0) {
			if (condition == 1) node->unSetForwardMark();
			else node->unSetBackwardMark();
		}
	}

	bool TetMeshDataStructure::isMarked(Vertex *u, Vertex *v, Vertex *w) const {
		Triangle f(u, v, w, true);
		TetVertexListNode *node = NULL;

		int condition = getMarkedNode(f, &node);
		if (condition > 0) {
			if (condition == 1) return node->isForwardMarked();
			else return node->isBackwardMarked();
		}
		return false;
	}

	bool TetMeshDataStructure::Adjacent(const Triangle &f, Vertex **z) const {
		TetVertexListNode *node = NULL;
		if (getAdjacentListNode(f, &node)) {
			*z = node->getVertex();
			return true;
		}
		return false;
	}

	bool TetMeshDataStructure::getAdjacentListNode(const Triangle& f, TetVertexListNode **z) const {
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
			if (x->hasEdgeList()) {
				EdgeListNode *linkHead = x->getEdgeList();
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
			if (w->hasEdgeList()) {
				EdgeListNode *linkHead = w->getEdgeList();
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


	int TetMeshDataStructure::getMarkedNode(const Triangle& f, TetVertexListNode **marked) const {
		int found = 0;
		bool ab = parityCheck(f.v[0], f.v[1]);
		bool ac = parityCheck(f.v[0], f.v[2]);

		if (ab)
			found = getMarkedNode(f.v[0], f.v[1], f.v[2], marked);
		else if (ac)
			found = getMarkedNode(f.v[2], f.v[0], f.v[1], marked);
		else
			found = getMarkedNode(f.v[1], f.v[2], f.v[0], marked);

		return found;
	}

	int TetMeshDataStructure::getMarkedNode(Vertex *w, Vertex *x, Vertex *y, TetVertexListNode **marked) const {
		int found = 0;
		bool swap = edgeOrderCheck(x, w);
		if (swap) std::swap(x, w);
		if (w->hasEdgeList()) {
			EdgeListNode *linkHead = w->getEdgeList();
			while (linkHead != NULL && linkHead->getEndVertex() != x) {
				linkHead = linkHead->getNextNode();
			}
			if (linkHead != NULL) {
				TetVertexListNode *node = linkHead->getLink();
				while (node != NULL && found == 0) {
					if (node->getVertex() == y) {
						*marked = node;
						found = swap ? 1 : 2;
					}
					node = node->getNextNode();
				}
			}
		}

		return found;
	}

	bool TetMeshDataStructure::adjacent2Vertex(Vertex *w, Tetrahedron *t) const {
		bool found = false;
		EdgeListNode *linkHead = NULL;
		bool hasList = w->hasEdgeList();
		if (!hasList) {
			Vertex *end = w->getSupplyVertex();
			if (end == NULL || !end->hasEdgeList()) return adjacent2VertexSlow(w, t);
			linkHead = end->getEdgeList();
			bool foundHead = false;
			while (linkHead != NULL && !foundHead) {
				if (linkHead->getEndVertex() == w)
					foundHead = true;
				else
					linkHead = linkHead->getNextNode();
			}
		}
		else
			linkHead = w->getEdgeList();

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
					*t = Tetrahedron(w, w->getSupplyVertex(), loop->getVertex(), parentLoop->getVertex());
				return true;
			}
		}
		else {
			if (!hasList) {
				Vertex *oppo = w->getSupplyVertex();
				Assert(oppo->hasEdgeList());
				linkHead = oppo->getEdgeList();
				while (linkHead) {
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
							return true;
						}
						else if (!childLoop->isPreFaceDeleted()) {
							*t = Tetrahedron(w, childLoop->getVertex(), oppo, linkHead->getEndVertex());
							return true;
						}
					}
					linkHead = linkHead->getNextNode();
				}
			}
		}

		return adjacent2VertexSlow(w, t);
	}

	void TetMeshDataStructure::addSegment(Vertex *a, Vertex *b) {
		if (!edgeOrderCheck(a, b)) std::swap(a, b);
		EdgeListNode *node = a->getEnforcedEdgeList();
		while (node != NULL && node->getEndVertex() != b)
			node = node->getNextNode();

		if (node == NULL) {
			EdgeListNode *oldNode = a->getEnforcedEdgeList();
			EdgeListNode *newNode = edgeNodePool->Alloc();
			newNode->setEndVertex(b);
			newNode->setNextNode(oldNode);
			newNode->setOtherVertex(NULL);
			a->setEnforcedEdgeList(newNode);
			if (oldNode == NULL && !a->isEnforcedEdgeMarked()) {
				a->setEnforcedEdgeMark();
				enforcedEdgeVertices.push_back(a);
			}
		}
	}

	void TetMeshDataStructure::deleteSegment(Vertex *a, Vertex *b) {
		if (!edgeOrderCheck(a, b)) std::swap(a, b);
		EdgeListNode *node = a->getEnforcedEdgeList();
		EdgeListNode *parentNode = NULL;
		while (node != NULL && node->getEndVertex() != b) {
			parentNode = node;
			node = node->getNextNode();
		}

		if (node) {
			if (parentNode)
				parentNode->setNextNode(node->getNextNode());
			else {
				EdgeListNode *next = node->getNextNode();
				a->setEnforcedEdgeList(next);
			}

			edgeNodePool->Dealloc(node);
		}
	}

	bool TetMeshDataStructure::adjacent2SegmentFast(const Segment &s, Tetrahedron *t) const {
		Vertex *a = s.v[0], *b = s.v[1];
		if (parityCheck(a, b)) {
			bool ordered = edgeOrderCheck(a, b);
			if (!ordered) std::swap(a, b);
			if (a->hasEdgeList()) {
				EdgeListNode *linkHead = a->getEdgeList();
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
			bool ordered = edgeOrderCheck(a, b);
			if (!ordered) std::swap(a, b);
			EdgeListNode *node = a->getEnforcedEdgeList();
			while (node != NULL && node->getEndVertex() != b) 
				node = node->getNextNode();

			if (node) {
				Vertex *c = node->getOtherVertex(), *d = NULL;
				if (c != NULL && Adjacent(Triangle(b, a, c), &d)) {
					if (ordered) *t = Tetrahedron(a, b, c, d);
					else  *t = Tetrahedron(b, a, d, c);
					return true;
				}
			}
		}

		return false;
	}

	void TetMeshDataStructure::setMark(Vertex *a, Vertex *b) {
		if (!matchVertexFlag(a->type, VertexType::Vertex_Segment) ||
			!matchVertexFlag(b->type, VertexType::Vertex_Segment))
			return;

		if (!edgeOrderCheck(a, b)) std::swap(a, b);
		EdgeListNode *node = a->getEnforcedEdgeList();
		while (node != NULL && node->getEndVertex() != b)
			node = node->getNextNode();

		if (node) node->setMark();
	}

	void TetMeshDataStructure::unSetMark(Vertex *a, Vertex *b) {
		if (!matchVertexFlag(a->type, VertexType::Vertex_Segment) ||
			!matchVertexFlag(b->type, VertexType::Vertex_Segment))
			return;

		if (!edgeOrderCheck(a, b)) std::swap(a, b);
		EdgeListNode *node = a->getEnforcedEdgeList();
		while (node != NULL && node->getEndVertex() != b)
			node = node->getNextNode();

		if (node) node->unSetMark();
	}

	bool TetMeshDataStructure::isMarked(Vertex *a, Vertex *b) const {
		bool ret = false;
		if (!matchVertexFlag(a->type, VertexType::Vertex_Segment) ||
			!matchVertexFlag(b->type, VertexType::Vertex_Segment))
			return ret;

		if (!edgeOrderCheck(a, b)) std::swap(a, b);
		EdgeListNode *node = a->getEnforcedEdgeList();
		while (node != NULL && node->getEndVertex() != b)
			node = node->getNextNode();

		if (node) ret = node->isMarked();

		return ret;
	}

	bool TetMeshDataStructure::testAndMark(Vertex *a, Vertex *b) {
		bool ret = false;
		if (!matchVertexFlag(a->type, VertexType::Vertex_Segment) ||
			!matchVertexFlag(b->type, VertexType::Vertex_Segment))
			return ret;

		if (!edgeOrderCheck(a, b)) std::swap(a, b);
		EdgeListNode *node = a->getEnforcedEdgeList();
		while (node != NULL && node->getEndVertex() != b)
			node = node->getNextNode();

		if (node) {
			ret = node->isMarked();
			node->setMark();
		}

		return ret;
	}

	bool TetMeshDataStructure::adjacent2VertexSlow(Vertex *w, Tetrahedron *t) const {
		for (auto vert : vertices) {
			if (matchVertexFlag(vert->type, VertexType::Vertex_Volume) 
				&& vert->hasEdgeList()) {
				EdgeListNode *linkHead = vert->getEdgeList();
				while (linkHead) {
					Vertex *endVert = linkHead->getEndVertex();
					TetVertexListNode *head = linkHead->getLink();
					TetVertexListNode *parentNode = head;
					TetVertexListNode *node = parentNode->getNextNode();
					while (node) {
						if (!node->isPreFaceDeleted()) {
							Vertex *c = parentNode->getVertex();
							Vertex *d = node->getVertex();
							if (endVert == w) {
								w->setSupplyVertexPointer(vert);
								*t = Tetrahedron(w, vert, d, c);
								return true;
							}
							else if (c == w) {
								w->setSupplyVertexPointer(vert);
								*t = Tetrahedron(w, d, vert, endVert);
								return true;
							}
							else if (d == w) {
								w->setSupplyVertexPointer(vert);
								*t = Tetrahedron(w, c, endVert, vert);
								return true;
							}
						}
						parentNode = node;
						node = node->getNextNode();
					}
					if (!head->isPreFaceDeleted()) {
						Vertex *c = parentNode->getVertex();
						Vertex *d = head->getVertex();
						if (endVert == w) {
							w->setSupplyVertexPointer(vert);
							*t = Tetrahedron(w, vert, d, c);
							return true;
						}
						if (c == w) {
							w->setSupplyVertexPointer(vert);
							*t = Tetrahedron(w, d, vert, endVert);
							return true;
						}
						else if (d == w) {
							w->setSupplyVertexPointer(vert);
							*t = Tetrahedron(w, c, endVert, vert);
							return true;
						}
					}
					linkHead = linkHead->getNextNode();
				}
			}
		}

		return false;
	}

	bool TetMeshDataStructure::adjacent2SegmentSlow(const Segment &s, Tetrahedron *t) const {
		Vertex *a = s.v[0], *b = s.v[1];
		for (auto vert : vertices) {
			if (matchVertexFlag(vert->type, VertexType::Vertex_Volume)
				&& vert->hasEdgeList()) {	
				EdgeListNode *linkHead = vert->getEdgeList();
				while (linkHead) {
					Vertex *endVert = linkHead->getEndVertex();

					int condition = 0;
					if (vert == a) condition = 1;
					else if (vert == b) condition = 2;
					else if (endVert == a) condition = 3;
					else if (endVert == b) condition = 4;

					if (condition != 0) {
						TetVertexListNode *head = linkHead->getLink();
						TetVertexListNode *parentNode = head;
						TetVertexListNode *node = parentNode->getNextNode();
						do {
							if (node == NULL) node = head;
							if (!node->isPreFaceDeleted()) {
								Vertex *c = parentNode->getVertex();
								Vertex *d = node->getVertex();

								//check all the case
								switch (condition) {
								case 1:
									if (c == b) {
										*t = Tetrahedron(a, b, d, endVert);
										return true;
									}
									else if (d == b) {
										*t = Tetrahedron(a, b, endVert, c);
										return true;
									}
									break;
								case 2:
									if (c == a) {
										*t = Tetrahedron(a, b, endVert, d);
										return true;
									}
									else if (d == a) {
										*t = Tetrahedron(a, b, c, endVert);
										return true;
									}
									break;
								case 3:
									if (c == b) {
										*t = Tetrahedron(a, b, vert, d);
										return true;
									}
									else if (d == b) {
										*t = Tetrahedron(a, b, c, vert);
										return true;
									}
									break;
								case 4:
									if (c == a) {
										*t = Tetrahedron(a, b, d, vert);
										return true;
									}
									else if (d == a) {
										*t = Tetrahedron(a, b, vert, c);
										return true;
									}
									break;
								default:
									break;
								}
							}
							parentNode = node;
							node = node->getNextNode();
						} while (parentNode != head);
					}
					linkHead = linkHead->getNextNode();
				}
			}
		}
		return false;
	}

	TetMeshDataStructure::TetMeshConstIterator::TetMeshConstIterator(const std::vector<Vertex *>::const_iterator& iter,
		const std::vector<Vertex *>::const_iterator& end) : vertIter(iter), vertEnd(end){
		linkHead = NULL;
		while (vertIter != vertEnd) {
			Vertex *ori = *vertIter;
			if (matchVertexFlag(ori->type, VertexType::Vertex_Volume) && ori->hasEdgeList()) 
				linkHead = ori->getEdgeList();

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
			++vertIter;
		}

		parent = NULL; child = NULL;
	}

	void TetMeshDataStructure::TetMeshConstIterator::findNext() {
		if (vertIter == vertEnd) return;
		do {
			if (child != linkHead->getLink()) {
				parent = child;
				child = parent->getNextNode();
				if (!child) child = linkHead->getLink();
			}
			else {
				linkHead = linkHead->getNextNode();
				if (!linkHead) {
					do {
						++vertIter;
					} while (vertIter != vertEnd && 
						(!matchVertexFlag((*vertIter)->type, VertexType::Vertex_Volume) || !((*vertIter)->hasEdgeList())));

					if (vertIter != vertEnd)
						linkHead = (*vertIter)->getEdgeList();
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
				Vertex *ori = *vertIter;
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

	TetMeshDataStructure::SegmentConstIterator::SegmentConstIterator(const std::vector<Vertex *>::const_iterator &iter, 
		const std::vector<Vertex *>::const_iterator &end): edgeVertIter(iter), edgeVertEnd(end) {
		edgeNode = NULL;
		while (edgeVertIter != edgeVertEnd) {
			Vertex *ori = *edgeVertIter;
			edgeNode = ori->getEnforcedEdgeList();
			if (edgeNode) {
				segment = Segment(ori, edgeNode->getEndVertex());
				break;
			}
			++edgeVertIter;
		}
	}

	void TetMeshDataStructure::SegmentConstIterator::findNext() {
		if (edgeVertIter == edgeVertEnd) return;

		edgeNode = edgeNode->getNextNode();
		if (edgeNode) segment.v[1] = edgeNode->getEndVertex();
		else {
			do {
				++edgeVertIter;
			} while (edgeVertIter != edgeVertEnd &&
				 !((*edgeVertIter)->getEnforcedEdgeList()));

			if (edgeVertIter != edgeVertEnd) {
				Vertex *ori = *edgeVertIter;
				edgeNode = ori->getEnforcedEdgeList();
				segment = Segment(ori, edgeNode->getEndVertex());
			}
			else {
				edgeNode = NULL;
				segment = Segment(NULL, NULL);
			}
		}
	}

	TetMeshDataStructure::TetMeshFacetConstCirculator::TetMeshFacetConstCirculator(const Segment& s, const TetMeshDataStructure *dataStructrue) {
		ori = s.v[0]; dest = s.v[1];
		ds = dataStructrue;
		linkListHead = NULL;
		listNode = NULL;
		initAdjacentVert = NULL;
		movingVert = NULL;
		reverse = false;

		bool flip = false;
		if (!ds->edgeOrderCheck(ori, dest)) {
			flip = true;
			std::swap(ori, dest);
		}

		//find the first facet
		if (ds->fastSegmentQueryCheck(s)) {
			if (ori->hasEdgeList()) {
				EdgeListNode *linkHead = ori->getEdgeList();
				while (linkHead != NULL && linkHead->getEndVertex() != dest)
					linkHead = linkHead->getNextNode();

				if (linkHead != NULL) {
					linkListHead = linkHead;
					listNode = linkHead->getLink();
					if (!flip) current = Triangle(ori, dest, listNode->getVertex());
					else current = Triangle(dest, ori, listNode->getVertex());
				}
			}
		}
		else {
			EdgeListNode *node = ori->getEnforcedEdgeList();
			while (node != NULL && node->getEndVertex() != dest)
				node = node->getNextNode();

			if (node) {
				Vertex *c = node->getOtherVertex();
				if (c && ds->Contain(Triangle(dest, ori, c))) {
					initAdjacentVert = c;
					movingVert = c;
					if (!flip) current = Triangle(ori, dest, c);
					else current = Triangle(dest, ori, c);
				}
			}
		}
	}

	void TetMeshDataStructure::TetMeshFacetConstCirculator::findNext() {
		if (linkListHead) {
			listNode = listNode->getNextNode();
			if (listNode == NULL) listNode = linkListHead->getLink();
			current.v[2] = listNode->getVertex();
		}
		else {
			Vertex *c = NULL;
			if (ds->Adjacent(Triangle(dest, ori, movingVert), &c))
				movingVert = c;
			else {
				reverse = !reverse;
				if (reverse) {
					if (ds->Adjacent(Triangle(ori, dest, initAdjacentVert), &movingVert))
						std::swap(dest, ori);
					else {
						movingVert = initAdjacentVert;
						reverse = !reverse;
					}
				}
				else {
					movingVert = initAdjacentVert;
					std::swap(dest, ori);
				}
			}

			current.v[2] = movingVert;
		}
	}

	void TetMeshDataStructure::getTetrahedrons(bool ghost, std::vector<Tetrahedron>& tets) const {
		for (auto vert : vertices) {
			if (matchVertexFlag(vert->type, VertexType::Vertex_Volume) && vert->hasEdgeList()) {
				EdgeListNode *linkHead = vert->getEdgeList();
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
	}

	std::vector<Tetrahedron> TetMeshDataStructure::getTetrahedrons(bool ghost) const {
		std::vector<Tetrahedron> tets;
		getTetrahedrons(ghost, tets);
		return tets;
	}

	void TetMeshDataStructure::Clear() {
		vertices.clear();
		nodePool->freeAll();
		edgeNodePool->freeAll();
		vertPool->freeAll();
		pointerList->freeAll();
		labeler.restartLable();
		ghost = NULL;
		deadVerticesStack = NULL;
	}

	void TetMeshDataStructure::insertToTopology(const Segment& s, Vertex *mayC, Vertex *mayD) {
		Vertex *a = s.v[0], *b = s.v[1], *c = mayC, *d = mayD;
		if (!edgeOrderCheck(a, b)) {
			std::swap(a, b);
			std::swap(c, d);
		}

		if (!a->hasEdgeList()) {
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
			a->setEdgeList(newNode);
		}
		else {
			EdgeListNode *parentLinkHead = NULL;
			EdgeListNode *linkHead = a->getEdgeList();
			bool foundHead = false;
			int distance = 0;
			while (linkHead != NULL && !foundHead) {
				if (linkHead->getEndVertex() == b)
					foundHead = true;
				else {
					parentLinkHead = linkHead;
					linkHead = linkHead->getNextNode();
					distance += 1;
				}
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
				newNode->setNextNode(a->getEdgeList());
				a->setEdgeList(newNode);
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
					Severe("Unexpted case in TetMeshDataStructure::insertToTopology");
					break;
				}

				//move the linkHead to the start of the edge list of the vertex
				constexpr int cacheLimit = 3;
				if (distance > cacheLimit) {
					parentLinkHead->setNextNode(linkHead->getNextNode());
					linkHead->setNextNode(a->getEdgeList());
					a->setEdgeList(linkHead);
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

		Assert(a->hasEdgeList());
		EdgeListNode *parentLinkHead = NULL;
		EdgeListNode *linkHead = a->getEdgeList();
		bool foundHead = false;
		int distance = 0;
		while (linkHead != NULL && !foundHead) {
			if (linkHead->getEndVertex() == b)
				foundHead = true;
			else {
				parentLinkHead = linkHead;
				linkHead = linkHead->getNextNode();
				distance += 1;
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
			constexpr int cacheLimit = 3;
			//clean the deleted node
			if (head->getNextNode() == NULL) {
				nodePool->Dealloc(head);
				linkHead->setLink(NULL);
				if (parentLinkHead != NULL)
					parentLinkHead->setNextNode(linkHead->getNextNode());
				else
					a->setEdgeList(linkHead->getNextNode());

				edgeNodePool->Dealloc(linkHead);
			}
			//move the linkHead to the start of the edge list of the vertex
			else if (distance > cacheLimit) {
				parentLinkHead->setNextNode(linkHead->getNextNode());
				linkHead->setNextNode(a->getEdgeList());
				a->setEdgeList(linkHead);
			}
		}
	}

	void TetMeshDataStructure::addSupplyVerts(Vertex *a, Vertex *b, Vertex *c, Vertex *d, int mode) {
		switch (mode) {
		case 0:
			if (!a->isGhost() && !a->hasEdgeList()) a->setSupplyVertexPointer(b);
			if (!b->isGhost() && !b->hasEdgeList()) b->setSupplyVertexPointer(a);
			if (!c->isGhost() && !c->hasEdgeList()) c->setSupplyVertexPointer(d);
			if (!d->isGhost() && !d->hasEdgeList()) d->setSupplyVertexPointer(c);
			break;
		case 1:
			if (!a->isGhost() && !a->hasEdgeList()) a->setSupplyVertexPointer(b);
			if (!b->isGhost() && !b->hasEdgeList()) b->setSupplyVertexPointer(c);
			if (!c->isGhost() && !c->hasEdgeList()) c->setSupplyVertexPointer(a);
			if (!d->isGhost() && !d->hasEdgeList()) d->setSupplyVertexPointer(edgeOrderCheck(a, b) ? a : b);
			break;
		case 2:
			if (!a->isGhost() && !a->hasEdgeList()) a->setSupplyVertexPointer(b);
			if (!b->isGhost() && !b->hasEdgeList()) b->setSupplyVertexPointer(c);
			if (!c->isGhost() && !c->hasEdgeList()) c->setSupplyVertexPointer(d);
			if (!d->isGhost() && !d->hasEdgeList()) d->setSupplyVertexPointer(a);
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
			if (!edgeOrderCheck(a, b)) {
				std::swap(a, b);
				std::swap(c, d);
			}

			EdgeListNode *node = a->getEnforcedEdgeList();
			while (node != NULL && node->getEndVertex() != b)
				node = node->getNextNode();

			if (node) node->setOtherVertex(c);
		}
	}

	TetMeshDataStructure::~TetMeshDataStructure() {
		delete nodePool;
		delete edgeNodePool;
		delete vertPool;
		delete pointerList;
	}
}
