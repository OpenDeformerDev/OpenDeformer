#include "stdafx.h"
#include "DelMesher.h"
#include <set>

namespace ODER{

Predicator<REAL> DelTriangulator::predicator;
Labeler Vertex::labeler;

void DelTriangulator::generateSubPolygons(std::vector<Vertex *> &vertices){
	topology.clear();
	nodePool->freeAll();

	Vertex *v0, *v1, *v2, *v3;

#ifdef ODER_DEBUG
	if (vertices.size() < 3)
		Severe("a subpolygon contains less than 3 vertices");
#endif

	int size = vertices.size();

	v0 = vertices[0]; v1 = vertices[1];
	for (int i = 2; i < size; i++){ 
		v2 = vertices[i];
		if (!predicator.coLine(v2->vert, v0->vert, v1->vert)){
			std::swap(vertices[2], vertices[i]);
			break;
		}
	}

	DelVector a = v1->vert - v0->vert;
	DelVector b = v2->vert - v0->vert;
	n = Normalize(a % b);

	Face f(v0, v1, v2);
	//start
	addTriangle(v0, v1, v2);
	addTriangle(v1, v0, ghost);
	addTriangle(v2, v1, ghost);
	addTriangle(v0, v2, ghost);
	for (int i = 3; i < size; i++){
		v3 = vertices[i];
		f = findPosition(v3, f, n);
		v0 = f.v[0]; v1 = f.v[1]; v2 = f.v[2];
		deleteTriangle(v0, v1, v2);
		digCavity(v3, Segment(v1, v0), n, &f);
		digCavity(v3, Segment(v0, v2), n);
		digCavity(v3, Segment(v2, v1), n);
	}
}

void DelTriangulator::outPut(DelMesher *mesher){
	std::set<Face, face_compare> output;
	VertexListNode *head = NULL;
	VertexListNode *parent = NULL;
	VertexListNode *child = NULL;
	
	for (auto entry : topology){
		Vertex *center = entry.first;
		if (center != ghost && entry.second){
			parent = entry.second;
			child = parent->getNextNode();
			while (child != NULL){
				Vertex *b = parent->getVertex();
				Vertex *c = child->getVertex();
				if (!child->isPreFaceDeleted()){
					Face f(center, b, c, true);
					if (f.v[0] != ghost)
					    output.insert(f);
				}
				parent = child;
				child = child->getNextNode();
			}
			if (!entry.second->isPreFaceDeleted()){
				Vertex *b = parent->getVertex();
				Vertex *c = entry.second->getVertex();
				Face f(center, b, c, true);
				if (f.v[0] != ghost)
					output.insert(f);
			}
		}
	}
	for (auto f : output){
		mesher->addTriangle(f.v[0], f.v[1], f.v[2]);
		mesher->mayEncroachedFaces.push_back(f);
	}
}

void DelTriangulator::addTriangle(Vertex *a, Vertex *b, Vertex *c){
	/*Face f(a, b, c, true);
	if (f.v[0] != ghost)
		polygons.insert(f);
	svHash[Segment(a, b)] = c;
	svHash[Segment(b, c)] = a;
	svHash[Segment(c, a)] = b;*/

	insertToTopology(a, b, c);
	insertToTopology(b, c, a);
	insertToTopology(c, a, b);
}

void DelTriangulator::deleteTriangle(Vertex *a, Vertex *b, Vertex *c){
	/*Face f(a, b, c, true);
	if (f.v[0] != ghost)
	    polygons.erase(f);
	svHash.erase(Segment(a, b));
	svHash.erase(Segment(b, c));
	svHash.erase(Segment(c, a));*/

	removeFromTopology(a, b, c);
	removeFromTopology(b, c, a);
	removeFromTopology(c, a, b);
}

Face DelTriangulator::findPosition(Vertex *u, const Face &f, const DelVector &n) const{
	size_t maxIter = std::max(topology.size(), size_t(128));
	Vertex *a = f.v[0], *b = f.v[1], *c = f.v[2];
	//ghost triangle
	if (a == ghost){
		if (predicator.inHalfSpace2D(u->vert, b->vert, c->vert, n))
			return f;
		Vertex *x = ghost;
		Adjacent(Segment(c, b), &x);
		//Face(x, c, b)
		a = x; std::swap(b, c);
	}

	Assert(a != ghost);
	//normal triangle
	while (maxIter--){
		//caution:because NULL->ghost,can't use NULL,
		//so it's might case memory corrupt above
		Vertex *x;
		bool ori0 = predicator.orient2d(b->vert, a->vert, u->vert, n) > 0.0;
		bool ori1 = predicator.orient2d(a->vert, c->vert, u->vert, n) > 0.0;

		int condition = ori0 + (ori1 << 1);
		switch (condition){
		case 0:
			if (predicator.orient2d(c->vert, b->vert, u->vert, n) > 0.0){
				Adjacent(Segment(c, b), &x);
				if (x == ghost)
					return Face(ghost, c, b);
				a = x; std::swap(b, c); //Face(x, c, b)
			}
			else
				return Face(a, b, c, true);
			break;
		case 1:
			Adjacent(Segment(b, a), &x);
			if (x == ghost)
				return Face(ghost, b, a);
			c = a; a = x; //Face(x, b, a)
			break;
		case 2:
			Adjacent(Segment(a, c), &x);
			if (x == ghost)
				return Face(ghost, a, c);
			b = a; a = x; //Face(x, a, c)
			break;
		case 3:
			if (Randomnation(2)){
				Adjacent(Segment(b, a), &x);
				if (x == ghost)
					return Face(ghost, b, a);
				c = a; a = x; //Face(x, b, a)
			}
			else{
				Adjacent(Segment(a, c), &x);
				if (x == ghost)
					return Face(ghost, a, c);
				b = a; a = x; //Face(x, a, c)
			}
			break;
		default:
			Severe("Unexpected Case in DelTriangulator::findPosition");
		}
	}
	Severe("Do not find Face in DelTriangulator::findPosition");
	//should never return this value
	return Face(NULL, NULL, NULL);
}

bool DelTriangulator::Adjacent(const Segment &s, Vertex **w) const{
	//auto wp = svHash.find(s);
	//if (wp == svHash.end())
		//return false;
	//*w = wp->second;
	//return true;

	auto entry = topology.find(s.v[0]);
	if (entry == topology.end())
		return false;
	VertexListNode *node = entry->second;
	bool found = false;
	while (node != NULL){
		if (node->getVertex() == s.v[1]){
			VertexListNode *nextNode = node->getNextNode();
			if (nextNode == NULL)
				nextNode = entry->second;
			if (!nextNode->isPreFaceDeleted()){
				*w = nextNode->getVertex();
				found = true;
			}
			break;
		}
		node = node->getNextNode();
	}
	return found;
}

void DelTriangulator::digCavity(Vertex *u, const Segment &s, const DelVector& n, Face *rf){
	Vertex *x;
	if (!Adjacent(s, &x))
		return;

	//Vertex *w = s.v[0], *v = s.v[1], *x;
	Vertex *verts[3];
	verts[0] = s.v[0];
	verts[1] = s.v[1];
	verts[2] = x;


	bool deleted = false;
	int i = -1;
	if ((i = findGhost(verts, 3)) >= 0){
		//ghost triangle
		Vertex *w = verts[NEXT_F(i)];
		Vertex *v = verts[NEXT_F(i + 1)];
		deleted = predicator.inHalfSpace2D(u->vert, w->vert, v->vert, n);
	}
	else{
		deleted = predicator.inOrthoCirclePerturbed(u->vert, u->weight, verts[1]->vert, verts[1]->weight, 
			                               verts[0]->vert, verts[0]->weight, verts[2]->vert, verts[2]->weight ,n) > 0.0;
	}
	if (deleted){
		deleteTriangle(verts[0], verts[1], verts[2]);
		digCavity(u, Segment(verts[2], verts[1]), n, rf);
		digCavity(u, Segment(verts[0], verts[2]), n, rf);
	}
	else{
		addTriangle(u, verts[1], verts[0]);
		if (rf)
			*rf = Face(u, verts[1], verts[0], true);
	}
}

void DelTriangulator::insertToTopology(Vertex *a, Vertex *b, Vertex *c){
	auto entry = topology.find(a);
	if (entry == topology.end()){
		VertexListNode *node = nodePool->Alloc();
		node->setVertex(b);
		node->setDeletedMark();
		VertexListNode *nextNode = nodePool->Alloc();
		nextNode->setVertex(c);
		node->setNextNode(nextNode);
		topology[a] = node;
	}
	else{
		VertexListNode *head = entry->second;
		VertexListNode *parent = head;
		VertexListNode *child = parent->getNextNode();
		VertexListNode *foundNode[2];
		VertexListNode *parentFoundNode;
		int found = 0;
		if (head->getVertex() == b || head->getVertex() == c){
			foundNode[found++] = head;
		}
		while (child != NULL && found < 2){
			Vertex *foundVert = child->getVertex();
			if (foundVert == b)
				foundNode[found++] = child;
			else if (foundVert == c){
				parentFoundNode = parent;
				foundNode[found++] = child;
			}
			parent = child;
			child = child->getNextNode();
		}

		switch (found){
		case 0:
		{
			VertexListNode *newNode = nodePool->Alloc();
			newNode->setVertex(b);
			newNode->setDeletedMark();
			VertexListNode *nextNode = nodePool->Alloc();
			nextNode->setVertex(c);
			newNode->setNextNode(nextNode);

			if (head->isPreFaceDeleted()){
				nextNode->setNextNode(head);
				head = newNode;
				entry->second = head;
			}
			else{
				parent = head;
				child = parent->getNextNode();
				while (child && !child->isPreFaceDeleted()){
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
			if (foundNode[0]->isPreFaceDeleted()){
				newNode->setVertex(b);
				if (foundNode[0] == head){
					head->unSetDeletedMark();
					newNode->setDeletedMark();
					newNode->setNextNode(head);
					head = newNode;
					entry->second = head;
				}
				else{
					foundNode[0]->unSetDeletedMark();
					newNode->setNextNode(foundNode[0]);
					newNode->setDeletedMark();
					parentFoundNode->setNextNode(newNode);
				}
			}
			else{
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
			if (unoredered){
				bNode = foundNode[1];  cNode = foundNode[0];
			}
			else{
				bNode = foundNode[0];  cNode = foundNode[1];
			}
			if (bNode->getNextNode() == cNode || (cNode == head && bNode->getNextNode() == NULL))
				cNode->unSetDeletedMark();
			else{
				//find the lost fan
				VertexListNode *start = bNode->getNextNode();
				VertexListNode *end = parentFoundNode;

				if (unoredered){
					if (head == cNode){
						VertexListNode *loop = start;
						while (loop != NULL){
							end = loop;
							loop = loop->getNextNode();
						}
						bNode->setNextNode(NULL);
					}
					else{
						//the lost fan incule tail->head
						VertexListNode *loop = start;
						if (loop == NULL){
							start = head;
							loop = start;
						}
						else{
							while (loop != NULL){
								end = loop;
								loop = loop->getNextNode();
							}
							end->setNextNode(head);
							loop = head;
						}

						while (loop != cNode){
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
				while (gapEnd && !gapEnd->isPreFaceDeleted()){
					gapStart = gapEnd;
					gapEnd = gapEnd->getNextNode();
				}
				//insert it
				if (gapEnd != NULL){
					if (gapStart != NULL){
						gapStart->setNextNode(start);
						end->setNextNode(gapEnd);
					}
					else{
						end->setNextNode(head);
						head = start;
						entry->second = head;
					}
				}
				//the lost fan lost
				else{
					VertexListNode *loop = start;
					while (loop && loop != end){
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

void DelTriangulator::removeFromTopology(Vertex *a, Vertex *b, Vertex *c){
	auto entry = topology.find(a);
	Assert(entry != topology.end());
	VertexListNode *head = entry->second;
	if (head){
		VertexListNode *grandparent = head;
		VertexListNode *parent = grandparent->getNextNode();
		VertexListNode *child = parent->getNextNode();

		if (grandparent->getVertex() == b && parent->getVertex() == c){
			//tail case
			if (grandparent->isPreFaceDeleted()){
				nodePool->Dealloc(grandparent);
				//check the node behide parent
				if (child && child->isPreFaceDeleted()){
					entry->second = child;
					nodePool->Dealloc(parent);
				}
				//normal case
				else{
					parent->setDeletedMark();
					entry->second = parent;
				}
			}
			//head case
			else{
				parent->setDeletedMark();
				//check the node behide parent
				if (child && child->isPreFaceDeleted()){
					grandparent->setNextNode(child);
					nodePool->Dealloc(parent);
				}
			}
		}
		else{
			bool found = false;
			while (child != NULL){
				if (parent->getVertex() == b && child->getVertex() == c){
					found = true;
					child->setDeletedMark();
					VertexListNode *grandchild = child->getNextNode();
					//head and tail case
					if (grandchild == NULL && head->isPreFaceDeleted()){
						if (parent->isPreFaceDeleted()){
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
					else if (parent->isPreFaceDeleted()){
						if (grandchild && grandchild->isPreFaceDeleted()){
							grandparent->setNextNode(grandchild);
							nodePool->Dealloc(child);
							nodePool->Dealloc(parent);
							child = NULL;
							parent = NULL;
						}
						else{
							grandparent->setNextNode(child);
							nodePool->Dealloc(parent);
							parent = NULL;
						}
					}
					//child and grandchild are all deletd
					else if (grandchild && grandchild->isPreFaceDeleted()){
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

			if (!found){
				if (parent->getVertex() == b && head->getVertex() == c){
					VertexListNode *next = head->getNextNode();
					//the node behide the head deleted
					if (next && next->isPreFaceDeleted()){
						next->setDeletedMark();
						entry->second = next;
						nodePool->Dealloc(head);
						head = next;
					}
					//normal case
					else{
						head->setDeletedMark();
						//check the tail
						if (parent->isPreFaceDeleted()){
							grandparent->setNextNode(NULL);
							nodePool->Dealloc(parent);
						}
					}
				}
			}
		}
		//clean the deleted node
		if (head->getNextNode() == NULL){
			nodePool->Dealloc(head);
			entry->second = NULL;
		}
	}
}

REAL DelMesher::Tetrahedron::maxREration = FLT_MAX;

DelMesher::Tetrahedron::Tetrahedron(Vertex *v0, Vertex *v1, Vertex *v2, Vertex *v3, bool ordered){
	reRation = r = 0.0;
	v[0] = v0; v[1] = v1; v[2] = v2; v[3] = v3;
	if (ordered)
		sortVertices();
}

void DelMesher::Tetrahedron::setRationAndRadius(const Vertex *ghost){
	if (v[0] == ghost || v[1] == ghost || v[2] == ghost || v[3] == ghost){
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

	predicator.Orthosphere(v[0]->vert, v[0]->weight, v[1]->vert, v[1]->weight, v[2]->vert, v[2]->weight, v[3]->vert, v[3]->weight, NULL, &r);
	
	REAL minEdgeLength = sqrt(std::min(da.length2(), std::min(db.length2(), std::min(dc.length2(), std::min(ca.length2(), std::min(cb.length2(), ba.length2()))))));
	reRation = r / minEdgeLength;
}

void DelMesher::Tetrahedron::sortVertices(){
	int min = 0, max = 0;
	for (int i = 1; i < 4; i++){
		if (v[i]->getLabel() < v[min]->getLabel())
			min = i;
		else if (v[i]->getLabel() > v[max]->getLabel())
			max = i;
	}
	if (min == 0){
		if (max == 1){
			std::swap(v[1], v[2]);
			std::swap(v[2], v[3]);
		}
		else if (max == 2){
			std::swap(v[2], v[3]);
			std::swap(v[1], v[2]);
		}
	}
	else if (min == 3){
		std::swap(v[0], v[3]);
		if (max == 0){
			std::swap(v[1], v[2]);
		}
		else{
			std::swap(v[max], v[3]);
		}
	}
	else{
		std::swap(v[0], v[min]);
		if (max == 3)
			std::swap(v[1], v[2]);
		else if (max == 0)
			std::swap(v[min], v[3]);
		else
			std::swap(v[max], v[3]);
	}
}

Predicator<REAL> DelMesher::predicator;

DelMesher::DelMesher(Vector *surfvs, int *segis, int *subpolygons, int numv, int numseg, int numpol, int *numsubpol, REAL maxRa, REAL maxR){
	maxRatio = maxRa;
	maxRadius = maxR;
	Tetrahedron::maxREration = maxRa;
	nodePool = new MemoryPool<VertexListNode>(256);
	edgeNodePool = new MemoryPool<EdgeListNode>();
	//init ghost vertex
	ghost = vertArena.Alloc();
	ghost->setGhost();
	//insert every vertex
	int i = 0;
	Vertex *verts = vertArena.Alloc(numv);
	vertices.reserve(numv);
	for (i = 0; i < numv; i++){
		boundBox.Insert(surfvs[i]);
		verts[i] = Vertex(surfvs[i]);
		verts[i].setLabel();
		vertices.push_back(verts + i);
	}
	//protect every vertex
	Protect();
	oriSegments.reserve(numseg);
	for (int i = 0; i < numseg; i++){
		Segment s = Segment(vertices[segis[2 * i]], vertices[segis[2 * i + 1]], true);
		//record orgin sgments
		oriSegments.push_back(s);
		//recorded
		segments.insert(s);
		mayEncroachedSegs.push_back(s);
	}
	//triangluate every polygon
	std::vector<Vertex *> polygonVerts;
	DelTriangulator t;
	for (int i = 0; i < numpol; i++){
		for (int j = 0; j < numsubpol[i]; j++){
			polygonVerts.push_back(vertices[subpolygons[j]]);
		}
		t.generateSubPolygons(polygonVerts);
		t.outPut(this);
		subpolygons += numsubpol[i];
		polygonVerts.clear();
	}

	triangulation3D();
}

void DelMesher::triangulation3D(){
	Vertex *a = vertices[0], *b = vertices[1], *c = vertices[2];
	Vertex *d;
	for (int i = 3; i < vertices.size(); i++){
		d = vertices[i];
		if (!predicator.coPlane(a->vert, b->vert, c->vert, d->vert)){
			std::swap(vertices[3], vertices[i]);
			break;
		}
	}
	if (predicator.orient3d(a->vert, b->vert, c->vert, d->vert) < 0.0)
		std::swap(a, b);

	Tetrahedron t(a, b, c, d, true);
	t.setRationAndRadius(ghost);
	if (t.getREration() > maxRatio || t.getRadius() > maxRadius)
		skinnyTets.push(t);

	addTet(a, b, c, d);
	addTet(ghost, b, d, c);
	addTet(ghost, a, c, d);
	addTet(ghost, a, b, c);
	addTet(ghost, a, d, b);

	Tetrahedron rt = t;
	for (int i = 4; i < vertices.size(); i++){
		Vertex *v = vertices[i];
		t = findPosition(vertices[i], rt);
		insertVertex(vertices[i], t, &rt, true);
	}
	std::sort(vertices.begin(), vertices.end(), [](const Vertex *left, const Vertex *right){
		return *left < *right;
	});
}

bool DelMesher::Adjacent(const Segment &s, Vertex **w) const{
	/*auto wp = svHash.find(s);
	if (wp == svHash.end())
		return false;
	*w = wp->second;
	return true;*/

	auto entry = polygonTopology.find(s.v[0]);
	if (entry == polygonTopology.end())
		return false;
	VertexListNode *node = entry->second;
	bool found = false;
	while (node != NULL && !found){
		if (node->getVertex() == s.v[1]){
			VertexListNode *nextNode = node->getNextNode();
			if (nextNode == NULL)
				nextNode = entry->second;
			if (!nextNode->isPreFaceDeleted()){
				*w = nextNode->getVertex();
				found = true;
			}
		}
		node = node->getNextNode();
	}
	return found;
}

bool DelMesher::Adjacent(const Face &f, Vertex **z) const{
	/*bool oldFound = false;
	auto zp = fvHash.find(f);
	if (zp != fvHash.end())
		oldFound = true;*/

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

bool DelMesher::Adjacent(Vertex *w, Vertex *x, Vertex *y, Vertex **z) const{
	bool found = false;
	if (w->getLabel() > x->getLabel()){
		if (x->hasList()){
			EdgeListNode *linkHead = x->getListHead();
			while (linkHead != NULL && linkHead->getEndVertex() != w){
				linkHead = linkHead->getNextNode();
			}
			if (linkHead != NULL){
				VertexListNode *node = linkHead->getLink();
				while (node != NULL && !found){
					if (node->getVertex() == y){
						VertexListNode *foundNode = node->getNextNode();
						if (foundNode == NULL)
							foundNode = linkHead->getLink();
						if (!foundNode->isPreFaceDeleted()){
							*z = foundNode->getVertex();
							found = true;
						}
					}
					node = node->getNextNode();
				}
			}
		}
	}
	else{
		if (w->hasList()){
			EdgeListNode *linkHead = w->getListHead();
			while (linkHead != NULL && linkHead->getEndVertex() != x){
				linkHead = linkHead->getNextNode();
			}
			if (linkHead != NULL){
				VertexListNode *parentNode = NULL;
				VertexListNode *node = linkHead->getLink();
				while (node != NULL && !found){
					if (node->getVertex() == y){
						if (!node->isPreFaceDeleted()){
							VertexListNode *foundNode = parentNode;
							if (foundNode == NULL){
								VertexListNode *loop = node;
								while (loop != NULL){
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

void DelMesher::addSupplyVerts(Vertex *a, Vertex *b, Vertex *c, Vertex *d, int mode){
	int aIndex = a->getLabel();
	int bIndex = b->getLabel();
	int cIndex = c->getLabel();

	int dIndex = -1;
	if (d)
	    dIndex = d->getLabel();

	switch (mode){
	case 0:
		if (aIndex > bIndex){
			if (!a->hasList())
				a->setEndVertexPoint(b);
		}
		else{
			if (!b->hasList())
				b->setEndVertexPoint(a);
		}
		if (cIndex > dIndex){
			if (!c->hasList())
				c->setEndVertexPoint(d);
		}
		else{
			if (!d->hasList())
				d->setEndVertexPoint(c);
		}
		break;
	case 1:
		if (aIndex > bIndex){
			if (aIndex > cIndex){
				if (!a->hasList())
					a->setEndVertexPoint(c);
			}
			else{
				if (!b->hasList())
					b->setEndVertexPoint(c);
			}
		}
		else if (bIndex > cIndex){
			if (!b->hasList())
				b->setEndVertexPoint(a);
		}
		else{
			if (!c->hasList())
				c->setEndVertexPoint(b);
		}
		break;
	case 2:
	{
		int v[4];
		v[0] = aIndex; v[1] = bIndex; v[2] = cIndex; v[3] = dIndex;
		int max = 0;
		for (int i = 1; i < 4; i++){
			if (v[i] > v[max])
				max = i;
		}
		if (max == 0){
			if (!a->hasList())
				a->setEndVertexPoint(b);
		}
		else if (max == 1){
			if (!b->hasList())
				b->setEndVertexPoint(c);
		}
		else if (max == 2){
			if (!c->hasList())
				c->setEndVertexPoint(d);
		}
		else{
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

bool DelMesher::adjacent2Vertex(Vertex *w, Face *f) const{
	/*auto zp = vsHash.find(w);
	Vertex *v = zp->second;
	Vertex *x;
	if (!Adjacent(Segment(w, v), &x))
		return false;
	*f = Face(w, v, x, true);
	return true;*/

	auto entry = polygonTopology.find(w);
	if (entry == polygonTopology.end())
		return false;
	VertexListNode *node = entry->second;
	VertexListNode *nextNode = NULL;
	bool found = false;
	while (node != NULL && !found){
		nextNode = node->getNextNode();
		if (nextNode && !nextNode->isPreFaceDeleted()){
			*f = Face(w, node->getVertex(), nextNode->getVertex());
			found = true;
		}
		node = nextNode;
	}
	return found;
}

bool DelMesher::adjacent2Vertex(Vertex *w, Tetrahedron *t) const{
	/*auto sp = vfHash.find(w);
	Segment s = sp->second;
	Vertex *x;
	if (!Adjacent(Face(w, s.v[0], s.v[1], true), &x))
		return false;
	*t = Tetrahedron(x, w, s.v[0], s.v[1], true);
	return true*/

	bool found = false;
	EdgeListNode *linkHead = NULL;
	if (!w->hasList()){
		Vertex *end = w->getEndVertex();
		linkHead = end->getListHead();
		bool foundHead = false;
		while (linkHead != NULL && !foundHead){
			if (linkHead->getEndVertex() == w)
				foundHead = true;
			else
				linkHead = linkHead->getNextNode();
		}
	}
	else
		linkHead = w->getListHead();

	if (linkHead != NULL){
		VertexListNode *parentLoop = linkHead->getLink();
		VertexListNode *loop = parentLoop->getNextNode();
		while (loop && loop->isPreFaceDeleted()){
			parentLoop = loop;
			loop = loop->getNextNode();
		}
		if (loop){
			*t = Tetrahedron(w, linkHead->getEndVertex(), parentLoop->getVertex(), loop->getVertex(), true);
			found = true;
		}
	}
	return found;
}

DelMesher::Tetrahedron DelMesher::findPosition(Vertex *u, const Tetrahedron &t, bool boundary) const{
	int maxIter = 16384;
	Vertex *a = t.v[0], *b = t.v[1], *c = t.v[2], *d = t.v[3];

	while (maxIter--){
		Vertex *x = NULL;
		if (a == ghost){
			if (predicator.inOrthoHalfSpace3D(u->vert, u->weight, b->vert, b->weight, c->vert, c->weight, d->vert, d->weight))
				return t;
			Adjacent(Face(b, d, c, true), &x);
			a = x; std::swap(d, c);
		}
		else{
			bool ori0 = predicator.orient3d(u->vert, a->vert, b->vert, c->vert) > 0.0;
			bool ori1 = predicator.orient3d(u->vert, a->vert, c->vert, d->vert) > 0.0;
			bool ori2 = predicator.orient3d(u->vert, a->vert, d->vert, b->vert) > 0.0;

			int condition = ori0 + (ori1 << 1) + (ori2 << 2);

			switch (condition){
			case 0:
				if (predicator.orient3d(u->vert, d->vert, c->vert, b->vert) > 0.0){
					Adjacent(Face(d, c, b, true), &x);
					if (x == ghost)
					    return Tetrahedron(ghost, d, c, b, true);
					//Tetrahedron(x, d, c, b)
				    a = x; std::swap(b, d);
				}
				else{
					return Tetrahedron(a, b, c, d, true);
				}
				break;
			case 1:
				Adjacent(Face(a, b, c, true), &x);
				if (x == ghost)
					return Tetrahedron(ghost, a, b, c, true);
				//Tetrahedron(x, a, b, c)
				d = c; c = b; b = a; a = x;
				break;
			case 2:
				Adjacent(Face(a, c, d, true), &x);
				if (x == ghost)
					return Tetrahedron(ghost, a, c, d, true);
				//Tetrahedron(x, a, c, d)
				b = a; a = x;
				break;
			case 3:
				if (Randomnation(2) == 0){
					Adjacent(Face(a, b, c, true), &x);
					if (x == ghost)
						return Tetrahedron(ghost, a, b, c, true);
					//Tetrahedron(x, a, b, c)
					d = c; c = b; b = a; a = x;
				}
				else{
					Adjacent(Face(a, c, d, true), &x);
					if (x == ghost)
						return Tetrahedron(ghost, a, c, d, true);
					//Tetrahedron(x, a, c, d)
					b = a; a = x;
				}
				break;
			case 4:
				Adjacent(Face(a, d, b, true), &x);
				if (x == ghost)
					return Tetrahedron(ghost, a, d, b, true);
				//Tetrahedron(x, a, d, b)
				c = d; d = b; b = a; a = x;
				break;
			case 5:
				if (Randomnation(2) == 0){
					Adjacent(Face(a, b, c, true), &x);
					if (x == ghost)
						return Tetrahedron(ghost, a, b, c, true);
					//Tetrahedron(x, a, b, c)
					d = c; c = b; b = a; a = x;
				}
				else{
					Adjacent(Face(a, d, b, true), &x);
					if (x == ghost)
						return Tetrahedron(ghost, a, d, b, true);
					//Tetrahedron(x, a, d, b)
					c = d; d = b; b = a; a = x;
				}
				break;
			case 6:
				if (Randomnation(2) == 0){
					Adjacent(Face(a, c, d, true), &x);
					if (x == ghost)
						return Tetrahedron(ghost, a, c, d, true);
					//Tetrahedron(x, a, c, d)
					b = a; a = x;
				}
				else{
					Adjacent(Face(a, d, b, true), &x);
					if (x == ghost)
						return Tetrahedron(ghost, a, d, b, true);
					//Tetrahedron(x, a, d, b)
					c = d; d = b; b = a; a = x;
				}
				break;
			case 7:
			{
				unsigned int term = Randomnation(3);
				if (term == 0){
					Adjacent(Face(a, b, c, true), &x);
					if (x == ghost)
						return Tetrahedron(ghost, a, b, c, true);
					//Tetrahedron(x, a, b, c)
					d = c; c = b; b = a; a = x;
				}
				else if (term == 1){
					Adjacent(Face(a, c, d, true), &x);
					if (x == ghost)
						return Tetrahedron(ghost, a, c, d, true);
					b = a; a = x;
				}
				else{
					Adjacent(Face(a, d, b, true), &x);
					if (x == ghost)
						return Tetrahedron(ghost, a, d, b, true);
					//Tetrahedron(x, a, d, b)
					c = d; d = b; b = a; a = x;
				}
				break;
			}
			default:
				Severe("Unexpected Case in 3d overloaded DelMesher::findPosition");
			}
		}
	}
	Severe("Do not find Tet in 3d overloaded DelMesher::findPosition");
	//should never return this value
	return Tetrahedron(NULL, NULL, NULL, NULL);
}

Face DelMesher::findPosition(Vertex *u, const Face& f) const{
	size_t maxIter = std::max(vertices.size(), size_t(16384));
	Vertex *a = f.v[0];
	Vertex *b = f.v[1], *c = f.v[2];

	while (maxIter--){
		Vertex *x = NULL;
		DelVector n = Normalize((a->vert - c->vert) % (b->vert - c->vert));
		bool ori0 = predicator.orient2d(b->vert, a->vert, u->vert, n) > 0.0;
		bool ori1 = predicator.orient2d(a->vert, c->vert, u->vert, n) > 0.0;

		int condition = ori0 + (ori1 << 1);
		switch (condition){
		case 0:
			if (predicator.orient2d(c->vert, b->vert, u->vert, n) > 0.0){
				Adjacent(Segment(c, b), &x);
				a = x; std::swap(b, c); //Face(x, c, b)
			}
			else{
				return Face(a, b, c, true);
			}
			break;
		case 1:
			Adjacent(Segment(b, a), &x);
			c = a; a = x; //Face(x, b, a)
			break;
		case 2:
			Adjacent(Segment(a, c), &x);
			b = a; a = x; //Face(x, a, c)
			break;
		case 3:
			if (Randomnation(2) == 0){
				Adjacent(Segment(b, a), &x);
				c = a; a = x; //Face(x, b, a)
			}
			else{
				Adjacent(Segment(a, c), &x);
				b = a; a = x; //Face(x, a, c)
			}
			break;
		default:
			Severe("Unexpected Case in 2d overloaded DelMesher::findPosition");
		}
	}
	Severe("Do not find Face in 2d overloaded DelMesher::findPosition");
	//should never return this value
	return Face(NULL, NULL, NULL);
}

void DelMesher::addTriangle(Vertex *a, Vertex *b, Vertex *c){
	/*Face f(a, b, c, true);
	polygons.insert(f);
	svHash[Segment(a, b)] = c;
	svHash[Segment(b, c)] = a;
	svHash[Segment(c, a)] = b;

	vsHash[a] = b;
	vsHash[b] = c;
	vsHash[c] = a;*/

	insertToFaceTopology(a, b, c);
	insertToFaceTopology(b, c, a);
	insertToFaceTopology(c, a, b);
}

void DelMesher::deleteTriangle(Vertex *a, Vertex *b, Vertex *c){
	/*Face f(a, b, c, true);
	polygons.erase(f);
	svHash.erase(Segment(a, b));
	svHash.erase(Segment(b, c));
	svHash.erase(Segment(c, a));*/

	removeFromFaceTopology(a, b, c);
	removeFromFaceTopology(b, c, a);
	removeFromFaceTopology(c, a, b);
}

void DelMesher::addTet(Vertex *a, Vertex *b, Vertex *c, Vertex *d){
	/*Tetrahedron t(a, b, c, d, true);
	if (t.v[0] != ghost)
		tets.insert(t);

	fvHash[Face(t.v[1], t.v[2], t.v[3], true)] = t.v[0];
	fvHash[Face(t.v[0], t.v[3], t.v[2], true)] = t.v[1];
	fvHash[Face(t.v[0], t.v[1], t.v[3], true)] = t.v[2];
	fvHash[Face(t.v[0], t.v[2], t.v[1], true)] = t.v[3];

	vfHash[t.v[0]] = Segment(t.v[2], t.v[1]);
	vfHash[t.v[1]] = Segment(t.v[2], t.v[3]);
	vfHash[t.v[2]] = Segment(t.v[0], t.v[3]);
	vfHash[t.v[3]] = Segment(t.v[0], t.v[1]);*/

	bool ab = parityCheck(a, b);
	bool ac = parityCheck(a, c);
	bool ad = parityCheck(a, d);

	if (ab){
		if (ac){
			if (ad){
				insertToTetTopology(Segment(a, b), c, d);
				insertToTetTopology(Segment(a, c), d, b);
				insertToTetTopology(Segment(a, d), b, c);
				insertToTetTopology(Segment(b, c), a, d);
				insertToTetTopology(Segment(b, d), c, a);
				insertToTetTopology(Segment(c, d), a, b);
				addSupplyVerts(a, b, c, d, 2);
			}
			else{
				insertToTetTopology(Segment(a, b), c, d);
				insertToTetTopology(Segment(a, c), d, b);
				insertToTetTopology(Segment(b, c), a, d);

				addSupplyVerts(a, b, c, NULL, 1);
			}
		}
		else if (ad){
			insertToTetTopology(Segment(a, b), c, d);
			insertToTetTopology(Segment(a, d), b, c);
			insertToTetTopology(Segment(b, d), c, a);

			addSupplyVerts(a, b, d, NULL, 1);
		}
		else{
			insertToTetTopology(Segment(a, b), c, d);
			insertToTetTopology(Segment(c, d), a, b);

			addSupplyVerts(a, b, c, d, 0);
		}
	}
	else if (ac){
		if (ad){
			insertToTetTopology(Segment(a, c), d, b);
			insertToTetTopology(Segment(a, d), b, c);
			insertToTetTopology(Segment(c, d), a, b);

			addSupplyVerts(a, c, d, NULL, 1);
		}
		else{
			insertToTetTopology(Segment(a, c), d, b);
			insertToTetTopology(Segment(b, d), c, a);

			addSupplyVerts(a, c, b, d, 0);
		}
	}
	else if (ad){
		insertToTetTopology(Segment(a, d), b, c);
		insertToTetTopology(Segment(b, c), a, d);

		addSupplyVerts(a, d, b, c, 0);
	}
	else{
		insertToTetTopology(Segment(b, c), a, d);
		insertToTetTopology(Segment(b, d), c, a);
		insertToTetTopology(Segment(c, d), a, b);

		addSupplyVerts(b, c, d, NULL, 1);
	}

}

void DelMesher::deleteTet(Vertex *a, Vertex *b, Vertex *c, Vertex *d){
	/*Tetrahedron t(a, b, c, d, true);
	if (t.v[0] != ghost)
		tets.erase(t);

	fvHash.erase(Face(t.v[3], t.v[1], t.v[2], true));
	fvHash.erase(Face(t.v[0], t.v[3], t.v[2], true));
	fvHash.erase(Face(t.v[0], t.v[1], t.v[3], true));
	fvHash.erase(Face(t.v[0], t.v[2], t.v[1], true));*/

	bool ab = parityCheck(a, b);
	bool ac = parityCheck(a, c);
	bool ad = parityCheck(a, d);

	if (ab){
		if (ac){
			if (ad){
				removeFromTetTopology(Segment(a, b), c, d);
				removeFromTetTopology(Segment(a, c), d, b);
				removeFromTetTopology(Segment(a, d), b, c);
				removeFromTetTopology(Segment(b, c), a, d);
				removeFromTetTopology(Segment(b, d), c, a);
				removeFromTetTopology(Segment(c, d), a, b);
			}
			else{
				removeFromTetTopology(Segment(a, b), c, d);
				removeFromTetTopology(Segment(a, c), d, b);
				removeFromTetTopology(Segment(b, c), a, d);
			}
		}
		else if (ad){
			removeFromTetTopology(Segment(a, b), c, d);
			removeFromTetTopology(Segment(a, d), b, c);
			removeFromTetTopology(Segment(b, d), c, a);
		}
		else{
			removeFromTetTopology(Segment(a, b), c, d);
			removeFromTetTopology(Segment(c, d), a, b);
		}
	}
	else if (ac){
		if (ad){
			removeFromTetTopology(Segment(a, c), d, b);
			removeFromTetTopology(Segment(a, d), b, c);
			removeFromTetTopology(Segment(c, d), a, b);
		}
		else{
			removeFromTetTopology(Segment(a, c), d, b);
			removeFromTetTopology(Segment(b, d), c, a);
		}
	}
	else if (ad){
		removeFromTetTopology(Segment(a, d), b, c);
		removeFromTetTopology(Segment(b, c), a, d);
	}
	else{
		removeFromTetTopology(Segment(b, c), a, d);
		removeFromTetTopology(Segment(b, d), c, a);
		removeFromTetTopology(Segment(c, d), a, b);
	}

}

void DelMesher::digCavity(Vertex *u, const Face &f, Tetrahedron *rt, bool insertToSkinny, bool trulyDeleteOrAdd, bool boundaryVert){
	Vertex *a;
	if (!Adjacent(f, &a))
		return;

	Vertex *b = f.v[0];
	Vertex *c = f.v[1];
	Vertex *d = f.v[2];

	bool deleted = false;
	if (a == ghost)
	    deleted = predicator.inOrthoHalfSpace3D(u->vert, u->weight, b->vert, b->weight, c->vert, c->weight, d->vert, d->weight, boundaryVert);
	else if (b == ghost){
		deleted = predicator.inOrthoHalfSpace3D(u->vert, u->weight, a->vert, a->weight, d->vert, d->weight, c->vert, c->weight, boundaryVert);
	}
	else
		deleted = predicator.inOrthoSpherePerturbed(u->vert, u->weight, b->vert, b->weight, d->vert, d->weight, c->vert, c->weight, a->vert, a->weight) > 0.0;

	if (deleted){
		if (!trulyDeleteOrAdd)
			tobeDeletedTets.push_back(Tetrahedron(a, b, c, d));
		deleteTet(a, b, c, d);
		digCavity(u, Face(a, c, d, true), rt, insertToSkinny, trulyDeleteOrAdd, boundaryVert);
		digCavity(u, Face(a, b, c, true), rt, insertToSkinny, trulyDeleteOrAdd, boundaryVert);
		digCavity(u, Face(a, d, b, true), rt, insertToSkinny, trulyDeleteOrAdd, boundaryVert);
	}
	else{
		if (trulyDeleteOrAdd){
			addTet(u, b, d, c);
		}
		else
			newFacesOfTets.push_back(Face(b, d, c, true));
		if (rt)
			*rt = Tetrahedron(u, b, d, c, true);
		if (insertToSkinny){
			Tetrahedron t(u, b, d, c, true);
			t.setRationAndRadius(ghost);
			if (t.getREration() > maxRatio || t.getRadius() > maxRadius)
				skinnyTets.push(t);

				Face faceOut(b, c, d, true), faceIn(b, d, c, true);
				if (findSubPolygons(faceOut))
					mayEncroachedFaces.push_back(faceOut);
				if (findSubPolygons(faceIn))
					mayEncroachedFaces.push_back(faceIn);
		}
	}
}

void DelMesher::digCavity(Vertex *u, const Vertex& aboveVert,const Segment &s, bool insertToQueue, bool trulyDeleteOrAdd){
	Vertex *x;
	if (!Adjacent(s, &x))
		return;

	Vertex *w = s.v[0];
	Vertex *v = s.v[1];
	bool deleted = false;
	if (!findSegment(Segment(w, v, true)))
		deleted = predicator.inOrthoSpherePerturbed(aboveVert.vert, aboveVert.weight, u->vert, u->weight, v->vert, v->weight, w->vert, w->weight, x->vert, x->weight) > 0.0;

	if (deleted){
		if (!trulyDeleteOrAdd)
			tobeDeletedFaces.push_back(Face(w, v, x, true));
		deleteTriangle(w, v, x);
		digCavity(u, aboveVert, Segment(x, v), insertToQueue, trulyDeleteOrAdd);
		digCavity(u, aboveVert, Segment(w, x), insertToQueue, trulyDeleteOrAdd);
	}
	else{
		if (trulyDeleteOrAdd)
			addTriangle(u, v, w);
		else
			newSegsOfFaces.push_back(Segment(v, w));
		if (insertToQueue){
			mayEncroachedFaces.push_back(Face(u, v, w, true));
			Segment s(v, w, true);
			if (findSegment(s))
				mayEncroachedSegs.push_back(s);
		}
	}
}

Face DelMesher::findFaceAroundOnPlane(const Vertex& origin, Vertex *center) const{
	int indice = -1;
	Face startFace;
	Face ret;
	Vertex *b, *c, *mark;

	if (!adjacent2Vertex(center, &startFace)){
		adjacent2Vertex(center, &startFace);
	}
	if (startFace.v[0] == center)
		indice = 0;
	else if (startFace.v[1] == center)
		indice = 1;
	else
		indice = 2;
	b = startFace.v[NEXT_F(indice)];
	c = startFace.v[NEXT_F(indice + 1)];
	mark = b;

	indice = 0;
	do{
		if (predicator.coPlane(center->vert, b->vert, c->vert, origin.vert)){
			ret = Face(center, b, c);
			break;
		}
		Adjacent(Segment(center, c), &b);
		std::swap(b, c);
	} while (b != mark);
	return ret;
}

void DelMesher::splitSubSegment(const Segment &s){
	Vertex *a = s.v[0], *b = s.v[1];
	Vertex *mid = Cover(s);

	Vertex *c = NULL, *d = NULL;
	//should be orient inward
	Face f0, f1;

	if (!Adjacent(s, &c) || !Adjacent(Segment(b, a), &d)){
		Severe("Unexpected error in DelMesher::splitSubSegment");
	}
	else{
		Vertex aboveVert0(a->vert + Normalize((a->vert - c->vert) % (b->vert - c->vert)));
		Vertex aboveVert1(b->vert + Normalize((b->vert - d->vert) % (a->vert - d->vert)));
		deleteTriangle(a, b, c);
		digCavity(mid, aboveVert0, Segment(c, b));
		digCavity(mid, aboveVert0, Segment(a, c));

		deleteTriangle(b, a, d);
		digCavity(mid, aboveVert1, Segment(d, a));
		digCavity(mid, aboveVert1, Segment(b, d));

		f0 = Face(a, b, c, true);
		f1 = Face(b, a, d, true);
	}

	Tetrahedron toBeDeleted;
	Vertex *e = NULL;
	if (Adjacent(f0, &e))
		toBeDeleted = Tetrahedron(e, f0.v[0], f0.v[1], f0.v[2], true);
	else if (Adjacent(f1, &e))
		toBeDeleted = Tetrahedron(e, f1.v[0], f1.v[1], f1.v[2], true);
	else{
		adjacent2Vertex(s.v[0], &toBeDeleted);
		toBeDeleted = findPosition(mid, toBeDeleted);
	}
	insertVertex(mid, toBeDeleted, NULL, true, true);
}

void DelMesher::splitSubPolygon(const Face &f){
	Vertex *a = f.v[0], *b = f.v[1], *c = f.v[2];
	Vertex center; 
	predicator.Orthocircle(a->vert, a->weight, b->vert, b->weight, c->vert, c->weight, &center.vert);
	Vertex aboveVert(a->vert + Normalize((a->vert - c->vert) % (b->vert - c->vert)));

	tobeDeletedFaces.push_back(f);
	deleteTriangle(a, b, c);

	digCavity(&center, aboveVert, Segment(b, a), false, false);
	digCavity(&center, aboveVert, Segment(a, c), false, false);
	digCavity(&center, aboveVert, Segment(c, b), false, false);

	bool encorach = false;
	Segment e;
	for (auto ns : newSegsOfFaces){
		Segment seg = Segment(ns.v[0], ns.v[1], true);
		if (findSegment(seg)){
			if (Encroached(seg, center)){
				encorach = true;
				e = seg;
				break;
			}
		}
	}

	if (!encorach){
		Vertex *newPoint = allocVertex(center);
		vertices.push_back(newPoint);
		for (auto segOfFace : newSegsOfFaces){
			addTriangle(newPoint, segOfFace.v[0], segOfFace.v[1]);
			mayEncroachedFaces.push_back(Face(newPoint, segOfFace.v[0], segOfFace.v[1], true));
		}
		Tetrahedron toBeDeleted;
		Vertex *d;
		if (Adjacent(f, &d))
			toBeDeleted = Tetrahedron(d, a, b, c, true);
		else{
			adjacent2Vertex(a, &toBeDeleted);
			toBeDeleted = findPosition(newPoint, toBeDeleted);
		}
		if (toBeDeleted.v[0] == NULL){
			bool test5 = predicator.inOrthoHalfSpace3D(center.vert, center.weight, toBeDeleted.v[1]->vert, toBeDeleted.v[1]->weight, toBeDeleted.v[2]->vert, toBeDeleted.v[2]->weight, toBeDeleted.v[3]->vert, toBeDeleted.v[3]->weight);
			float x = 1.f;
		}
		insertVertex(newPoint, toBeDeleted, NULL, true, true);
	}
	else{
		for (auto i : tobeDeletedFaces){
			addTriangle(i.v[0], i.v[1], i.v[2]);
		}
		splitSubSegment(e);
	}
	tobeDeletedFaces.clear();
	newSegsOfFaces.clear();
}

void DelMesher::splitTetrahedron(const Tetrahedron& tet){
	Vertex *a = tet.v[0], *b = tet.v[1], *c = tet.v[2], *d = tet.v[3];
	Vertex center;
	predicator.Orthosphere(a->vert, a->weight, b->vert, b->weight, 
		c->vert, c->weight, d->vert, d->weight, &center.vert);

	tobeDeletedTets.push_back(tet);
	deleteTet(a, b, c, d);

	digCavity(&center, Face(b, d, c, true), NULL, false, false);
	digCavity(&center, Face(a, c, d, true), NULL, false, false);
	digCavity(&center, Face(a, d, b, true), NULL, false, false);
	digCavity(&center, Face(a, b, c, true), NULL, false, false);

	bool encrochment = false;
	Face encrocachedFace;
	bool result = false;
	for (auto nf : newFacesOfTets){
		result = findSubPolygons(nf);
		if(!result){
			std::swap(nf.v[2], nf.v[1]);
			result = findSubPolygons(nf);
			if (!result)
				continue;
		}
		if (Encroached(nf, center)){
			encrocachedFace = nf;
			encrochment = true;
			break;
		}
	}

	if (!encrochment){
		Vertex *newVert = allocVertex(center);
		vertices.push_back(newVert);
		for (auto f : newFacesOfTets){
			addTet(newVert, f.v[0], f.v[1], f.v[2]);
			Tetrahedron t(newVert, f.v[0], f.v[1], f.v[2], true);
			t.setRationAndRadius(ghost);
			if (t.getREration() > maxRatio || t.getRadius() > maxRadius)
			    skinnyTets.push(t);
		}
	}
	else{
		for (auto t : tobeDeletedTets)
			addTet(t.v[0], t.v[1], t.v[2], t.v[3]);
		splitSubPolygon(encrocachedFace);
	}
	tobeDeletedTets.clear();
	newFacesOfTets.clear();
}

Vertex* DelMesher::allocVertex(const DelVector &vert, REAL weight){
	Vertex *newVertex = vertArena.Alloc();
	newVertex->vert = vert;
	newVertex->weight = weight;
	newVertex->setLabel();
	Assert(boundBox.Inside(newVertex->vert));
	return newVertex;
}

Vertex* DelMesher::allocVertex(const Vertex &vertex){
	Vertex *newVertex = vertArena.Alloc();
	*newVertex = vertex;
	newVertex->setLabel();
	Assert(boundBox.Inside(newVertex->vert));
	return newVertex;
}

void DelMesher::insertVertex(Vertex *u, const Tetrahedron& tet, Tetrahedron *rt, bool insertToSkinny, bool boundaryVert){
	Vertex *a = tet.v[0], *b = tet.v[1], *c = tet.v[2], *d = tet.v[3];
	deleteTet(a, b, c, d);
	digCavity(u, Face(b, d, c, true), rt, insertToSkinny, true, boundaryVert);
	digCavity(u, Face(a, c, d, true), NULL, insertToSkinny, true, boundaryVert);
	digCavity(u, Face(a, d, b, true), NULL, insertToSkinny, true, boundaryVert);
	digCavity(u, Face(a, b, c, true), NULL, insertToSkinny, true, boundaryVert);
}

void DelMesher::insertSurfaceVertex(Vertex *u, const Face& f, bool insertToQueue){
	Vertex *a = f.v[0], *b = f.v[1], *c = f.v[2];
	Vertex aboveVert(a->vert + Normalize((a->vert - c->vert) % (b->vert - c->vert)));
	deleteTriangle(a, b, c);
	digCavity(u, aboveVert, Segment(b, a), insertToQueue);
	digCavity(u, aboveVert, Segment(c, b), insertToQueue);
	digCavity(u, aboveVert, Segment(a, c), insertToQueue);
}


bool DelMesher::Encroached(const Segment& s) const{
	Segment s2(s.v[1], s.v[0]);
	Vertex *x, *y;
	if (!Adjacent(s, &x) || !Adjacent(s2, &y))
		return true;

	DelVector center;
	REAL r = 0.0;
	predicator.OrthoCenter(s.v[0]->vert, s.v[0]->weight, s.v[1]->vert, s.v[1]->weight, &center, &r);
	REAL r2 = r*r;

	if ((x->vert - center).length2() - x->weight <= r2)
		return true;
	if ((y->vert - center).length2() - y->weight <= r2)
		return true;

	return false;
}

bool DelMesher::Encroached(const Face& f) const{
	Face f2(f.v[0], f.v[2], f.v[1], true);
	Vertex *x, *y;
	if (!Adjacent(f, &x) || !Adjacent(f2, &y))
		return true;
	Vertex c;
	REAL r = 0.0;
	predicator.Orthocircle(f.v[0]->vert, f.v[0]->weight, f.v[1]->vert, f.v[1]->weight, f.v[2]->vert, f.v[2]->weight, &c.vert, &r);
	REAL r2 = r*r;

	bool encroached0 = false;
	bool encroached1 = false;
	if (x != ghost)
		encroached0 = ((x->vert - c.vert).length2() - x->weight) <= r2;
	if (y != ghost)
		encroached1 = ((y->vert - c.vert).length2() - y->weight) <= r2;

	if (encroached0 || encroached1)
		return true;
	return false;
}

bool DelMesher::Encroached(const Segment &s, const Vertex &v) const{
	DelVector center;
	REAL r = 0.0;
	predicator.OrthoCenter(s.v[0]->vert, s.v[0]->weight, s.v[1]->vert, s.v[1]->weight, &center, &r);

	if ((v.vert - center).length2() - v.weight <= r*r)
		return true;

	return false;
}

bool DelMesher::Encroached(const Face &f, const Vertex &v) const{
	Vertex center;
	REAL r;
	predicator.Orthocircle(f.v[0]->vert, f.v[0]->weight, f.v[1]->vert, f.v[1]->weight, f.v[2]->vert, f.v[2]->weight, &center.vert, &r);
	REAL r2 = r*r;

	return ((v.vert - center.vert).length2() - v.weight) <= r2;
}

void DelMesher::Protect(){
	int size = vertices.size();
	REAL maxRadius2 = maxRadius*maxRadius;
	REAL *lfses = new REAL[size];
	for (int i = 0; i < size; i++){
		lfses[i] = FLT_MAX;
	}
	for (int i = 0; i < size; i++){
		Vertex *v = vertices[i];
		for (int j = i + 1; j < size; j++){
			REAL lfs = (v->vert - vertices[j]->vert).length2();
			if (lfs < lfses[i]) lfses[i] = lfs;
			if (lfs < lfses[j]) lfses[j] = lfs;
		}
		v->weight = std::min(maxRadius2, lfses[i] * REAL(0.125));
	}
	delete[] lfses;
}

//segment s should not in DelMesher::segments before calling Cover
Vertex* DelMesher::Cover(const Segment &s){
	DelVector a = s.v[0]->vert;
	DelVector b = s.v[1]->vert;
	REAL aWeight = s.v[0]->weight, bWeight = s.v[1]->weight;
	REAL maxRadius2 = maxRadius*maxRadius;
	//compute the center on s of the ball orthogonal to both end
	REAL weight = 0.0;
	DelVector v;
	predicator.OrthoCenter(a, aWeight, b, bWeight, &v, &weight);
	weight *= weight;

	REAL localGapSize2 = estimateLocalGapSize2(v);

	Segment s0, s1;
	Vertex *newVert = NULL;
	//determine whether the protection ends here
	if (weight <= std::min(maxRadius2, localGapSize2*REAL(0.25))){
		newVert = allocVertex(v, weight);
		vertices.push_back(newVert);
		s0 = Segment(newVert, s.v[0], true);
		s1 = Segment(newVert, s.v[1], true);
	}
	//add weighted point to protect
	else{
		DelVector ab = b - a;
		REAL abLen = ab.length();
		REAL ra = sqrt(aWeight);
		REAL rb = sqrt(bWeight);
		REAL extraLength = (ra - rb) / (REAL(2.0) * abLen);
		REAL gapLength = abLen - ra - rb;
		newVert = allocVertex(a + (REAL(0.5) + extraLength)*ab, std::min(maxRadius2, std::min(localGapSize2 * REAL(0.25), gapLength*gapLength / (REAL(36.0) - REAL(16.0)*SQRTF_2))));
		vertices.push_back(newVert);
		s0 = Segment(newVert, s.v[0], true);
		s1 = Segment(newVert, s.v[1], true);
	}

	segments.erase(s);
	segments.insert(s0);
	segments.insert(s1);
	mayEncroachedSegs.push_back(s0);
	mayEncroachedSegs.push_back(s1);

	return newVert;
}

REAL DelMesher::estimateLocalGapSize2(const DelVector &c) const{
	REAL localGapSize2 = FLT_MAX;
	int oriSegmentCounts = oriSegments.size();
	for (int i = 0; i < oriSegments.size();i++){
		Segment segment = oriSegments[i];
		DelVector a = segment.v[0]->vert, b = segment.v[1]->vert;
		REAL size = FLT_MAX;
		if (predicator.coLine(a, b, c)){
			size = std::max((c - a).length2(), (c - b).length2());
		}
		else{
			DelVector ab = b - a;
			REAL t = (c - a)*ab / ab.length2();
			if (t <= 0.0) size = (c - a).length2();
			else if (t >= 1.f) size = (c - b).length2();
			else{
				DelVector d = a + t*ab;
				size = (c - d).length2();
			}
		}
		if (size < localGapSize2){
			localGapSize2 = size;
		}
	}
	return localGapSize2;
}


void DelMesher::insertToFaceTopology(Vertex *a, Vertex *b, Vertex *c){
	auto entry = polygonTopology.find(a);
	if (entry == polygonTopology.end() || entry->second == NULL){
		VertexListNode *node = nodePool->Alloc();
		node->setVertex(b);
		node->setDeletedMark();
		VertexListNode *nextNode = nodePool->Alloc();
		nextNode->setVertex(c);
		node->setNextNode(nextNode);
		polygonTopology[a] = node;
	}
	else{
		VertexListNode *head = entry->second;
		VertexListNode *parent = head;
		VertexListNode *child = parent->getNextNode();
		VertexListNode *foundNode[2];
		VertexListNode *parentFoundNode;
		int found = 0;
		if (head->getVertex() == b || head->getVertex() == c){
			foundNode[found++] = head;
		}
		while (child != NULL && found < 2){
			Vertex *foundVert = child->getVertex();
			if (foundVert == b)
				foundNode[found++] = child;
			else if (foundVert == c){
				parentFoundNode = parent;
				foundNode[found++] = child;
			}
			parent = child;
			child = child->getNextNode();
		}

		switch (found){
		case 0:
		{
			VertexListNode *newNode = nodePool->Alloc();
			newNode->setVertex(b);
			newNode->setDeletedMark();
			VertexListNode *nextNode = nodePool->Alloc();
			nextNode->setVertex(c);
			newNode->setNextNode(nextNode);

			if (head->isPreFaceDeleted()){
				nextNode->setNextNode(head);
				head = newNode;
				entry->second = head;
			}
			else{
				parent = head;
				child = parent->getNextNode();
				while (child && !child->isPreFaceDeleted()){
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
			if (foundNode[0]->isPreFaceDeleted()){
				newNode->setVertex(b);
				if (foundNode[0] == head){
					head->unSetDeletedMark();
					newNode->setDeletedMark();
					newNode->setNextNode(head);
					head = newNode;
					entry->second = head;
				}
				else{
					foundNode[0]->unSetDeletedMark();
					newNode->setNextNode(foundNode[0]);
					newNode->setDeletedMark();
					parentFoundNode->setNextNode(newNode);
				}
			}
			else{
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
			if (unoredered){
				bNode = foundNode[1];  cNode = foundNode[0];
			}
			else{
				bNode = foundNode[0];  cNode = foundNode[1];
			}
			if (bNode->getNextNode() == cNode || (cNode == head && bNode->getNextNode() == NULL))
				cNode->unSetDeletedMark();
			else{
				//find the lost fan
				VertexListNode *start = bNode->getNextNode();
				VertexListNode *end = parentFoundNode;

				if (unoredered){
					if (head == cNode){
						VertexListNode *loop = start;
						while (loop != NULL){
							end = loop;
							loop = loop->getNextNode();
						}
						bNode->setNextNode(NULL);
					}
					else{
						//the lost fan incule tail->head
						VertexListNode *loop = start;
						if (loop == NULL){
							start = head;
							loop = start;
						}
						else{
							while (loop != NULL){
								end = loop;
								loop = loop->getNextNode();
							}
							end->setNextNode(head);
							loop = head;
						}

						while (loop != cNode){
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
				while (gapEnd && !gapEnd->isPreFaceDeleted()){
					gapStart = gapEnd;
					gapEnd = gapEnd->getNextNode();
				}
				//insert it
				if (gapEnd != NULL){
					if (gapStart != NULL){
						gapStart->setNextNode(start);
						end->setNextNode(gapEnd);
					}
					else{
						end->setNextNode(head);
						head = start;
						entry->second = head;
					}
				}
				//the lost fan lost
				else{
					VertexListNode *loop = start;
					while (loop && loop != end){
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

void DelMesher::removeFromFaceTopology(Vertex *a, Vertex *b, Vertex *c){
	auto entry = polygonTopology.find(a);
	Assert(entry != polygonTopology.end());
	VertexListNode *head = entry->second;
	if (head){
		VertexListNode *grandparent = head;
		VertexListNode *parent = grandparent->getNextNode();
		VertexListNode *child = parent->getNextNode();

		if (grandparent->getVertex() == b && parent->getVertex() == c){
			//tail case
			if (grandparent->isPreFaceDeleted()){
				nodePool->Dealloc(grandparent);
				//check the node behide parent
				if (child && child->isPreFaceDeleted()){
					head = child;
					entry->second = head;
					nodePool->Dealloc(parent);
				}
				else{
					parent->setDeletedMark();
					head = parent;
					entry->second = head;
				}
			}
			//head case
			else{
				parent->setDeletedMark();
				//check the node behide parent
				if (child && child->isPreFaceDeleted()){
					grandparent->setNextNode(child);
					nodePool->Dealloc(parent);
				}
			}
		}
		else{
			bool found = false;
			while (child != NULL){
				if (parent->getVertex() == b && child->getVertex() == c){
					found = true;
					child->setDeletedMark();
					VertexListNode *grandchild = child->getNextNode();
					//head and tail case
					if (grandchild == NULL && head->isPreFaceDeleted()){
						if (parent->isPreFaceDeleted()){
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
					else if (parent->isPreFaceDeleted()){
						if (grandchild && grandchild->isPreFaceDeleted()){
							grandparent->setNextNode(grandchild);
							nodePool->Dealloc(child);
							nodePool->Dealloc(parent);
							child = NULL;
							parent = NULL;
						}
						else{
							grandparent->setNextNode(child);
							nodePool->Dealloc(parent);
							parent = NULL;
						}
					}
					//child and grandchild are all deletd
					else if (grandchild && grandchild->isPreFaceDeleted()){
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

			if (!found){
				if (parent->getVertex() == b && head->getVertex() == c){
					//check the tail
					if (parent->isPreFaceDeleted()){
						grandparent->setNextNode(NULL);
						nodePool->Dealloc(parent);
					}

					VertexListNode *next = head->getNextNode();
					//the node behide the head deleted
					if (next && next->isPreFaceDeleted()){
						next->setDeletedMark();
						entry->second = next;
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
		if (head->getNextNode() == NULL){
			nodePool->Dealloc(head);
			entry->second = NULL;
		}
	}
}

void DelMesher::insertToTetTopology(const Segment& s, Vertex *mayC, Vertex *mayD){
	Vertex *a = s.v[0], *b = s.v[1], *c = mayC, *d = mayD;
	if (*a > *b){
		std::swap(a, b);
		std::swap(c, d);
	}

	if (!a->hasList()){
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

	else{
		EdgeListNode *linkHead = a->getListHead();
		bool foundHead = false;
		while (linkHead != NULL && !foundHead){
			if (linkHead->getEndVertex() == b)
				foundHead = true;
			else
			    linkHead = linkHead->getNextNode();
		}
		if (!foundHead){
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
		else{
			VertexListNode *head = linkHead->getLink();
			VertexListNode *parent = head;
			VertexListNode *child = parent->getNextNode();
			VertexListNode *foundNode[2];
			VertexListNode *parentFoundNode = NULL;
			int found = 0;
			if (head->getVertex() == c || head->getVertex() == d){
				foundNode[found++] = head;
			}
			while (child != NULL && found < 2){
				Vertex *foundVert = child->getVertex();
				if (foundVert == c)
					foundNode[found++] = child;
				else if (foundVert == d){
					parentFoundNode = parent;
					foundNode[found++] = child;
				}
				parent = child;
				child = child->getNextNode();
			}

			switch (found){
			case 0:
			{
				VertexListNode *newNode = nodePool->Alloc();
				newNode->setVertex(c);
				newNode->setDeletedMark();
				VertexListNode *nextNode = nodePool->Alloc();
				nextNode->setVertex(d);
				newNode->setNextNode(nextNode);

				if (head->isPreFaceDeleted()){
					nextNode->setNextNode(head);
					head = newNode;
					linkHead->setLink(head);
				}
				else{
					parent = head;
					child = parent->getNextNode();
				    while (child && !child->isPreFaceDeleted()){
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
				if (foundNode[0]->isPreFaceDeleted()){
					Assert(foundNode[0]->getVertex() == d);
					newNode->setVertex(c);
					if (foundNode[0] == head){
						head->unSetDeletedMark();
						newNode->setDeletedMark();
						newNode->setNextNode(head);
						head = newNode;
						linkHead->setLink(head);
					}
					else{
						foundNode[0]->unSetDeletedMark();
						newNode->setNextNode(foundNode[0]);
						newNode->setDeletedMark();
						parentFoundNode->setNextNode(newNode);
					}
				}
				else{
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
				if (unoredered){
					cNode = foundNode[1];  dNode = foundNode[0];
				}
				else{
					cNode = foundNode[0];  dNode = foundNode[1];
				}
				if (cNode->getNextNode() == dNode || (dNode == head && cNode->getNextNode() == NULL))
					dNode->unSetDeletedMark();
				else{
					//find the lost fan
					VertexListNode *start = cNode->getNextNode();
					VertexListNode *end = parentFoundNode;

					if (unoredered){
                        if (head == dNode){
							VertexListNode *loop = start;
							while (loop != NULL){
								end = loop;
								loop = loop->getNextNode();
							}
							cNode->setNextNode(NULL);
						}
						else{
							//the lost fan incule tail->head
							VertexListNode *loop = start;
							if (loop == NULL){
								start = head;
								loop = start;
							}
							else{
								while (loop != NULL){
									end = loop;
									loop = loop->getNextNode();
								}
								end->setNextNode(head);
								loop = head;
							}

							while (loop != dNode){
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
					while (gapEnd && !gapEnd->isPreFaceDeleted()){
						gapStart = gapEnd;
						gapEnd = gapEnd->getNextNode();
					}
					//insert it
					if (gapEnd != NULL){
						if (gapStart != NULL){
							gapStart->setNextNode(start);
							end->setNextNode(gapEnd);
						}
						else{
							end->setNextNode(head);
							head = start;
							linkHead->setLink(head);
						}
					}
					//the lost fan lost forever
					else{
						VertexListNode *loop = start;
						while (loop && loop != end){
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

void DelMesher::removeFromTetTopology(const Segment &s, Vertex *mayC, Vertex *mayD){
	Vertex *a = s.v[0], *b = s.v[1], *c = mayC, *d = mayD;
	if (*a > *b){
		std::swap(a, b);
		std::swap(c, d);
	}

	Assert(a->hasList());
	EdgeListNode *parentLinkHead = NULL;
	EdgeListNode *linkHead = a->getListHead();
	bool foundHead = false;
	while (linkHead != NULL && !foundHead){
		if (linkHead->getEndVertex() == b)
			foundHead = true;
		else{
			parentLinkHead = linkHead;
			linkHead = linkHead->getNextNode();
		}
	}
	Assert(foundHead);
	if (foundHead){
		VertexListNode *loop = linkHead->getLink();
		int startnum = 0;
		while (loop){
			loop = loop->getNextNode();
			startnum++;
		}

		VertexListNode *head = linkHead->getLink();
		VertexListNode *grandparent = head;
		VertexListNode *parent = grandparent->getNextNode();
		VertexListNode *child = parent->getNextNode();

		if (grandparent->getVertex() == c && parent->getVertex() == d){
			//tail case
			if (grandparent->isPreFaceDeleted()){
				nodePool->Dealloc(grandparent);
				//check the node behide parent
				if (child && child->isPreFaceDeleted()){
					head = child;
					linkHead->setLink(head);
					nodePool->Dealloc(parent);
				}
				else{
					parent->setDeletedMark();
					head = parent;
					linkHead->setLink(head);
				}
			}
			//head case
			else{
				parent->setDeletedMark();
				//check the node behide parent
				if (child && child->isPreFaceDeleted()){
					grandparent->setNextNode(child);
					nodePool->Dealloc(parent);
				}
			}
		}
		else{
			bool found = false;
			while (child != NULL){
				if (parent->getVertex() == c && child->getVertex() == d){
					found = true;
					child->setDeletedMark();
					VertexListNode *grandchild = child->getNextNode();
					//head and tail case
					if (grandchild == NULL && head->isPreFaceDeleted()){
						if (parent->isPreFaceDeleted()){
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
					else if (parent->isPreFaceDeleted()){
						if (grandchild && grandchild->isPreFaceDeleted()){
							grandparent->setNextNode(grandchild);
							nodePool->Dealloc(child);
							nodePool->Dealloc(parent);
							child = NULL;
							parent = NULL;
						}
						else{
							grandparent->setNextNode(child);
							nodePool->Dealloc(parent);
							parent = NULL;
						}
					}
					//child and grandchild are all deletd
					else if (grandchild && grandchild->isPreFaceDeleted()){
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

			if (!found){
				if (parent->getVertex() == c && head->getVertex() == d){
					//check the tail
					if (parent->isPreFaceDeleted()){
						grandparent->setNextNode(NULL);
						nodePool->Dealloc(parent);
					}
					//the node behide the head deleted
					VertexListNode *next = head->getNextNode();
					if (next && next->isPreFaceDeleted()){
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
		if (head->getNextNode() == NULL){
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

bool DelMesher::findSubPolygons(const Face& f) const{
	Vertex *x = NULL;
	Adjacent(Segment(f.v[0], f.v[1]), &x);
	return x == f.v[2];
}

bool DelMesher::findTet(const Tetrahedron& t) const{
	Vertex *x = NULL;
	Adjacent(Face(t.v[0], t.v[1], t.v[3]), &x);
	return x == t.v[2];
}


Reference<Mesh> DelMesher::generateMesh(){
	tobeDeletedFaces.reserve(16);
	newSegsOfFaces.reserve(32);
	tobeDeletedTets.reserve(16);
	newFacesOfTets.reserve(32);

start:
	while (!mayEncroachedSegs.empty()){
		Segment s = mayEncroachedSegs.front();
		if (!Encroached(s) || !findSegment(s)){
			mayEncroachedSegs.pop_front();
			continue;
		}
		splitSubSegment(s);
		mayEncroachedSegs.pop_front();
	}
	while (!mayEncroachedFaces.empty()){
		Face f = mayEncroachedFaces.front();
		if (!findSubPolygons(f) || !Encroached(f)){
			mayEncroachedFaces.pop_front();
			continue;
		}
		splitSubPolygon(f);
		if (!mayEncroachedSegs.empty())
			goto start;
		mayEncroachedFaces.pop_front();
	}
	while (!skinnyTets.empty()){
		Tetrahedron t = skinnyTets.top();
		if (!findTet(t)){
			skinnyTets.pop();
			continue;
		}
		splitTetrahedron(t);
		if (!mayEncroachedSegs.empty() || !mayEncroachedFaces.empty())
			goto start;
		skinnyTets.pop();
	}

	//output to mesh
	std::set<Tetrahedron, tet_compare> tets;
	for (auto vert : vertices){
		if (vert->hasList()){
			EdgeListNode *linkHead = vert->getListHead();
			while (linkHead){
				Vertex *endVert = linkHead->getEndVertex();
				VertexListNode *head = linkHead->getLink();
				VertexListNode *parentNode = head;
				VertexListNode *node = parentNode->getNextNode();
				while (node){
					if (!node->isPreFaceDeleted()){
						Tetrahedron t(vert, endVert, parentNode->getVertex(), node->getVertex(), true);
						if (t.v[0] != ghost)
							tets.insert(t);
					}
					parentNode = node;
					node = node->getNextNode();
				}
				if (!head->isPreFaceDeleted()){
					Tetrahedron t(vert, endVert, parentNode->getVertex(), head->getVertex(), true);
					if (t.v[0] != ghost)
						tets.insert(t);
				}
				linkHead = linkHead->getNextNode();
			}
		}
	}

	std::set<Face, face_compare> polygons;
	for (auto entry : polygonTopology){
		Vertex *center = entry.first;
		if (entry.second){
			VertexListNode *parent = entry.second;
			VertexListNode *child = parent->getNextNode();
			while (child != NULL){
				Vertex *b = parent->getVertex();
				Vertex *c = child->getVertex();
				if (!child->isPreFaceDeleted()){
					polygons.insert(Face(center, b, c, true));
				}
				parent = child;
				child = child->getNextNode();
			}
			if (!entry.second->isPreFaceDeleted()){
				Vertex *b = parent->getVertex();
				Vertex *c = entry.second->getVertex();
				polygons.insert(Face(center, b, c, true));
			}
		}
	}	

	std::map<Vertex*, int> vi;
	for (int i = 0; i < vertices.size(); i++){
		vi[vertices[i]] = i;
	}

	int maxGap = INT_MIN;
	constexpr int numTetNode = 4;
	int vertCount = vertices.size();
	Mesh *mesh = new TetMesh(vertCount, tets.size(), polygons.size());

	int tetNodeIndex[numTetNode];
	int i = 0;
	for (auto t : tets){
		tetNodeIndex[0] = vi[t.v[0]];
		tetNodeIndex[1] = vi[t.v[1]];
		tetNodeIndex[2] = vi[t.v[2]];
		tetNodeIndex[3] = vi[t.v[3]];
		mesh->setElement(i++, tetNodeIndex);
	}
	int *vertIndices = allocAligned<int>(vertCount);
	MeshRelabeler labler(vertCount);
	labler.getNewLables(vertIndices, *mesh);

	for (int i = 0; i < vertCount; i++){
		Vertex *vertex = vertices[i];
		mesh->setVertex(vertIndices[i], Vector{ float(vertex->vert.x), float(vertex->vert.y), float(vertex->vert.z) });
	}

	i = 0;
	int surfVertIndices[3];
	for (auto f : polygons){
		surfVertIndices[0] = vertIndices[vi[f.v[0]]];
		surfVertIndices[1] = vertIndices[vi[f.v[1]]];
		surfVertIndices[2] = vertIndices[vi[f.v[2]]];
		mesh->setFacet(i++, surfVertIndices);
	}

	freeAligned(vertIndices);

	return mesh;
}
}

