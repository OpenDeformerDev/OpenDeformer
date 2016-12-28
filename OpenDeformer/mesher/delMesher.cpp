#include "stdafx.h"
#include "delMesher.h"
#include "geometer.h"
#include <memory>
#include <set>
#include <random>
#include <stack>

namespace ODER{

Predicator<REAL> DelTriangulator::predicator;

void DelTriangulator::generateSubPolygons(Vertex **vertices, int *segments, int vertexCount, int segmentCount, 
	const Face& ref, bool boundaryOnly, std::vector<Face>& output){
	meshRep.Clear();
	meshRep.Reserve(vertexCount);
	this->segments.clear();

	Vertex *v0, *v1, *v2, *v3;

#ifdef ODER_DEBUG
	if (vertexCount < 3)
		Severe("a subpolygon contains less than 3 vertices");
#endif

	//alloc ghost veritces
	ghost = meshRep.getGhostVertex();

	//construct vertices map
	for (int i = 0; i < vertexCount; i++) {
		Vertex *oriVert = vertices[i];
		Vertex *newVert = meshRep.allocVertex(oriVert->vert, oriVert->weight);
		newVert->setVertexPointer(oriVert->getPointedVertex());
		oriVert->setVertexPointer(newVert);
	}

	//shuffle input segments
	std::random_device rng;
	std::default_random_engine randomEngine(rng());
	for (int i = segmentCount - 1; i > 0; i--) {
		std::uniform_int_distribution<int> indexDistrubution(0, i);
		int newIndex = indexDistrubution(randomEngine);
		std::swap(segments[i * 2 + 0], segments[newIndex * 2 + 0]);
		std::swap(segments[i * 2 + 1], segments[newIndex * 2 + 1]);
	}

	for (int i = 0; i < segmentCount; i++)
		this->segments.insert(Segment(vertices[segments[i * 2 + 0]]->getPointedVertex(), vertices[segments[i * 2 + 1]]->getPointedVertex()));

	ghost->vert = calculateAbovePoint(vertexCount, vertices, ref);
	v0 = vertices[0]->getPointedVertex();
	v1 = vertices[1]->getPointedVertex();
	v2 = vertices[2]->getPointedVertex();

	vertices[0]->setVertexPointer(v0->getPointedVertex());
	v0->setVertexPointer(vertices[0]);
	vertices[1]->setVertexPointer(v1->getPointedVertex());
	v1->setVertexPointer(vertices[1]);
	vertices[2]->setVertexPointer(v2->getPointedVertex());
	v2->setVertexPointer(vertices[2]);

	Face f(v0, v1, v2);
	//start
	meshRep.addTriangle(v0, v1, v2);
	meshRep.addTriangle(v1, v0, ghost);
	meshRep.addTriangle(v2, v1, ghost);
	meshRep.addTriangle(v0, v2, ghost);
	for (int i = 3; i < vertexCount; i++){
		v3 = vertices[i]->getPointedVertex();
		f = findPosition(v3, f);
		v0 = f.v[0]; v1 = f.v[1]; v2 = f.v[2];
		meshRep.deleteTriangle(v0, v1, v2);
		digCavity(v3, Segment(v1, v0), &f, 0);
		digCavity(v3, Segment(v0, v2), NULL, 0);
		digCavity(v3, Segment(v2, v1), NULL, 0);


		vertices[i]->setVertexPointer(v3->getPointedVertex());
		v3->setVertexPointer(vertices[i]);
	}

	//insert segments
	for (auto s : this->segments) {
		if (!meshRep.Contain(s))
			insertSegment(s);
	}

	//clean out of region triangles
	if (boundaryOnly) {
		for (auto s : this->segments) {
			std::swap(s.v[0], s.v[1]);
			propagateClean(s, 0);
		}
		for (auto f : meshRep)
			output.push_back(Face(f.v[0]->getPointedVertex(), f.v[1]->getPointedVertex(), f.v[2]->getPointedVertex()));
	}
	else {
		for (auto f : meshRep) {
			if (!f.v[0]->isGhost())
				output.push_back(Face(f.v[0]->getPointedVertex(), f.v[1]->getPointedVertex(), f.v[2]->getPointedVertex()));
		}
	}
}

void DelTriangulator::generateSubPolygons(Vertex **vertices, Segment *segments, int vertexCount, int segmentCount, 
	const Face& ref, bool boundaryOnly, std::vector<Face>& output) {
	meshRep.Clear();
	meshRep.Reserve(vertexCount);
	this->segments.clear();

	Vertex *v0, *v1, *v2, *v3;

#ifdef ODER_DEBUG
	if (vertexCount < 3)
		Severe("a subpolygon contains less than 3 vertices");
#endif

	//alloc ghost veritces
	ghost = meshRep.getGhostVertex();

	for (int i = 0; i < vertexCount; i++) {
		Vertex *oriVert = vertices[i];
		Vertex *newVert = meshRep.allocVertex(oriVert->vert, oriVert->weight);
		newVert->setVertexPointer(oriVert->getPointedVertex());
		oriVert->setVertexPointer(newVert);
	}

	//shuffle input segments
	std::random_device rng;
	std::default_random_engine randomEngine(rng());
	std::shuffle(segments, segments + segmentCount, randomEngine);

	for (int i = 0; i < segmentCount; i++) {
		Segment seg = segments[i];
		this->segments.insert(Segment(seg.v[0]->getPointedVertex(), seg.v[1]->getPointedVertex()));
	}

	ghost->vert = calculateAbovePoint(vertexCount, vertices, ref);
	v0 = vertices[0]->getPointedVertex();
	v1 = vertices[1]->getPointedVertex();
	v2 = vertices[2]->getPointedVertex();

	vertices[0]->setVertexPointer(v0->getPointedVertex());
	v0->setVertexPointer(vertices[0]);
	vertices[1]->setVertexPointer(v1->getPointedVertex());
	v1->setVertexPointer(vertices[1]);
	vertices[2]->setVertexPointer(v2->getPointedVertex());
	v2->setVertexPointer(vertices[2]);

	Face f(v0, v1, v2);
	//start
	meshRep.addTriangle(v0, v1, v2);
	meshRep.addTriangle(v1, v0, ghost);
	meshRep.addTriangle(v2, v1, ghost);
	meshRep.addTriangle(v0, v2, ghost);
	for (int i = 3; i < vertexCount; i++) {
		v3 = vertices[i]->getPointedVertex();
		f = findPosition(v3, f);
		v0 = f.v[0]; v1 = f.v[1]; v2 = f.v[2];
		meshRep.deleteTriangle(v0, v1, v2);
		digCavity(v3, Segment(v1, v0), &f, 0);
		digCavity(v3, Segment(v0, v2), NULL, 0);
		digCavity(v3, Segment(v2, v1), NULL, 0);

		vertices[i]->setVertexPointer(v3->getPointedVertex());
		v3->setVertexPointer(vertices[i]);
	}

	//insert segments
	for (auto s : this->segments) {
		if (!meshRep.Contain(s))
			insertSegment(s);
	}

	//clean out of region triangles
	if (boundaryOnly) {
		for (auto s : this->segments) {
			std::swap(s.v[0], s.v[1]);
			propagateClean(s, 0);
		}
		for (auto f : meshRep) 
			output.push_back(Face(f.v[0]->getPointedVertex(), f.v[1]->getPointedVertex(), f.v[2]->getPointedVertex()));
	}
	else {
		for (auto f : meshRep) {
			if (!f.v[0]->isGhost())
				output.push_back(Face(f.v[0]->getPointedVertex(), f.v[1]->getPointedVertex(), f.v[2]->getPointedVertex()));
		}
	}
}

void DelTriangulator::insertSegments(const Face *triangles, const Segment *segments, int triangleCount, int segmentCount, std::vector<Face>& output) {
	meshRep.Clear();
	meshRep.Reserve(triangleCount);
	this->segments.clear();

	//clean pointers
	std::vector<Vertex *> oldPointers;
	oldPointers.reserve(3 * triangleCount);
	for (int i = 0; i < triangleCount; i++) {
		Face f = triangles[i];
		for (int j = 0; j < 3; j++) {
			oldPointers.push_back(f.v[j]->getPointedVertex());
			f.v[j]->setVertexPointer(NULL);
		}
	}

	REAL maxArea2 = 0;
	Face largetsFace;
	for (int i = 0; i < triangleCount; i++) {
		Face f = triangles[i];
		for (int j = 0; j < 3; j++) {
			if (f.v[j]->getPointedVertex() == NULL) {
				Vertex *newVert = meshRep.allocVertex(f.v[j]->vert, f.v[j]->weight);
				f.v[j]->setVertexPointer(newVert);
				newVert->setVertexPointer(f.v[j]);
			}
		}
		meshRep.addTriangle(f.v[0]->getPointedVertex(), f.v[1]->getPointedVertex(), f.v[2]->getPointedVertex());
		REAL area2 = (Geometer::triangleNormal(f.v[0]->vert, f.v[1]->vert, f.v[2]->vert)).length2();
		if (area2 > maxArea2) {
			maxArea2 = area2;
			largetsFace = f;
		}
	}
	//calculate above point
	ghost = meshRep.getGhostVertex();
	ghost->vert = Geometer::calculateAbovePoint(largetsFace.v[0]->vert, largetsFace.v[1]->vert, largetsFace.v[2]->vert);

	for (int i = 0; i < segmentCount; i++) {
		Segment oriSeg = segments[i];
		Segment newSeg(oriSeg.v[0]->getPointedVertex(), oriSeg.v[1]->getPointedVertex());
		this->segments.insert(newSeg);
		if (!meshRep.Contain(newSeg))
			insertSegment(newSeg);
	}

	//restore pointers
	for (int i = 0; i < triangleCount; i++) {
		Face f = triangles[i];
		for (int j = 0; j < 3; j++) 
			f.v[j]->setVertexPointer(oldPointers[i * 3 + j]);
	}

	for (auto f : meshRep)
		output.push_back(Face(f.v[0]->getPointedVertex(), f.v[1]->getPointedVertex(), f.v[2]->getPointedVertex()));
}

void DelTriangulator::propagateClean(const Segment& s, int depth) {
	Vertex *w;
	if (!meshRep.Adjacent(s, &w)) return;

	Vertex *u = s.v[0], *v = s.v[1];

	meshRep.deleteTriangle(w, u, v);

	Segment uw(u, w), wv(w, v);
	if (segments.find(uw) == segments.end()) propagateClean(uw, depth + 1);
	if (segments.find(wv) == segments.end()) propagateClean(wv, depth + 1);
}

Face DelTriangulator::findPosition(Vertex *u, const Face &f) const{
	size_t maxIter = std::numeric_limits<size_t>::max();
	Vertex *a = f.v[0], *b = f.v[1], *c = f.v[2];
	//ghost triangle
	if (a->isGhost()){
		if (predicator.inHalfSpace2D(u->vert, b->vert, c->vert, ghost->vert))
			return f;
		Vertex *x = ghost;
		meshRep.Adjacent(Segment(c, b), &x);
		//Face(x, c, b)
		a = x; std::swap(b, c);
	}

	Assert(!a->isGhost());
	//normal triangle

	while (maxIter--){
		//caution:because NULL->ghost,can't use NULL,
		//so it's might case memory corrupt above
		Vertex *x;
		bool ori0 = predicator.orient2d(b->vert, a->vert, u->vert, ghost->vert) > 0;
		bool ori1 = predicator.orient2d(a->vert, c->vert, u->vert, ghost->vert) > 0;

		int condition = ori0 + (ori1 << 1);
		switch (condition){
		case 0:
			if (predicator.orient2d(c->vert, b->vert, u->vert, ghost->vert) > 0){
				meshRep.Adjacent(Segment(c, b), &x);
				if (x->isGhost())
					return Face(ghost, c, b);
				a = x; std::swap(b, c); //Face(x, c, b)
			}
			else
				return Face(a, b, c, true);
			break;
		case 1:
			meshRep.Adjacent(Segment(b, a), &x);
			if (x->isGhost())
				return Face(ghost, b, a);
			c = a; a = x; //Face(x, b, a)
			break;
		case 2:
			meshRep.Adjacent(Segment(a, c), &x);
			if (x->isGhost())
				return Face(ghost, a, c);
			b = a; a = x; //Face(x, a, c)
			break;
		case 3:
			if (Randomnation(2)){
				meshRep.Adjacent(Segment(b, a), &x);
				if (x->isGhost())
					return Face(ghost, b, a);
				c = a; a = x; //Face(x, b, a)
			}
			else{
				meshRep.Adjacent(Segment(a, c), &x);
				if (x->isGhost())
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

void DelTriangulator::insertSegment(const Segment& s) {
	std::vector<Vertex *> positive, negtive;
	findCavity(s, positive, negtive);
	triangulateHalfHole(positive);
	triangulateHalfHole(negtive);
}

void DelTriangulator::findCavity(const Segment& s, std::vector<Vertex *>& positive, std::vector<Vertex *>& negtive) {
	Vertex *a = s.v[0], *b = s.v[1];
	if (Randomnation(2)) std::swap(a, b);
	DelVector aa = a->vert, bb = b->vert, above = ghost->vert;
	Face f;

	positive.push_back(a);
	negtive.push_back(a);

	bool find = meshRep.findIntersectedFace(a, bb, &f);
	Assert(find);

	Vertex *c = f.v[1], *d = f.v[2];
	meshRep.deleteTriangle(a, c, d);
	negtive.push_back(c);
	positive.push_back(d);

	//transverse the intersected triangle
	Vertex *opposite = NULL;
	meshRep.Adjacent(Segment(d, c), &opposite);
	meshRep.deleteTriangle(d, c, opposite);

	while (opposite != b) {
		REAL ori = predicator.orient2d(aa, bb, opposite->vert, above);
		Assert(ori != 0);
		if (ori > 0) {
			positive.push_back(opposite);
			d = opposite;
		}
		else {
			negtive.push_back(opposite);
			c = opposite;
		}
		meshRep.Adjacent(Segment(d, c), &opposite);
		meshRep.deleteTriangle(d, c, opposite);
	}

	positive.push_back(b);
	negtive.push_back(b);

	std::reverse(positive.begin(), positive.end());
}

void DelTriangulator::triangulateHalfHole(const std::vector<Vertex *>& vertices){
	int size = vertices.size();
	int *prev = new int[size];
	int *next = new int[size];
	REAL *distance = new REAL[size];
	cavityRep.Clear();

	//costruct vertices map
	std::vector<Vertex *> cavityVertices;
	cavityVertices.reserve(size);
	for (auto vert : vertices) {
		Vertex *cavityVert = cavityRep.allocVertex(vert->vert, vert->weight);
		cavityVert->setVertexPointer(vert);
		cavityVertices.push_back(cavityVert);
	}

	Vertex *start = cavityVertices[size - 1], *end = cavityVertices[0];
	prev[0] = size - 1; prev[size - 1] = size - 2;
	next[0] = 1; next[size - 1] = 0;
	distance[0] = distance[size - 1] = 0;

	DelVector ab = end->vert - start->vert;
	for (int i = 1; i < size - 1; i++) {
		prev[i] = i - 1;
		next[i] = i + 1;
		distance[i] = ((cavityVertices[i]->vert - start->vert) % ab).length2();
	}

	int indexCount = size - 2;
	int *indices = new int[indexCount];
	for (int i = 0; i < indexCount; i++)
		indices[i] = i + 1;
	std::random_device rng;
	std::default_random_engine randomEngine(rng());
	std::shuffle(indices, indices + indexCount, randomEngine);
	
	for (int i = indexCount - 1; i > 0; i--) {
		std::uniform_int_distribution<int> indexDistrubution(0, i);
		int j = 0, trueIndex = 0;
		do {
			j = indexDistrubution(randomEngine);
			trueIndex = indices[j];
		}while (distance[trueIndex] < distance[prev[trueIndex]] &&
			distance[trueIndex] < distance[next[trueIndex]]);
		std::swap(indices[i], indices[j]);
		next[prev[trueIndex]] = next[trueIndex];
		prev[next[trueIndex]] = prev[trueIndex];
	}
	cavityRep.addTriangle(start, end, cavityVertices[indices[0]]);


	std::deque<Vertex *> fanVertices;
	std::vector<Vertex *> convexPoly;
	for (int i = 1; i < indexCount; i++) {
		int trueIndex = indices[i];
		bool convexPolyOri = true;
		insertVertexToCavity(cavityVertices[trueIndex], cavityVertices[next[trueIndex]], cavityVertices[prev[trueIndex]],
			false, true, fanVertices, 0);

		Assert(fanVertices.size() != 1);
		if (fanVertices.size() > 2) {
			convexPoly.reserve(fanVertices.size());
			for (auto vert : fanVertices)
				convexPoly.push_back(vert);
			triangulateConvexPoly(cavityVertices[trueIndex], convexPoly);

			fanVertices.clear();
			convexPoly.clear();
		}
	}

	for (auto f : cavityRep)
		meshRep.addTriangle(f.v[0]->getPointedVertex(), f.v[1]->getPointedVertex(), f.v[2]->getPointedVertex());

	delete[] indices;
	delete[] prev;
	delete[] next;
	delete[] distance;
}

void DelTriangulator::insertVertexToCavity(Vertex *u, Vertex *v, Vertex *w, 
	bool mark, bool oriTest, std::deque<Vertex *>& marked, int depth) {
	Vertex *x = NULL;
	bool exited = cavityRep.Adjacent(Segment(w, v), &x);
	bool outCircle = true;
	bool deleted = false;
	bool negativeOri = false;
	if (exited) {
		outCircle = predicator.inOrthoCirclePerturbed(u->vert, u->weight, v->vert, v->weight, w->vert, w->weight,
			x->vert, x->weight, ghost->vert) > 0;

		if (oriTest) negativeOri = predicator.orient2d(u->vert, v->vert, w->vert, ghost->vert) <= 0;
		deleted = outCircle | negativeOri;
	}

	if (deleted) {
		cavityRep.deleteTriangle(x, w, v);
		bool nextMark = outCircle ? false : true;
		insertVertexToCavity(u, v, x, nextMark, negativeOri, marked, depth + 1);
		insertVertexToCavity(u, x, w, nextMark, negativeOri, marked, depth + 1);
	}
	else {
		if (mark) {
			if (marked.empty()) {
				marked.push_back(v);
				marked.push_back(w);
			}
			else {
				if (marked.back() == v)
					marked.push_back(w);
				else {
					Assert(marked.front() == w);
					marked.push_front(v);
				}
			}
		}
		else
			cavityRep.addTriangle(u, v, w);
	}
}

void DelTriangulator::insertVertexToCavity(Vertex *u, Vertex *v, Vertex *w, int depth) {
	Vertex *x = NULL;
	bool exited = cavityRep.Adjacent(Segment(w, v), &x);
	bool deleted = false;

	if (exited) {
		deleted = predicator.inOrthoCirclePerturbed(u->vert, u->weight, v->vert, v->weight,
			w->vert, w->weight, x->vert, x->weight, ghost->vert) > 0;
	}

	if (deleted) {
		cavityRep.deleteTriangle(x, w, v);
		insertVertexToCavity(u, v, x, depth + 1);
		insertVertexToCavity(u, x, w, depth + 1);
	}
	else
		cavityRep.addTriangle(u, v, w);
}

void DelTriangulator::triangulateConvexPoly(Vertex *u, const std::vector<Vertex *>& convexPoly) {
	int size = convexPoly.size() - 2;
	int *indices = new int[size];
	int *prev = new int[size + 2];
	int *next = new int[size + 2];
	bool *marked = new bool[size + 2];

	//isolate convex poly
	u->setMark();
	for (auto v : convexPoly) v->setMark();

	prev[0] = -1; prev[size + 1] = size;
	next[0] = 1; next[size + 1] = -1;
	for (int i = 0; i < size; i++) {
		prev[i + 1] = i;
		next[i + 1] = i + 2;
		indices[i] = i + 1;
	}
	std::random_device rng;
	std::default_random_engine randomEngine(rng());
	std::shuffle(indices, indices + size, randomEngine);

	for (int i = size - 1; i >= 0; i--) {
		int trueIndex = indices[i];
		next[prev[trueIndex]] = next[trueIndex];
		prev[next[trueIndex]] = prev[trueIndex];
	}
	cavityRep.addTriangle(u, convexPoly[0], convexPoly[size + 1]);

	for (int i = 0; i < size; i++) {
		int trueIndex = indices[i];
		insertVertexToCavity(convexPoly[trueIndex], convexPoly[next[trueIndex]], convexPoly[prev[trueIndex]], 0);
	}

	//restore canvity
	u->unSetMark();
	for (auto v : convexPoly) v->unSetMark();

	delete[] prev;
	delete[] next;
	delete[] indices;
}

void DelTriangulator::digCavity(Vertex *u, const Segment &s, Face *rf, int depth){
	Vertex *x;
	if (!meshRep.Adjacent(s, &x))
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
		deleted = predicator.inHalfSpace2D(u->vert, w->vert, v->vert, ghost->vert);
	}
	else{
		deleted = predicator.inOrthoCirclePerturbed(u->vert, u->weight, verts[1]->vert, verts[1]->weight, 
			                               verts[0]->vert, verts[0]->weight, verts[2]->vert, verts[2]->weight, ghost->vert) > 0;
	}
	if (deleted){
		meshRep.deleteTriangle(verts[0], verts[1], verts[2]);
		digCavity(u, Segment(verts[2], verts[1]), rf, depth + 1);
		digCavity(u, Segment(verts[0], verts[2]), rf, depth + 1);
	}
	else{
		meshRep.addTriangle(u, verts[1], verts[0]);
		if (rf)
			*rf = Face(u, verts[1], verts[0], true);
	}
}

DelVector DelTriangulator::calculateAbovePoint(int vertexCount, Vertex** vertices, const Face& ref) {
	Vertex *a = vertices[0];
	Vertex *b = NULL, *c = NULL;
	int bIndex = 0, cIndex = 0;

	REAL longest = 0;
	for (int i = 1; i < vertexCount; i++) {
		REAL length = (vertices[i]->vert - a->vert).length2();
		if (length > longest) {
			b = vertices[i];
			bIndex = i;
			longest = length;
		}
	}
	std::swap(vertices[1], vertices[bIndex]);

	DelVector vec = b->vert - a->vert;
	REAL largest = 0;
	for (int i = 2; i < vertexCount; i++) {
		REAL area = (vec % (vertices[i]->vert - a->vert)).length2();
		if (area > largest) {
			c = vertices[i];
			cIndex = i;
			largest = area;
		}
	}
	std::swap(vertices[2], vertices[cIndex]);

	DelVector n = sqrt(longest) * REAL(0.5) * Normalize(vec % (c->vert - a->vert));
	DelVector above = a->vert + n;
	if (predicator.orient3d(above, ref.v[0]->vert, ref.v[1]->vert, ref.v[2]->vert) < 0) {
		above = a->vert - n;
		std::swap(vertices[1], vertices[2]);
	}

	return above;
}

Predicator<REAL> DelMesher::predicator;

DelMesher::DelMesher(Vector *surfvs, int *triangles, int numv, int numtri, REAL maxR, REAL maxRa, REAL facetAngleTol) {
	maxRatio = maxRa;
	maxRadius = maxR;

	REAL maxREration = maxRatio;
	auto compare = [maxREration](const Tetrahedron& left, const Tetrahedron& right) {
		if (left.getREration() <= maxREration && right.getREration() <= maxREration)
			return left.getRadius() < right.getRadius();
		return left.getREration() < right.getREration();
	};
	skinnyTets = std::priority_queue<Tetrahedron, std::vector<Tetrahedron>, 
		std::function<bool(const Tetrahedron&, const Tetrahedron&)>>(compare, std::vector<Tetrahedron>());

	oriVertices.reserve(numv);
	TriMeshDataStructure oriSurface;
	for (int i = 0; i < numv; i++) {
		boundBox.Insert(surfvs[i]);
		Vertex *vert = oriSurface.allocVertex(DelVector{ surfvs[i].x, surfvs[i].y, surfvs[i].z }, REAL(0));
		oriVertices.push_back(vert);
	}

	std::vector<Face> oriTriangles;
	oriTriangles.reserve(numtri);
	for (int i = 0; i < numtri; i++) {
		oriSurface.addTriangle(oriVertices[triangles[3 * i]], oriVertices[triangles[3 * i + 1]], oriVertices[triangles[3 * i + 2]]);
		oriTriangles.push_back(Face(oriVertices[triangles[3 * i]], oriVertices[triangles[3 * i + 1]], oriVertices[triangles[3 * i + 2]]));
	}

	REAL facetRadianTol = facetAngleTol / REAL(180) * REAL(M_PI);
	std::vector<Segment> boundaries;
	std::vector<Vertex *> coplanarVertices;
	DelTriangulator triangulator;
	std::vector<Face> faces;

	int faceIndex = 0;
	for (auto f : oriTriangles) {
		if (oriSurface.Contain(f)) {
			detectCoplanarFaces(f, facetRadianTol, coplanarVertices, boundaries, oriSurface);
			for (auto s : boundaries) {
				for (int i = 0; i < 2; i++) {
					if (s.v[i]->getPointedVertex() == NULL) {
						Vertex *newVert = meshRep.allocVertex(s.v[i]->vert, s.v[i]->weight, 
							VertexType(Vertex_Fixed | Vertex_LowDimension));
						surfaceRep.bindVolumeVertex(newVert);
						s.v[i]->setVertexPointer(newVert);
					}
				}
				meshRep.addSegment(s.v[0]->getPointedVertex(), s.v[1]->getPointedVertex());
			}

			for (auto v : coplanarVertices) {
				if (v->getPointedVertex() == NULL) {
					Vertex *newVert = meshRep.allocVertex(v->vert, v->weight, Vertex_FixedFacet);
					surfaceRep.bindVolumeVertex(newVert);
					v->setVertexPointer(newVert);
				}
			}

			//re-triangulation
			if (coplanarVertices.size() > 3) {
				triangulator.generateSubPolygons(&coplanarVertices[0], &boundaries[0], coplanarVertices.size(), boundaries.size(), f, true, faces);
				for (auto face : faces) 
					surfaceRep.addTriangle(face.v[0]->getPointedVertex(), face.v[1]->getPointedVertex(), face.v[2]->getPointedVertex(), faceIndex);
			}
			else
				surfaceRep.addTriangle(f.v[0]->getPointedVertex(), f.v[1]->getPointedVertex(), f.v[2]->getPointedVertex(), faceIndex);

			//store vertices per polygon
			size_t coplanarVerticesCount = coplanarVertices.size();
			uintptr_t *verticesArray = new uintptr_t[coplanarVerticesCount + 1];
			verticesArray[0] = (uintptr_t)coplanarVerticesCount;
			for (size_t i = 0; i < coplanarVerticesCount; i++)
				verticesArray[i + 1] = (uintptr_t)coplanarVertices[i]->getPointedVertex();

			verticesPerPolygon.push_back(verticesArray);
			coplanarVertices.clear();
			faces.clear();
			boundaries.clear();
			faceIndex += 1;
		}
	}
	for (auto &v : oriVertices) v = v->getPointedVertex();

	detectAcuteVertices();

	//allocate space
	tobeDeletedFaces.reserve(16);
	newSegsOfFaces.reserve(32);
	tobeDeletedTets.reserve(16);
	newTets.reserve(32);

	constrainedTriangulation();
}

DelMesher::~DelMesher() {
	for (auto p : verticesPerPolygon) delete[] p;
}

void DelMesher::triangulation3D(std::vector<Vertex *>& vertices, TetMeshDataStructure& meshRep, bool insertToSkinny){
	Vertex *a = vertices[0], *b = vertices[1];
	Vertex *c = NULL, *d = NULL;
	for (std::vector<Vertex *>::size_type i = 2; i < vertices.size(); i++) {
		c = vertices[i];
		if (!predicator.fastCoLine(a->vert, b->vert, c->vert)) {
			std::swap(vertices[2], vertices[i]);
			break;
		}
	}
	for (std::vector<Vertex *>::size_type i = 3; i < vertices.size(); i++){
		d = vertices[i];
		if (!predicator.fastCoPlane(a->vert, b->vert, c->vert, d->vert)){
			std::swap(vertices[3], vertices[i]);
			break;
		}
	}
	if (predicator.orient3d(a->vert, b->vert, c->vert, d->vert) < 0)
		std::swap(a, b);

	Tetrahedron t(a, b, c, d, true);
	if (insertToSkinny) {
		t.setRationAndRadius();
		if (t.getREration() > maxRatio || t.getRadius() > maxRadius)
			skinnyTets.push(t);
	}

	Vertex *ghost = meshRep.getGhostVertex();

	meshRep.addTetrahedron(a, b, c, d);
	meshRep.addTetrahedron(ghost, b, d, c);
	meshRep.addTetrahedron(ghost, a, c, d);
	meshRep.addTetrahedron(ghost, a, b, c);
	meshRep.addTetrahedron(ghost, a, d, b);

	REAL vvv = predicator.orient3d(a->vert, b->vert, c->vert, d->vert);
	VolumeVertexInsertionFlags vifs;
	vifs.skinnyTetTest = insertToSkinny;
	Tetrahedron rt = t;
	for (std::vector<Vertex *>::size_type i = 4; i < vertices.size(); i++) {
		Vertex *v = vertices[i];
		t = findPosition(v, rt, meshRep);
		insertVertex(v, t, meshRep, vifs, &rt);
	}
}

void DelMesher::constrainedTriangulation() {
	triangulation3D(oriVertices, meshRep, false);
	for (auto iter = meshRep.segmentBegin(); iter != meshRep.segmentEnd(); ++iter) mayMissingSegs.push_back(*iter);
	std::sort(mayMissingSegs.begin(), mayMissingSegs.end(), 
		[](const Segment& left, const Segment& right) {
		return (left.v[0]->vert - left.v[1]->vert).length2() < (right.v[0]->vert - right.v[1]->vert).length2();
	});
	segmentsRecovery(false);
	for (auto f : surfaceRep) mayMissingFaces.push_back(f);
	facesRecovery();
}

void DelMesher::detectAcuteVertices() const {
	constexpr int maxCount = 5;
	constexpr REAL acuteRadian = REAL(89) / REAL(180) * REAL(M_PI);
	std::unordered_map<Vertex *, std::vector<Vertex *>> corners;

	for (auto iter = meshRep.segmentBegin(); iter != meshRep.segmentEnd(); ++iter) {
		Segment seg = *iter;
		for (int i = 0; i < 2; i++) {
			Vertex *a = seg.v[i], *b = seg.v[(i + 1) % 2];
			if (!matchVertexFlag(a->getVertexType(), VertexType::Vertex_Acute)) {
				auto found = corners.find(a);
				if (found != corners.end()) {
					const auto &endVertices = found->second;
					bool detected = false;
					for (auto v : endVertices) {
						if (Geometer::interiorAngle(a->vert, b->vert, v->vert) < acuteRadian) {
							detected = true;
							break;
						}
					}
					if (detected) {
						a->setAcute();
						corners.erase(found);
					}
					else {
						if (endVertices.size() < maxCount)
							found->second.push_back(b);
						else {
							a->setAcute();
							corners.erase(found);
						}
					}
				}
				else {
					std::vector<Vertex *> endVerts;
					endVerts.reserve(maxCount);
					endVerts.push_back(b);
					corners.insert(std::make_pair(a, std::move(endVerts)));
				}
			}
		}
	}
}

Tetrahedron DelMesher::findPosition(Vertex *u, const Tetrahedron &t, const TetMeshDataStructure& meshRep) const{
	size_t maxIter = std::numeric_limits<size_t>::max();
	Vertex *a = t.v[0], *b = t.v[1], *c = t.v[2], *d = t.v[3];
	if (b->isGhost()) {
		std::swap(a, b); std::swap(c, d);
	}
	else if (c->isGhost()) {
		std::swap(a, c); std::swap(b, d);
	}
	else if (d->isGhost()) {
		std::swap(a, d); std::swap(b, c);
	}

	while (maxIter--) {
		Vertex *x = NULL;
		if (a->isGhost()){
			if (predicator.inOrthoHalfSpace3D(u->vert, u->weight, b->vert, b->weight, c->vert, c->weight, d->vert, d->weight))
				return t;
			meshRep.Adjacent(Face(b, d, c, true), &x);
			a = x; std::swap(d, c);
		}
		else{
			bool ori0 = predicator.orient3d(u->vert, a->vert, b->vert, c->vert) > 0;
			bool ori1 = predicator.orient3d(u->vert, a->vert, c->vert, d->vert) > 0;
			bool ori2 = predicator.orient3d(u->vert, a->vert, d->vert, b->vert) > 0;

			int condition = ori0 + (ori1 << 1) + (ori2 << 2);

			switch (condition){
			case 0:
				if (predicator.orient3d(u->vert, d->vert, c->vert, b->vert) > 0){
					meshRep.Adjacent(Face(d, c, b), &x);
					if (x->isGhost())
					    return Tetrahedron(x, d, c, b, true);
					//Tetrahedron(x, d, c, b)
				    a = x; std::swap(b, d);
				}
				else{
					return Tetrahedron(a, b, c, d, true);
				}
				break;
			case 1:
				meshRep.Adjacent(Face(a, b, c), &x);
				if (x->isGhost())
					return Tetrahedron(x, a, b, c, true);
				//Tetrahedron(x, a, b, c)
				d = c; c = b; b = a; a = x;
				break;
			case 2:
				meshRep.Adjacent(Face(a, c, d), &x);
				if (x->isGhost())
					return Tetrahedron(x, a, c, d, true);
				//Tetrahedron(x, a, c, d)
				b = a; a = x;
				break;
			case 3:
				if (Randomnation(2) == 0){
					meshRep.Adjacent(Face(a, b, c), &x);
					if (x->isGhost())
						return Tetrahedron(x, a, b, c, true);
					//Tetrahedron(x, a, b, c)
					d = c; c = b; b = a; a = x;
				}
				else{
					meshRep.Adjacent(Face(a, c, d), &x);
					if (x->isGhost())
						return Tetrahedron(x, a, c, d, true);
					//Tetrahedron(x, a, c, d)
					b = a; a = x;
				}
				break;
			case 4:
				meshRep.Adjacent(Face(a, d, b), &x);
				if (x->isGhost())
					return Tetrahedron(x, a, d, b, true);
				//Tetrahedron(x, a, d, b)
				c = d; d = b; b = a; a = x;
				break;
			case 5:
				if (Randomnation(2) == 0){
					meshRep.Adjacent(Face(a, b, c), &x);
					if (x->isGhost())
						return Tetrahedron(x, a, b, c, true);
					//Tetrahedron(x, a, b, c)
					d = c; c = b; b = a; a = x;
				}
				else{
					meshRep.Adjacent(Face(a, d, b), &x);
					if (x->isGhost())
						return Tetrahedron(x, a, d, b, true);
					//Tetrahedron(x, a, d, b)
					c = d; d = b; b = a; a = x;
				}
				break;
			case 6:
				if (Randomnation(2) == 0){
					meshRep.Adjacent(Face(a, c, d), &x);
					if (x->isGhost())
						return Tetrahedron(x, a, c, d, true);
					//Tetrahedron(x, a, c, d)
					b = a; a = x;
				}
				else{
					meshRep.Adjacent(Face(a, d, b), &x);
					if (x->isGhost())
						return Tetrahedron(x, a, d, b, true);
					//Tetrahedron(x, a, d, b)
					c = d; d = b; b = a; a = x;
				}
				break;
			case 7:
			{
				unsigned int term = Randomnation(3);
				if (term == 0){
					meshRep.Adjacent(Face(a, b, c), &x);
					if (x->isGhost())
						return Tetrahedron(x, a, b, c, true);
					//Tetrahedron(x, a, b, c)
					d = c; c = b; b = a; a = x;
				}
				else if (term == 1){
					meshRep.Adjacent(Face(a, c, d), &x);
					if (x->isGhost())
						return Tetrahedron(x, a, c, d, true);
					b = a; a = x;
				}
				else{
					meshRep.Adjacent(Face(a, d, b), &x);
					if (x->isGhost())
						return Tetrahedron(x, a, d, b, true);
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

Face DelMesher::findPosition(Vertex *u, const DelVector& above, const Face& f) const{
	size_t maxIter = std::numeric_limits<size_t>::max();
	Vertex *a = f.v[0];
	Vertex *b = f.v[1], *c = f.v[2];

	while (maxIter--){
		Vertex *x = NULL;
		bool ori0 = predicator.orient2d(b->vert, a->vert, u->vert, above) > 0;
		bool ori1 = predicator.orient2d(a->vert, c->vert, u->vert, above) > 0;

		int condition = ori0 + (ori1 << 1);
		switch (condition){
		case 0:
			if (predicator.orientCoplane(c->vert, b->vert, u->vert) > 0){
				surfaceRep.Adjacent(Segment(c, b), &x);
				a = x; std::swap(b, c); //Face(x, c, b)
			}
			else{
				return Face(a, b, c, true);
			}
			break;
		case 1:
			surfaceRep.Adjacent(Segment(b, a), &x);
			c = a; a = x; //Face(x, b, a)
			break;
		case 2:
			surfaceRep.Adjacent(Segment(a, c), &x);
			b = a; a = x; //Face(x, a, c)
			break;
		case 3:
			if (Randomnation(2) == 0){
				surfaceRep.Adjacent(Segment(b, a), &x);
				c = a; a = x; //Face(x, b, a)
			}
			else{
				surfaceRep.Adjacent(Segment(a, c), &x);
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

bool DelMesher::digCavity(Vertex *u, const Face& f, TetMeshDataStructure& meshRep,
	const VolumeVertexInsertionFlags& vifs, Tetrahedron *rt){
	Vertex *a;
	if (!meshRep.Adjacent(f, &a)) return true;

	bool success = true;

	Vertex *b = f.v[0];
	Vertex *c = f.v[1];
	Vertex *d = f.v[2];

	bool deleted = false;
	bool constrained = false;
	if (vifs.cdt) constrained = surfaceRep.Contain(f) || surfaceRep.Contain(Face(b, d, c));
	if (!constrained) {
		if (a->isGhost())
			deleted = predicator.inOrthoHalfSpace3D(u->vert, u->weight, b->vert, b->weight, c->vert, c->weight, d->vert, d->weight);
		else if (b->isGhost())
			deleted = predicator.inOrthoHalfSpace3D(u->vert, u->weight, a->vert, a->weight, d->vert, d->weight, c->vert, c->weight);
		else
			deleted = predicator.inOrthoSpherePerturbed(u->vert, u->weight, b->vert, b->weight, d->vert, d->weight, c->vert, c->weight, a->vert, a->weight) > 0;
	}

	if (deleted) {
		if (!vifs.trueInsertion) tobeDeletedTets.push_back(Tetrahedron(a, b, c, d));
		meshRep.deleteTetrahedron(a, b, c, d);

		if (vifs.missingFaceTest || vifs.encroachFaceTest) {
			Face front(b, c, d), back(b, d, c);
			if (surfaceRep.Contain(front)) {
				if (vifs.missingFaceTest) mayMissingFaces.push_back(front);
				if (vifs.encroachFaceTest) mayEncroachedFaces.push_back(front);
			}
			if (surfaceRep.Contain(back)) {
				if (vifs.missingFaceTest) mayMissingFaces.push_back(back);
				if (vifs.encroachFaceTest) mayEncroachedFaces.push_back(back);
			}
		}
		if (vifs.encroachSegTest || vifs.missingSegTest) {
			Segment ab(a, b), ac(a, c), ad(a, d);
			if (meshRep.isSegment(ab) && !meshRep.testAndMark(a, b)) markedSegments.push_back(ab);
			if (meshRep.isSegment(ac) && !meshRep.testAndMark(a, c)) markedSegments.push_back(ac);
			if (meshRep.isSegment(ad) && !meshRep.testAndMark(a, d)) markedSegments.push_back(ad);
		}

		success &= digCavity(u, Face(a, c, d, true), meshRep, vifs, rt);
		success &= digCavity(u, Face(a, b, c, true), meshRep, vifs, rt);
		success &= digCavity(u, Face(a, d, b, true), meshRep, vifs, rt);
	}
	else {
		if (vifs.trueInsertion) meshRep.addTetrahedron(u, b, d, c);
		else newTets.push_back(Tetrahedron(u, b, d, c));

		if (constrained) success = predicator.orient3d(u->vert, b->vert, d->vert, c->vert) > REAL(0);
		if (rt) *rt = Tetrahedron(u, b, d, c, true);
		if (vifs.encroachFaceTest) {
			Face front(b, d, c), back(b, c, d);
			if (surfaceRep.Contain(front)) mayEncroachedFaces.push_back(front);
			if (surfaceRep.Contain(back)) mayEncroachedFaces.push_back(back);
		}
		if (vifs.skinnyTetTest){
			Tetrahedron t(u, b, d, c, true);
			t.setRationAndRadius();
			if (t.getREration() > maxRatio || t.getRadius() > maxRadius)
				skinnyTets.push(t);
		}
	}

	return success;
}

void DelMesher::digCavity(Vertex *u, const DelVector& above,const Segment &s, int index, const SurfaceVertexInsertionFlags& vifs){
	Vertex *x;
	int newIndex;
	if (!surfaceRep.Adjacent(s, &x, &newIndex))
		return;

	Vertex *w = s.v[0];
	Vertex *v = s.v[1];

	if (!findSegment(Segment(w, v, true)) &&
		predicator.inOrthoCirclePerturbed(u->vert, u->weight, v->vert, v->weight, w->vert, w->weight, x->vert, x->weight, above) > 0){
		if (!vifs.trueInsertion)
			tobeDeletedFaces.push_back(Face(w, v, x, index, true));
		surfaceRep.deleteTriangle(w, v, x);
		digCavity(u, above, Segment(x, v), newIndex, vifs);
		digCavity(u, above, Segment(w, x), newIndex, vifs);
	}
	else{
		if (vifs.trueInsertion)
			surfaceRep.addTriangle(u, v, w, index);
		else
			newSegsOfFaces.push_back(Segment(v, w));
		if (vifs.missingTest) mayMissingFaces.push_back(Face(u, v, w, index, true));
		if (vifs.encroachmentTest) mayEncroachedFaces.push_back(Face(u, v, w, index, true));
	}
}

bool DelMesher::findIntersectedTetrahedron(Vertex *a, const DelVector& bb, Tetrahedron *t) const {
	if (!meshRep.adjacent2Vertex(a, t)) return false;

	Vertex *b = t->v[1], *c = t->v[2], *d = t->v[3];
	if (b->isGhost()) {
		bool found = meshRep.Adjacent(Face(a, c, d), &b);
		if (!found || b->isGhost()) return false;
		std::swap(b, c);
	}
	if (c->isGhost()) {
		bool found = meshRep.Adjacent(Face(a, d, b), &c);
		if (!found || c->isGhost()) return false;
		std::swap(b, c);
	}
	if (d->isGhost()) {
		bool found = meshRep.Adjacent(Face(a, b, c), &d);
		if (!found || d->isGhost()) return false;
		std::swap(b, c);
	}

	const Predicator<REAL> predicator;
	enum Move { Left, Right, Horizon };
	bool terminate = false;
	do {
		bool ori0 = predicator.orient3d(bb, a->vert, d->vert, c->vert) >= 0;
		bool ori1 = predicator.orient3d(bb, a->vert, b->vert, d->vert) >= 0;
		bool ori2 = predicator.orient3d(bb, a->vert, c->vert, b->vert) >= 0;

		int condition = ori0 + (ori1 << 1) + (ori2 << 2);
		Move nextMove;
		switch (condition) {
		case 0:
		{
			switch (Randomnation(3)) {
			case 0:
				nextMove = Move::Left;
				break;
			case 1:
				nextMove = Move::Right;
				break;
			default:
				nextMove = Move::Horizon;
				break;
			}
			break;
		}
		case 1:
			if (Randomnation(2)) nextMove = Move::Right;
			else nextMove = Move::Horizon;
			break;
		case 2:
			if (Randomnation(2)) nextMove = Move::Left;
			else nextMove = Move::Horizon;
			break;
		case 3:
			nextMove = Move::Horizon;
			break;
		case 4:
			if (Randomnation(2)) nextMove = Move::Left;
			else nextMove = Move::Right;
			break;
		case 5:
			nextMove = Move::Right;
			break;
		case 6:
			nextMove = Move::Left;
			break;
		case 7:
			*t = Tetrahedron(a, b, c, d);
			return true;
		default:
			Severe("Unexpected Condition in DelMesher::findIntersectedTet");
			break;
		}

		switch (nextMove) {
		case Move::Left:
		{
			bool found = meshRep.Adjacent(Face(a, c, d), &b);
			terminate = !found || b->isGhost();
			break;
		}
		case Move::Right:
		{
			bool found = meshRep.Adjacent(Face(a, d, b), &c);
			terminate = !found || c->isGhost();
			break;
		}
		case Move::Horizon:
		{
			bool found = meshRep.Adjacent(Face(a, b, c), &d);
			terminate = !found || d->isGhost();
			break;
		}
		default:
			Severe("Unexpected Move Case in DelMesher::findIntersectedTet");
			break;
		}
		std::swap(b, c);

	} while (!terminate);

	return false;
}

Vertex* DelMesher::findSegmentEncroachedReference(Vertex *end, const Tetrahedron& intersected) const {
	Vertex *a = intersected.v[0], *b = intersected.v[1], *c = intersected.v[2], *d = intersected.v[3];

	DelVector startVert = a->vert, endVert = end->vert;
	Vertex *ref = b;
	REAL maxAngle = Geometer::interiorAngle(ref->vert, startVert, endVert);

	REAL angle = Geometer::interiorAngle(c->vert, startVert, endVert);
	if (angle > maxAngle) {
		ref = c;
		maxAngle = angle;
	}
	angle = Geometer::interiorAngle(d->vert, startVert, endVert);
	if (angle > maxAngle) {
		ref = d;
		maxAngle = angle;
	}

	enum Intersection { AcrossFace, AcrossEdge };
	Intersection intersect = Intersection::AcrossFace;
	if (predicator.orient3d(endVert, a->vert, c->vert, b->vert) == 0)
		intersect = Intersection::AcrossEdge;
	else if (predicator.orient3d(endVert, a->vert, d->vert, c->vert) == 0) {
		std::swap(d, b);
		std::swap(b, c);
		intersect = Intersection::AcrossEdge;
	}
	else if (predicator.orient3d(endVert, a->vert, b->vert, d->vert) == 0) {
		std::swap(b, c);
		std::swap(d, b);
		intersect = Intersection::AcrossEdge;
	}

	bool terminate = false;
	while (!terminate) {
		switch (intersect) {
		case Intersection::AcrossFace: //bcd is the intersected face
		{
			bool found = meshRep.Adjacent(Face(b, d, c), &a);
			Assert(found && !a->isGhost());
			if (a != end) {
				angle = Geometer::interiorAngle(a->vert, startVert, endVert);
				if (angle > maxAngle) {
					maxAngle = angle;
					ref = a;
				}
				std::swap(a, d);
				bool dbc = predicator.Intersection(d->vert, b->vert, c->vert, startVert, endVert);
				bool dab = predicator.Intersection(d->vert, a->vert, b->vert, startVert, endVert);
				bool dca = predicator.Intersection(d->vert, c->vert, a->vert, startVert, endVert);
				int condition = dbc + (dab << 1) + (dca << 2);
				switch (condition) {
				case 1:
					intersect = Intersection::AcrossFace;
					break;
				case 2:
					intersect = Intersection::AcrossFace;
					std::swap(c, a); std::swap(b, c);
					break;
				case 3:
					intersect = Intersection::AcrossEdge;
					std::swap(b, d); std::swap(c, d);
					break;
				case 4:
					intersect = Intersection::AcrossFace;
					std::swap(a, b); std::swap(b, c);
					break;
				case 5:
					intersect = Intersection::AcrossEdge;
					std::swap(c, d); std::swap(b, d);
					break;
				case 6:
					intersect = Intersection::AcrossEdge;
					std::swap(c, a); std::swap(b, d);
					break;
				default:
					Severe("Unexpected Case in DelMesher::findSegmentEncroachedReference");
					terminate = true;
					break;
				}
			}
			else
				terminate = true;

			break;
		}
		case Intersection::AcrossEdge: //bc is the intersected edge and bcd is the intersected face
		{
			bool hasIntersect = false;
			do {
				bool found = meshRep.Adjacent(Face(b, d, c), &a);
				Assert(found);
				if (a == end) {
					terminate = true;
					break;
				}
				if (!a->isGhost()) {
					angle = Geometer::interiorAngle(a->vert, startVert, endVert);
					if (angle > maxAngle) {
						maxAngle = angle;
						ref = a;
					}
					std::swap(a, d);
					bool dab = predicator.Intersection(d->vert, a->vert, b->vert, startVert, endVert);
					bool dca = predicator.Intersection(d->vert, c->vert, a->vert, startVert, endVert);
					bool coplane = predicator.orient3d(endVert, b->vert, c->vert, d->vert) == REAL(0);
					int condition = dab + (dca << 1);
					switch (condition) {
					case 0:
						break;
					case 1:
						if (!coplane) intersect = Intersection::AcrossFace;
						std::swap(a, c); std::swap(c, d);
						hasIntersect = true;
						break;
					case 2:
						if (!coplane) intersect = Intersection::AcrossFace;
						std::swap(a, b); std::swap(b, d);
						hasIntersect = true;
						break;
					case 3:
						std::swap(a, b); std::swap(c, d);
						hasIntersect = true;
						break;
					default:
						break;
					}
				}
				else {
					found = meshRep.Adjacent(Face(a, c, b), &d);
					Assert(found);
				}
			} while (!hasIntersect);
			break;
		}
		default:
			Severe("Unexpected Case in DelMesher::findSegmentEncroachedReference");
			terminate = true;
			break;
		}
	}

	return ref;
}

void DelMesher::splitSubSegment(const Segment& s, Vertex* ref, bool missingFaceTest){
	Vertex *a = s.v[0], *b = s.v[1];
	DelVector aa = a->vert, bb = b->vert;
	DelVector ab = bb - aa;
	REAL abLen = ab.length();
	REAL invabLen = REAL(1) / abLen;

	//calulate steiner point
	DelVector steinerPoint;
	Vertex *acuteVert = NULL, *otherVert = NULL;
	bool aIsAcute = matchVertexFlag(a->getVertexType(), VertexType::Vertex_Acute);
	bool bIsAcute = matchVertexFlag(b->getVertexType(), VertexType::Vertex_Acute);
	if (aIsAcute != bIsAcute) {
		if (aIsAcute) {
			acuteVert = a; otherVert = b; 
		}
		else if (bIsAcute) {
			acuteVert = b; otherVert = a;
		}
	}

	bool acuteSplit = acuteVert && meshRep.Contain(Segment(acuteVert, ref));
	if (acuteSplit) {
		DelVector acutePoint = acuteVert->vert;
		DelVector otherPoint = otherVert->vert;
		DelVector refPoint = ref->vert;
		REAL r = (refPoint - acutePoint).length();
		steinerPoint = acutePoint + (r * invabLen) * (otherPoint - acutePoint);

		REAL steOtherLen = (steinerPoint - otherPoint).length();
		REAL steRefLen = (steinerPoint - refPoint).length();

		if (steOtherLen < steRefLen) {
			r = nearestPowerOfTwo(REAL(0.5) * abLen);
			steinerPoint = acutePoint + (r * invabLen) * (otherPoint - acutePoint);
		}
	}
	else {
		DelVector ar = ref->vert - aa;
		DelVector br = ref->vert - bb;
		REAL arLen = ar.length();
		REAL brLen = br.length();
		if (arLen < REAL(0.5) * abLen)
			steinerPoint = aa + (arLen * invabLen) * ab;
		else if (brLen < REAL(0.5) * abLen)
			steinerPoint = bb - (brLen * invabLen) * ab;
		else
			steinerPoint = REAL(0.5) * (aa + bb);
	}

	//check to see steiner point isn't too close to both end point
	REAL t = (steinerPoint - aa).length() * invabLen;
	if (t < REAL(0.2) || t > REAL(0.8)) {
		if (!acuteSplit) 
			steinerPoint = REAL(0.5) * (aa + bb);
		else {
			REAL r = nearestPowerOfTwo(REAL(0.5) * abLen);
			steinerPoint = acuteVert->vert + (r * invabLen) * (otherVert->vert - acuteVert->vert);
		}
	}

	VertexType type = VertexType::Vertex_FreeSegment;
	bool onSurface = surfaceRep.Contain(s) || surfaceRep.Contain(Segment(s.v[1], s.v[0]));
	if (onSurface) type = VertexType(type | VertexType::Vertex_Facet);
	Vertex *steinerVert = allocVertex(steinerPoint, REAL(0), type);

	//insert two sub segments
	Segment s0(a, steinerVert);
	Segment s1(steinerVert, b);
	meshRep.deleteSegment(s.v[0], s.v[1]);
	meshRep.addSegment(s0.v[0], s0.v[1]);
	meshRep.addSegment(s1.v[0], s1.v[1]);
	mayMissingSegs.push_front(s0);
	mayMissingSegs.push_front(s1);

	//set original segment index
	if (matchVertexFlag(a->getVertexType(), VertexType::Vertex_FreeSegment))
		steinerVert->setOriSegmentIndex(a->getOriSegmentIndex());
	else if (matchVertexFlag(b->getVertexType(), VertexType::Vertex_FreeSegment))
		steinerVert->setOriSegmentIndex(b->getOriSegmentIndex());
	else {
		steinerVert->setOriSegmentIndex(oriSegments.size());
		oriSegments.push_back(s);
	}

	Tetrahedron toBeDeleted;
	meshRep.adjacent2Vertex(a, &toBeDeleted);
	toBeDeleted = findPosition(steinerVert, toBeDeleted, meshRep);
	REAL ori = predicator.orient3d(steinerVert->vert, toBeDeleted.v[0]->vert, toBeDeleted.v[1]->vert, toBeDeleted.v[2]->vert);

	VolumeVertexInsertionFlags vifs;
	vifs.missingSegTest = true;
	insertVertex(steinerVert, toBeDeleted, meshRep, vifs);

	//inserte it to surface
	if (onSurface) {
		SurfaceVertexInsertionFlags svifs;
		svifs.missingTest = missingFaceTest;
		insertSurfaceSegmentVertex(steinerVert, s, svifs);
	}
}

//all segments to be recovered must be in mayEncroachedSegs
void DelMesher::segmentsRecovery(bool missingTest) {
	while (!mayMissingSegs.empty()) {
		Segment s = mayMissingSegs.front();
		mayMissingSegs.pop_front();
		if (meshRep.isSegment(s) && !meshRep.Contain(s)) {
			Tetrahedron t;
			if (!meshRep.fastVertexQueryCheck(s.v[0])) std::swap(s.v[0], s.v[1]);
			if (findIntersectedTetrahedron(s.v[0], s.v[1]->vert, &t)) {
				Assert(s.v[1] != t.v[1] && s.v[1] != t.v[2] && s.v[1] != t.v[3]);
				splitSubSegment(s, findSegmentEncroachedReference(s.v[1], t), missingTest);
			}
			else
				Assert(false);
		}
	}
}

void DelMesher::splitSubSegment(const Segment &s){
	Vertex *a = s.v[0], *b = s.v[1];
	Vertex *mid = allocVertex((a->vert + b->vert) * REAL(0.5), REAL(0), VertexType::Vertex_FreeSegment);

	Vertex *c = NULL, *d = NULL;
	//should be orient inward
	Face f;
	SurfaceVertexInsertionFlags svifs;
	svifs.encroachmentTest = true;
	insertSurfaceSegmentVertex(mid, s, svifs, &f);

	Tetrahedron toBeDeleted;
	Vertex *e = NULL;
	if (matchVertexFlag(mid->getVertexType(), VertexType::Vertex_Facet) && meshRep.Adjacent(f, &e))
		toBeDeleted = Tetrahedron(e, f.v[0], f.v[1], f.v[2], true);
	else{
		meshRep.adjacent2Vertex(s.v[0], &toBeDeleted);
		toBeDeleted = findPosition(mid, toBeDeleted, meshRep);
	}
	VolumeVertexInsertionFlags vifs;
	vifs.encroachSegTest = true;
	insertVertex(mid, toBeDeleted, meshRep, vifs);
}

void DelMesher::splitSubPolygon(const Face &f){
	Vertex *a = f.v[0], *b = f.v[1], *c = f.v[2];
	Vertex center; 
	Geometer::Orthocircle(a->vert, a->weight, b->vert, b->weight, c->vert, c->weight, &center.vert);
	DelVector above = Geometer::calculateAbovePoint(a->vert, b->vert, c->vert);

	tobeDeletedFaces.push_back(f);
	surfaceRep.deleteTriangle(a, b, c);

	SurfaceVertexInsertionFlags svifs;
	svifs.trueInsertion = false;
	insertSurfaceVertex(&center, f, svifs);

	bool encorach = false;
	Segment e;
	for (auto ns : newSegsOfFaces){
		Segment seg = Segment(ns.v[0], ns.v[1], true);
		if (findSegment(seg)){
			if (Encroached(seg, &center)){
				encorach = true;
				e = seg;
				break;
			}
		}
	}

	if (!encorach){
		Vertex *newPoint = allocVertex(center.vert, center.weight, VertexType::Vertex_FreeFacet);
		for (auto segOfFace : newSegsOfFaces){
			surfaceRep.addTriangle(newPoint, segOfFace.v[0], segOfFace.v[1]);
			mayEncroachedFaces.push_back(Face(newPoint, segOfFace.v[0], segOfFace.v[1], true));
		}
		Tetrahedron toBeDeleted;
		Vertex *d;
		if (meshRep.Adjacent(f, &d))
			toBeDeleted = Tetrahedron(d, a, b, c, true);
		else{
			meshRep.adjacent2Vertex(a, &toBeDeleted);
			toBeDeleted = findPosition(newPoint, toBeDeleted, meshRep);
		}
		VolumeVertexInsertionFlags vifs;
		vifs.encroachSegTest = true; vifs.encroachFaceTest = true;
		insertVertex(newPoint, toBeDeleted, meshRep, vifs);
	}
	else{
		for (auto i : tobeDeletedFaces){
			surfaceRep.addTriangle(i.v[0], i.v[1], i.v[2]);
		}
		splitSubSegment(e);
	}
	tobeDeletedFaces.clear();
	newSegsOfFaces.clear();
}

void DelMesher::splitTetrahedron(const Tetrahedron& tet){
	Vertex *a = tet.v[0], *b = tet.v[1], *c = tet.v[2], *d = tet.v[3];
	Vertex center;
	Geometer::Orthosphere(a->vert, a->weight, b->vert, b->weight,
		c->vert, c->weight, d->vert, d->weight, &center.vert);

	tobeDeletedTets.push_back(tet);
	meshRep.deleteTetrahedron(a, b, c, d);

	VolumeVertexInsertionFlags vifs;
	vifs.trueInsertion = false;
	digCavity(&center, Face(b, d, c, true), meshRep, vifs, NULL);
	digCavity(&center, Face(a, c, d, true), meshRep, vifs, NULL);
	digCavity(&center, Face(a, d, b, true), meshRep, vifs, NULL);
	digCavity(&center, Face(a, b, c, true), meshRep, vifs, NULL);

	bool encrochment = false;
	Face encrocachedFace;
	bool result = false;
	for (auto t : newTets){
		Face nf(t.v[1], t.v[2], t.v[3]);
		result = surfaceRep.Contain(nf);
		if(!result){
			std::swap(nf.v[2], nf.v[1]);
			result = surfaceRep.Contain(nf);
			if (!result)
				continue;
		}
		if (Encroached(nf, &center)){
			encrocachedFace = nf;
			encrochment = true;
			break;
		}
	}

	if (!encrochment){
		Vertex *newVert = allocVertex(center.vert, center.weight, VertexType::Vertex_FreeVolume);
		for (auto t: newTets){
			meshRep.addTetrahedron(newVert, t.v[1], t.v[2], t.v[3]);
			Tetrahedron t(newVert, t.v[1], t.v[2], t.v[3], true);
			t.setRationAndRadius();
			if (t.getREration() > maxRatio || t.getRadius() > maxRadius)
			    skinnyTets.push(t);
		}
	}
	else{
		for (auto t : tobeDeletedTets)
			meshRep.addTetrahedron(t.v[0], t.v[1], t.v[2], t.v[3]);
		splitSubPolygon(encrocachedFace);
	}
	tobeDeletedTets.clear();
	newTets.clear();
}

//all faces to be recovered must be in mayMissingFaces
void DelMesher::facesRecovery() {
	//working vectors
	std::vector<Vertex *> regionVertices, positiveVertices, negativeVertices;
	std::vector<Segment> regionBoundaries;
	std::vector<Face> regionFaces, positiveFaces, negativeFaces;
	std::vector<Tetrahedron> deleted, inserted;

	while (!mayMissingFaces.empty()) {
		Face f = mayMissingFaces.front();
		mayMissingFaces.pop_front();
		if (surfaceRep.Contain(f) && !meshRep.Contain(f)) {
			if (!faceRecovery(f, regionVertices, regionBoundaries, regionFaces,
				positiveVertices, positiveFaces, negativeVertices, negativeFaces,
				deleted, inserted)) {
				refineRegion(f, false);
				segmentsRecovery(true);
			}
		}
	}
}

//also marking region vertices
void DelMesher::findMissingRegion(const Face& missed, std::vector<Vertex *> &regionVertices,
	std::vector<Segment>& regionBoundary, std::vector<Face>& regionFaces) {
	for (int i = 0; i < 3; i++) {
		regionVertices.push_back(missed.v[i]);
		missed.v[i]->setMark();
	}
	regionFaces.push_back(missed);
	surfaceRep.setDeletedMark(missed.v[0], missed.v[1], missed.v[2]);
	propagateFindRegion(Segment(missed.v[1], missed.v[0]),
		regionVertices, regionBoundary, regionFaces, 0);
	propagateFindRegion(Segment(missed.v[2], missed.v[1]),
		regionVertices, regionBoundary, regionFaces, 0);
	propagateFindRegion(Segment(missed.v[0], missed.v[2]),
		regionVertices, regionBoundary, regionFaces, 0);

	for (auto f : regionFaces) surfaceRep.unSetDeletedMark(f.v[0], f.v[1], f.v[2]);
}

void DelMesher::propagateFindRegion(const Segment& edge, std::vector<Vertex *> &regionVertices,
	std::vector<Segment>& regionBoundary, std::vector<Face>& regionFaces, int depth) {
	if (findSegment(Segment(edge.v[0], edge.v[1], true))) {
		regionBoundary.push_back(Segment(edge.v[1], edge.v[0]));
		return;
	}

	Vertex *w = NULL;
	int index = -1;
	if (!surfaceRep.Adjacent(edge, &w, &index)) return;

	Vertex *u = edge.v[0], *v = edge.v[1];
	Face f(u, v, w, index);
	if (!meshRep.Contain(f)) {
		if (!w->isMarked()) {
			regionVertices.push_back(w);
			w->setMark();
		}
		regionFaces.push_back(f);
		surfaceRep.setDeletedMark(u, v, w);
		propagateFindRegion(Segment(u, w), regionVertices, regionBoundary, regionFaces, depth + 1);
		propagateFindRegion(Segment(w, v), regionVertices, regionBoundary, regionFaces, depth + 1);
	}
	else
		regionBoundary.push_back(Segment(v, u));
}

bool DelMesher::findCrossEdge(const Segment& boundary, const std::vector<Face>& regionFaces, Segment& cross) const {
	Vertex *a = boundary.v[0], *b = boundary.v[1], *c = NULL, *d = NULL;
	Tetrahedron t;
	if (!findIntersectedTetrahedron(a, b->vert, &t)) return false;

	b = t.v[1]; c = t.v[2]; d = t.v[3];
	if (c == boundary.v[1]) {
		std::swap(b, c);
		std::swap(c, d);
	}
	else if (d == boundary.v[1]) {
		std::swap(b, d);
		std::swap(c, d);
	}

	if (b != boundary.v[1]) return false;

	bool hasIntersection = false;
	Vertex *anchor = c;
	do {
		if (!c->isGhost() && !d->isGhost()) {
			for (auto f : regionFaces) {
				if (predicator.Intersection(f.v[0]->vert, f.v[1]->vert, f.v[2]->vert, d->vert, c->vert) &&
					predicator.orient3d(c->vert, f.v[0]->vert, f.v[1]->vert, f.v[2]->vert) > REAL(0)) {
					hasIntersection = true;
					cross = Segment(c, d);
					break;
				}
			}
		}
		bool found = meshRep.Adjacent(Face(a, b, c), &d);
		Assert(found);
		std::swap(c, d);
	} while (c != anchor);

	return hasIntersection;
}

//require: all region vertices are marked
//ensure: all vertices in positiveVertices and negativeVertices are marked
bool DelMesher::findCavity(const std::vector<Segment>& regionBoundaries, std::vector<Face>& regionFaces,
	std::vector<Vertex *>& positiveVertices, std::vector<Face>& positiveFaces,
	std::vector<Vertex *>& negativeVertices, std::vector<Face>& negativeFaces,
	std::vector<Tetrahedron>& deleted) {
	enum CavityTriangleType { TwoPositive, TwoNegative, TwoCoplanar, PositiveNegativeCoplanar, NegativePositiveCoplanar };
	std::stack<std::pair<Face, CavityTriangleType>> triangleStack;
	std::vector<Segment> newSegments;
	std::vector<Face> toBeErased;
	toBeErased.reserve(regionFaces.size());

	for (auto s : regionBoundaries) surfaceRep.setMark(s.v[0], s.v[1]);

	bool positive = true;
	Segment seg = regionBoundaries[0];
	//get a cross edge
	Segment cross;
	bool found = findCrossEdge(seg, regionFaces, cross);
	Assert(found);
	//detect whether one of the vertex of the cross edge is on face
	if (!cross.v[1]->isMarked()) {
		//normal routine
		deleted.push_back(Tetrahedron(cross.v[1], seg.v[0], seg.v[1], cross.v[0]));
		meshRep.deleteTetrahedron(cross.v[1], seg.v[0], seg.v[1], cross.v[0]);

		positiveVertices.push_back(cross.v[0]);
		positiveFaces.push_back(Face(seg.v[1], seg.v[0], cross.v[0]));
		negativeVertices.push_back(cross.v[1]);
		negativeFaces.push_back(Face(seg.v[0], seg.v[1], cross.v[1]));

		triangleStack.push(std::make_pair(Face(cross.v[0], cross.v[1], seg.v[0]),
			CavityTriangleType::PositiveNegativeCoplanar));
		triangleStack.push(std::make_pair(Face(cross.v[1], cross.v[0], seg.v[1]),
			CavityTriangleType::NegativePositiveCoplanar));

		cross.v[0]->setMark(); cross.v[1]->setMark();
	}
	else {
		Vertex *posi = cross.v[0], *onFace = cross.v[1];
		toBeErased.push_back(Face(seg.v[0], seg.v[1], onFace));

		if (!surfaceRep.isMarked(onFace, seg.v[0])) {
			newSegments.push_back(Segment(onFace, seg.v[0]));
			triangleStack.push(std::make_pair(Face(onFace, seg.v[0], posi), CavityTriangleType::TwoCoplanar));
		}
		if (!surfaceRep.isMarked(seg.v[1], onFace)) {
			newSegments.push_back(Segment(seg.v[1], onFace));
			triangleStack.push(std::make_pair(Face(seg.v[1], onFace, posi), CavityTriangleType::TwoCoplanar));
		}
	}

	//find other cross tets
	while (!triangleStack.empty()) {
		auto pair = triangleStack.top();
		triangleStack.pop();

		Face intersect = pair.first;
		CavityTriangleType faceType = pair.second;

		Vertex *a = intersect.v[0], *b = intersect.v[1], *c = intersect.v[2], *d = NULL;
		if (!meshRep.Adjacent(intersect, &d)) continue;

		Assert(!d->isGhost());

		switch (faceType) {
		case CavityTriangleType::TwoPositive:
		{
			deleted.push_back(Tetrahedron(d, a, b, c));
			meshRep.deleteTetrahedron(d, a, b, c);
			Face intersected;
			bool hasIntersection = false;
			for (auto f : regionFaces) {
				if (predicator.Intersection(f.v[0]->vert, f.v[1]->vert, f.v[2]->vert, a->vert, b->vert, d->vert)) {
					intersected = f;
					hasIntersection = true;
					break;
				}
			}

			if (d != intersected.v[0] && d != intersected.v[1] && d != intersected.v[2]) {
				if (hasIntersection) {
					if (!d->isMarked()) {
						negativeVertices.push_back(d);
						d->setMark();
					}
					triangleStack.push(std::make_pair(Face(a, b, d), CavityTriangleType::TwoPositive));
					triangleStack.push(std::make_pair(Face(c, d, b), CavityTriangleType::TwoNegative));
					triangleStack.push(std::make_pair(Face(d, c, a), CavityTriangleType::TwoNegative));
				}
				else {
					if (!d->isMarked()) {
						positiveVertices.push_back(d);
						d->setMark();
					}
					positiveFaces.push_back(Face(b, a, d));
					triangleStack.push(std::make_pair(Face(a, d, c), CavityTriangleType::TwoPositive));
					triangleStack.push(std::make_pair(Face(d, b, c), CavityTriangleType::TwoPositive));
				}
			}
			else {
				Assert(hasIntersection);
				positiveFaces.push_back(Face(b, a, d));
				triangleStack.push(std::make_pair(Face(b, c, d), CavityTriangleType::PositiveNegativeCoplanar));
				triangleStack.push(std::make_pair(Face(c, a, d), CavityTriangleType::NegativePositiveCoplanar));
			}
			break;
		}
		case CavityTriangleType::TwoNegative:
		{
			deleted.push_back(Tetrahedron(d, a, b, c));
			meshRep.deleteTetrahedron(d, a, b, c);
			Face intersected;
			bool hasIntersection = false;
			for (auto f : regionFaces) {
				if (predicator.Intersection(f.v[0]->vert, f.v[1]->vert, f.v[2]->vert, a->vert, b->vert, d->vert)) {
					intersected = f;
					hasIntersection = true;
					break;
				}
			}

			if (d != intersected.v[0] && d != intersected.v[1] && d != intersected.v[2]) {
				if (hasIntersection) {
					if (!d->isMarked()) {
						positiveVertices.push_back(d);
						d->setMark();
					}
					triangleStack.push(std::make_pair(Face(a, b, d), CavityTriangleType::TwoNegative));
					triangleStack.push(std::make_pair(Face(c, d, b), CavityTriangleType::TwoPositive));
					triangleStack.push(std::make_pair(Face(d, c, a), CavityTriangleType::TwoPositive));
				}
				else {
					if (!d->isMarked()) {
						negativeVertices.push_back(d);
						d->setMark();
					}
					negativeFaces.push_back(Face(b, a, d));
					triangleStack.push(std::make_pair(Face(a, d, c), CavityTriangleType::TwoNegative));
					triangleStack.push(std::make_pair(Face(d, b, c), CavityTriangleType::TwoNegative));
				}
			}
			else {
				Assert(hasIntersection);
				negativeFaces.push_back(Face(b, a, d));
				triangleStack.push(std::make_pair(Face(b, c, d), CavityTriangleType::NegativePositiveCoplanar));
				triangleStack.push(std::make_pair(Face(c, a, d), CavityTriangleType::PositiveNegativeCoplanar));
			}
			break;
		}
		case CavityTriangleType::TwoCoplanar:
		{
			Face intersected;
			bool hasIntersection = false;
			do {
				for (auto f : regionFaces) {
					if (predicator.Intersection(f.v[0]->vert, f.v[1]->vert, f.v[2]->vert, c->vert, d->vert)) {
						intersected = f;
						hasIntersection = true;
						break;
					}
				}
				if (hasIntersection || !meshRep.Adjacent(Face(a, b, d), &c)) break;
				std::swap(c, d);
			} while (true);

			if (hasIntersection) {
				if (d != intersected.v[0] && d != intersected.v[1] && d != intersected.v[2]) {
					deleted.push_back(Tetrahedron(d, a, b, c));
					meshRep.deleteTetrahedron(d, a, b, c);
					if (!d->isMarked()) {
						negativeVertices.push_back(d);
						d->setMark();
					}
					if (!c->isMarked()) {
						positiveVertices.push_back(c);
						c->setMark();
					}
					negativeFaces.push_back(Face(d, b, a));
					positiveFaces.push_back(Face(c, a, b));
					triangleStack.push(std::make_pair(Face(d, c, a), CavityTriangleType::NegativePositiveCoplanar));
					triangleStack.push(std::make_pair(Face(c, d, b), CavityTriangleType::PositiveNegativeCoplanar));
				}
				else {
					toBeErased.push_back(Face(b, a, d));
					if (!surfaceRep.isMarked(a, d)) {
						newSegments.push_back(Segment(a, d));
						triangleStack.push(std::make_pair(Face(a, d, c), CavityTriangleType::TwoCoplanar));
					}
					if (!surfaceRep.isMarked(d, b)) {
						newSegments.push_back(Segment(d, b));
						triangleStack.push(std::make_pair(Face(d, b, c), CavityTriangleType::TwoCoplanar));
					}
				}
			}
			break;
		}
		case CavityTriangleType::PositiveNegativeCoplanar:
		{
			deleted.push_back(Tetrahedron(d, a, b, c));
			meshRep.deleteTetrahedron(d, a, b, c);
			DelVector pp = a->vert, qq = b->vert, rr = d->vert;
			Face intersected;
			bool hasIntersection = false;
			for (auto f : regionFaces) {
				if (predicator.Intersection(f.v[0]->vert, f.v[1]->vert, f.v[2]->vert, pp, rr) ||
					predicator.Intersection(f.v[0]->vert, f.v[1]->vert, f.v[2]->vert, qq, rr)) {
					intersected = f;
					hasIntersection = true;
					break;
				}
			}
			Assert(hasIntersection);

			REAL ori = predicator.orient3d(rr, intersected.v[0]->vert, intersected.v[1]->vert, intersected.v[2]->vert);
			if (ori > 0) {
				if (!d->isMarked()) {
					positiveVertices.push_back(d);
					d->setMark();
				}
				positiveFaces.push_back(Face(a, c, d));
				triangleStack.push(std::make_pair(Face(d, a, b), CavityTriangleType::TwoPositive));
				triangleStack.push(std::make_pair(Face(d, b, c), CavityTriangleType::PositiveNegativeCoplanar));
			}
			else if (ori < 0) {
				if (!d->isMarked()) {
					negativeVertices.push_back(d);
					d->setMark();
				}
				negativeFaces.push_back(Face(c, b, d));
				triangleStack.push(std::make_pair(Face(b, d, a), CavityTriangleType::TwoNegative));
				triangleStack.push(std::make_pair(Face(a, d, c), CavityTriangleType::PositiveNegativeCoplanar));
			}
			else {
				Assert(d->isMarked());
				//two vertices on face
				triangleStack.push(std::make_pair(Face(a, b, d), CavityTriangleType::PositiveNegativeCoplanar));
				if (!surfaceRep.isMarked(d, c)) {
					if (!surfaceRep.Contain(Segment(d, c))) newSegments.push_back(Segment(d, c));
					triangleStack.push(std::make_pair(Face(d, c, a), CavityTriangleType::TwoCoplanar));
				}
				positiveFaces.push_back(Face(a, c, d));
				negativeFaces.push_back(Face(c, b, d));
			}
			break;
		}
		case CavityTriangleType::NegativePositiveCoplanar:
		{
			deleted.push_back(Tetrahedron(d, a, b, c));
			meshRep.deleteTetrahedron(d, a, b, c);
			DelVector pp = a->vert, qq = b->vert, rr = d->vert;
			Face intersected;
			bool hasIntersection = false;
			for (auto f : regionFaces) {
				if (predicator.Intersection(f.v[0]->vert, f.v[1]->vert, f.v[2]->vert, pp, rr) ||
					predicator.Intersection(f.v[0]->vert, f.v[1]->vert, f.v[2]->vert, qq, rr)) {
					intersected = f;
					hasIntersection = true;
					break;
				}
			}
			Assert(hasIntersection);

			REAL ori = predicator.orient3d(rr, intersected.v[0]->vert, intersected.v[1]->vert, intersected.v[2]->vert);
			if (ori < 0) {
				if (!d->isMarked()) {
					negativeVertices.push_back(d);
					d->setMark();
				}
				negativeFaces.push_back(Face(a, c, d));
				triangleStack.push(std::make_pair(Face(d, a, b), CavityTriangleType::TwoNegative));
				triangleStack.push(std::make_pair(Face(d, b, c), CavityTriangleType::NegativePositiveCoplanar));
			}
			else if (ori > 0) {
				if (!d->isMarked()) {
					positiveVertices.push_back(d);
					d->setMark();
				}
				positiveFaces.push_back(Face(c, b, d));
				triangleStack.push(std::make_pair(Face(b, d, a), CavityTriangleType::TwoPositive));
				triangleStack.push(std::make_pair(Face(a, d, c), CavityTriangleType::NegativePositiveCoplanar));
			}
			else {
				Assert(d->isMarked());
				//two vertices on face
				triangleStack.push(std::make_pair(Face(a, b, d), CavityTriangleType::NegativePositiveCoplanar));
				if (!surfaceRep.isMarked(c, d)) {
					if (!surfaceRep.Contain(Segment(c, d))) newSegments.push_back(Segment(c, d));
					triangleStack.push(std::make_pair(Face(c, d, b), CavityTriangleType::TwoCoplanar));
				}
				negativeFaces.push_back(Face(a, c, d));
				positiveFaces.push_back(Face(c, b, d));
			}
			break;
		}
		default:
			Severe("Unexpected case in DelMesher::findCavity");
			break;
		}
	}

	for (auto s : regionBoundaries) surfaceRep.unSetMark(s.v[0], s.v[1]);

	if (newSegments.size() > 0) {
		int faceIndex = regionFaces[0].index;
		for (auto f : regionFaces) surfaceRep.deleteTriangle(f.v[0], f.v[1], f.v[2]);
		if (toBeErased.size() != regionFaces.size()) {
			DelTriangulator triangulator;
			for (auto f : toBeErased) surfaceRep.addTriangle(f.v[0], f.v[1], f.v[2], faceIndex);
			toBeErased.clear();
			triangulator.insertSegments(&regionFaces[0], &newSegments[0], regionFaces.size(), newSegments.size(), toBeErased);
			regionFaces.clear();
			for (auto f : toBeErased) {
				if (!surfaceRep.Contain(f)) {
					surfaceRep.addTriangle(f.v[0], f.v[1], f.v[2], faceIndex);
					mayMissingFaces.push_back(f);
					regionFaces.push_back(f);
				}
			}
		}
		else {
			for (auto f : toBeErased) 
				surfaceRep.addTriangle(f.v[0], f.v[1], f.v[2], faceIndex);

			return false;
		}
	}

	return true;
}

bool DelMesher::faceRecovery(Face& f, std::vector<Vertex *>& regionVertices,
	std::vector<Segment>& regionBoundaries, std::vector<Face>& regionFaces,
	std::vector<Vertex *>& positiveVertices, std::vector<Face>& positiveFaces,
	std::vector<Vertex *>& negativeVertices, std::vector<Face>& negativeFaces,
	std::vector<Tetrahedron>& deleted, std::vector<Tetrahedron>& inserted) {

	Face seed = f;
	regionFaces.clear(); regionVertices.clear(); regionBoundaries.clear();
	positiveVertices.clear(); positiveFaces.clear();
	negativeVertices.clear(); negativeFaces.clear();
	deleted.clear(); inserted.clear();

	//form region
	findMissingRegion(seed, regionVertices, regionBoundaries, regionFaces);

	bool success = true;
	//form cavity
	if (findCavity(regionBoundaries, regionFaces, 
		positiveVertices, positiveFaces, negativeVertices, negativeFaces,
		deleted)) {
		bool posSuccess = false, negSuccess = false;
		positiveVertices.reserve(positiveVertices.size() + regionVertices.size());
		for (auto v : regionVertices) positiveVertices.push_back(v);
		posSuccess = triangulateCavity(regionFaces, positiveFaces, positiveVertices, deleted, inserted, f);
		if (posSuccess) {
			for (auto& face : regionFaces) std::swap(face.v[1], face.v[2]);
			negativeVertices.reserve(negativeVertices.size() + regionVertices.size());
			for (auto v : regionVertices) negativeVertices.push_back(v);
			negSuccess = triangulateCavity(regionFaces, negativeFaces, negativeVertices, deleted, inserted, f);
		}

		success = posSuccess & negSuccess;

		if (success) {
			for (auto t : inserted) 
				meshRep.addTetrahedron(t.v[0], t.v[1], t.v[2], t.v[3]);
		}
		else {
			for (auto t : deleted)
				meshRep.addTetrahedron(t.v[0], t.v[1], t.v[2], t.v[3]);
		}
	}

	for (auto v : regionVertices) v->unSetMark();
	for (auto v : positiveVertices) v->unSetMark();
	for (auto v : negativeVertices) v->unSetMark();

	return success;
}

//require: all vertices in boundaryVertices are marked
bool DelMesher::triangulateCavity(const std::vector<Face>& regionFaces, std::vector<Face>& boundaryFaces,
	std::vector<Vertex *>& cavityVertices, std::vector<Tetrahedron>& deleted, std::vector<Tetrahedron>& inserted, Face& encroached) {
	size_t boundaryFacesSize = boundaryFaces.size();

	std::vector<Face> newBoundaryFaces;
	std::vector<Face> shellFaces;
	shellFaces.reserve(boundaryFacesSize);
	newBoundaryFaces.reserve(boundaryFacesSize);

	cavityRep.Clear();
	cavityRep.Reserve(boundaryFacesSize);

	//counstruct vertices map
	std::vector<Vertex *> newVertices;
	newVertices.reserve(boundaryFacesSize);
	for (auto oldVert : cavityVertices) {
		Vertex *newVert = cavityRep.allocVertex(oldVert->vert, oldVert->weight);
		newVertices.push_back(newVert);
		oldVert->setVertexPointer(newVert);
		newVert->setVertexPointer(oldVert);
	}

	//construct triangulation
	std::random_device rd;
	std::default_random_engine randomEngine(rd());
	std::shuffle(newVertices.begin(), newVertices.end(), randomEngine);
	triangulation3D(newVertices, cavityRep, false);

	int swapCount = 0;
	do {
		shellFaces.clear();
		bool enlarged = false;
		for (auto f : boundaryFaces) {
			Face face(f.v[0]->getPointedVertex(), f.v[1]->getPointedVertex(), f.v[2]->getPointedVertex());
			if (!cavityRep.Contain(face)) {
				Vertex *oppo = NULL;
				std::swap(f.v[1], f.v[2]);
				if (meshRep.Adjacent(f, &oppo)) {
					enlarged = true;
					meshRep.deleteTetrahedron(oppo, f.v[0], f.v[1], f.v[2]);
					deleted.push_back(Tetrahedron(oppo, f.v[0], f.v[1], f.v[2]));
					if (!oppo->isMarked()) {
						//alloc vertex
						Vertex *newOppo = cavityRep.allocVertex(oppo->vert, oppo->weight, oppo->getVertexType());
						oppo->setVertexPointer(newOppo);
						newOppo->setVertexPointer(oppo);

						//insert it
						Tetrahedron hint;
						bool found = cavityRep.adjacent2Vertex(face.v[0], &hint);
						Assert(found);
						hint = findPosition(newOppo, hint, cavityRep);
						insertVertex(newOppo, hint, cavityRep);

						oppo->setMark();
						cavityVertices.push_back(oppo);
					}
					shellFaces.push_back(Face(oppo, f.v[0], f.v[2]));
					shellFaces.push_back(Face(oppo, f.v[1], f.v[0]));
					shellFaces.push_back(Face(oppo, f.v[2], f.v[1]));
				}

				Face back(f.v[0], f.v[2], f.v[1]);
				if (surfaceRep.Contain(f)) mayMissingFaces.push_front(f);
				if (surfaceRep.Contain(back)) mayMissingFaces.push_front(back);
			}
			else
				shellFaces.push_back(face);
		}

		std::swap(shellFaces, boundaryFaces);
		swapCount += 1;
		if (!enlarged) break;
	} while (true);

	bool success = true;
	for (auto f : regionFaces) {
		Face face(f.v[0]->getPointedVertex(), f.v[1]->getPointedVertex(), f.v[2]->getPointedVertex());
		//check region faces in cavity and mark it
		if (cavityRep.Contain(face))
			cavityRep.setMark(face.v[0], face.v[1], face.v[2]);
		else {
			//if not, find missing triangle edge put to v[1]-v[2]
			if (!surfaceRep.Contain(f)) std::swap(f.v[1], f.v[2]);
			for (int i = 0; i < 3; i++) {
				Vertex *a = f.v[i], *b = f.v[NEXT_F(i)], *c = f.v[NEXT_F(NEXT_F(i))];
				if (!findSegment(Segment(b, c, true))) {
					encroached = Face(a, b, c);
					break;
				}
			}
			success = false;
			break;
		}
	}

	if (success) {
		//mark shell
		for (auto f : shellFaces) 
			cavityRep.setMark(f.v[0]->getPointedVertex(), f.v[1]->getPointedVertex(), f.v[2]->getPointedVertex());
		//clean tets out of cavity
		for (auto f : shellFaces) 
			propagateCleanCavity(Face(f.v[0]->getPointedVertex(), f.v[2]->getPointedVertex(), f.v[1]->getPointedVertex()), 0);
		//collect tets
		for (auto t : cavityRep)
			inserted.push_back(Tetrahedron(t.v[0]->getPointedVertex(), t.v[1]->getPointedVertex(), 
				t.v[2]->getPointedVertex(), t.v[3]->getPointedVertex()));
	}

	if (swapCount & 1) std::swap(shellFaces, boundaryFaces);

	for (auto v : cavityVertices) v->setVertexPointer(NULL);
	return success;
}

void DelMesher::propagateCleanCavity(const Face& f, int depth) {
	Vertex *a = NULL;
	if (!cavityRep.Adjacent(f, &a) || cavityRep.isMarked(f.v[0], f.v[1], f.v[2]))
		return;

	Vertex *b = f.v[0], *c = f.v[1], *d = f.v[2];
	cavityRep.deleteTetrahedron(a, b, c, d);

	propagateCleanCavity(Face(a, c, d), depth + 1);
	propagateCleanCavity(Face(a, d, b), depth + 1);
	propagateCleanCavity(Face(a, b, c), depth + 1);
}

void DelMesher::refineRegion(const Face& regionFace, bool encroachTest) {
	Vertex *steinerVert = allocVertex((regionFace.v[1]->vert + regionFace.v[2]->vert) * REAL(0.5), 
		REAL(0), VertexType::Vertex_FreeFacet);

	Tetrahedron toBeDeleted;
	meshRep.adjacent2Vertex(regionFace.v[1], &toBeDeleted);
	toBeDeleted = findPosition(regionFace.v[1], toBeDeleted, meshRep);

	Assert(mayMissingSegs.size() == 0);
	VolumeVertexInsertionFlags vifs;
	vifs.cdt = true; vifs.trueInsertion = false;
	vifs.missingSegTest = true; vifs.missingFaceTest = true; 
	vifs.encroachSegTest = encroachTest; vifs.encroachFaceTest = encroachTest;
	insertVertex(steinerVert, toBeDeleted, meshRep, vifs);
	//detect encroachment
	bool encroach = false;
	Segment encrochedSeg;
	while (!mayMissingSegs.empty()) {
		encrochedSeg = mayMissingSegs.back();
		if (Encroached(encrochedSeg, steinerVert)) {
			encroach = true;
			break;
		}
		mayMissingSegs.pop_back();
	}

	if (!encroach) {
		for (auto t : newTets)
			meshRep.addTetrahedron(t.v[0], t.v[1], t.v[2], t.v[3]);

		SurfaceVertexInsertionFlags svifs;
		svifs.missingTest = true; svifs.encroachmentTest = encroachTest;
		insertSurfaceVertex(steinerVert, regionFace, svifs);
	}
	else {
		//fall back
		for (auto t : tobeDeletedTets)
			meshRep.addTetrahedron(t.v[0], t.v[1], t.v[2], t.v[3]);

		mayMissingSegs.clear();

		//split encroch segment in mid point and update related information
		deallocVertex(steinerVert);
		bool onSurface = surfaceRep.Contain(encrochedSeg) || surfaceRep.Contain(Segment(encrochedSeg.v[1], encrochedSeg.v[0]));
		VertexType type = VertexType::Vertex_FreeSegment;
		if (onSurface) type = VertexType(type | VertexType::Vertex_Facet);
		steinerVert = allocVertex((encrochedSeg.v[0]->vert + encrochedSeg.v[1]->vert) * REAL(0.5), REAL(0), type);

		Segment s0(encrochedSeg.v[0], steinerVert), s1(encrochedSeg.v[1], steinerVert);
		meshRep.deleteSegment(encrochedSeg.v[0], encrochedSeg.v[1]);
		meshRep.addSegment(s0.v[0], s0.v[1]);
		meshRep.addSegment(s1.v[0], s1.v[1]);
		if (encroachTest) {
			mayEncroachedSegs.push_back(s0);
			mayEncroachedSegs.push_back(s1);
		}

		Vertex *a = encrochedSeg.v[0], *b = encrochedSeg.v[1];
		if (matchVertexFlag(a->getVertexType(), VertexType::Vertex_FreeSegment))
			steinerVert->setOriSegmentIndex(a->getOriSegmentIndex());
		else if (matchVertexFlag(b->getVertexType(), VertexType::Vertex_FreeSegment))
			steinerVert->setOriSegmentIndex(b->getOriSegmentIndex());
		else {
			steinerVert->setOriSegmentIndex(oriSegments.size());
			oriSegments.push_back(encrochedSeg);
		}

		bool found = meshRep.adjacent2SegmentFast(encrochedSeg, &toBeDeleted);
		Assert(found);

		VolumeVertexInsertionFlags newVifs;
		newVifs.cdt = true; newVifs.missingSegTest = true; newVifs.missingFaceTest = true;
		newVifs.encroachSegTest = encroachTest; newVifs.encroachFaceTest = encroachTest;
		insertVertex(steinerVert, toBeDeleted, meshRep, newVifs);

		//insert it to surface
		if (onSurface) {
			SurfaceVertexInsertionFlags newSvifs;
			newSvifs.missingTest = true; newSvifs.encroachmentTest = encroachTest;
			insertSurfaceSegmentVertex(steinerVert, encrochedSeg, newSvifs);
		}
	}
	newTets.clear();
	tobeDeletedTets.clear();
}

void DelMesher::insertVertex(Vertex *u, const Tetrahedron& tet, TetMeshDataStructure& meshRep,
	const VolumeVertexInsertionFlags& vifs, Tetrahedron *rt){
	Vertex *a = tet.v[0], *b = tet.v[1], *c = tet.v[2], *d = tet.v[3];

	if (!vifs.trueInsertion) tobeDeletedTets.push_back(tet);
	meshRep.deleteTetrahedron(a, b, c, d);

	if (vifs.encroachSegTest || vifs.missingSegTest) {
		markedSegments.clear();
		Segment segs[6];
		segs[0] = Segment(a, b); segs[1] = Segment(a, c); segs[2] = Segment(a, d); 
		segs[3] = Segment(b, c); segs[4] = Segment(b, d); segs[5] = Segment(c, d);
		for (int i = 0; i < 6; i++) {
			if (meshRep.isSegment(segs[i])) {
				meshRep.setMark(segs[i].v[0], segs[i].v[1]);
				markedSegments.push_back(segs[i]);
			}
		}
	}

	if (!vifs.cdt) {
		digCavity(u, Face(b, d, c, true), meshRep, vifs, rt);
		digCavity(u, Face(a, c, d, true), meshRep, vifs);
		digCavity(u, Face(a, d, b, true), meshRep, vifs);
		digCavity(u, Face(a, b, c, true), meshRep, vifs);
	}
	else {
		bool success = true;
		VolumeVertexInsertionFlags nvifs = vifs;
		nvifs.trueInsertion = false;
		success &= digCavity(u, Face(b, d, c, true), meshRep, nvifs);
		success &= digCavity(u, Face(a, c, d, true), meshRep, nvifs);
		success &= digCavity(u, Face(a, d, b, true), meshRep, nvifs);
		success &= digCavity(u, Face(a, b, c, true), meshRep, nvifs);

		if (!success) {
			std::vector<Vertex *> cavityVertices;
			cavityVertices.reserve(newTets.size());
			std::vector<Face> boundaryFaces;
			boundaryFaces.reserve(newTets.size());

			cavityVertices.push_back(u);
			u->setMark();
			for (auto t : newTets) {
				Face f = Face(t.v[1], t.v[2], t.v[3]);
				boundaryFaces.push_back(f);
				for (int i = 0; i < 3; i++) {
					if (!f.v[i]->isMarked()) {
						f.v[i]->setMark();
						cavityVertices.push_back(f.v[i]);
					}
				}
			}
			newTets.clear();

			int initDeletedCount = tobeDeletedTets.size();
			//triangulateCavity(std::vector<Face>(0), vifs.missingFaceTest, vifs.encroachFaceTest, 
				//boundaryFaces, cavityVertices, tobeDeletedTets, newTets, Face());
			for (auto v : cavityVertices) v->unSetMark();

			if (vifs.encroachSegTest || vifs.missingSegTest) {
				int finalDeletedCount = tobeDeletedTets.size();
				for (int i = initDeletedCount; i < finalDeletedCount; i++) {
					Tetrahedron t = tobeDeletedTets[i];
					Segment segs[6];
					segs[0] = Segment(t.v[0], t.v[1]); segs[1] = Segment(t.v[0], t.v[2]); segs[2] = Segment(t.v[0], t.v[3]);
					segs[3] = Segment(t.v[1], t.v[2]); segs[4] = Segment(t.v[1], t.v[3]); segs[5] = Segment(t.v[2], t.v[3]);
					for (int j = 0; j < 6; j++) {
						if (meshRep.isSegment(segs[j]) && !meshRep.testAndMark(segs[j].v[0], segs[j].v[1])) 
							markedSegments.push_back(segs[j]);
					}
				}
			}
		}

		if (vifs.trueInsertion) {
			for (auto t : newTets) meshRep.addTetrahedron(t.v[0], t.v[1], t.v[2], t.v[3]);
			newTets.clear();
			tobeDeletedTets.clear();
		}

		if (rt) meshRep.adjacent2Vertex(u, rt);
	}

	if (vifs.encroachSegTest || vifs.missingSegTest)
		for (auto s : markedSegments) meshRep.unSetMark(s.v[0], s.v[1]);

	if (vifs.encroachSegTest)
		for (auto s : markedSegments) mayEncroachedSegs.push_back(s);
	if (vifs.missingSegTest)
		for (auto s : markedSegments) mayMissingSegs.push_back(s);
}

void DelMesher::insertSurfaceVertex(Vertex *u, const Face& f, const SurfaceVertexInsertionFlags& vifs){
	if (vifs.trueInsertion) surfaceRep.bindVolumeVertex(u);
	Vertex *a = f.v[0], *b = f.v[1], *c = f.v[2];
	DelVector aboveVert = Geometer::calculateAbovePoint(a->vert, b->vert, c->vert);
	surfaceRep.deleteTriangle(a, b, c);
	digCavity(u, aboveVert, Segment(b, a), f.index, vifs);
	digCavity(u, aboveVert, Segment(c, b), f.index, vifs);
	digCavity(u, aboveVert, Segment(a, c), f.index, vifs);
}

void DelMesher::insertSurfaceSegmentVertex(Vertex *u, const Segment &s, const SurfaceVertexInsertionFlags& vifs, Face *inFace) {
	if (vifs.trueInsertion) surfaceRep.bindVolumeVertex(u);
	Vertex *a = s.v[0], *b = s.v[1];
	bool onFacet = false;
	Vertex *oppo = NULL;
	int index = -1;
	if (surfaceRep.Adjacent(s, &oppo, &index)) {
		onFacet = true;
		surfaceRep.deleteTriangle(a, b, oppo);
		DelVector above = Geometer::calculateAbovePoint(a->vert, b->vert, oppo->vert);
		digCavity(u, above, Segment(a, oppo), index, vifs);
		digCavity(u, above, Segment(oppo, b), index, vifs);

		if (inFace) *inFace = Face(b, a, oppo);
	}
	if (surfaceRep.Adjacent(Segment(b, a), &oppo, &index)) {
		onFacet = true;
		DelVector above = Geometer::calculateAbovePoint(b->vert, a->vert, oppo->vert);
		digCavity(u, above, Segment(b, oppo), index, vifs);
		digCavity(u, above, Segment(oppo, a), index, vifs);

		if (inFace) *inFace = Face(a, b, oppo);
	}
}

//Also mark all vertices
void DelMesher::detectCoplanarFaces(const Face& f, REAL facetRadianTol,
	std::vector<Vertex *>& coplanarVertices, std::vector<Segment>& boundaries, TriMeshDataStructure& surfRep) const {
	Vertex *a = f.v[0], *b = f.v[1], *c = f.v[2];
	for (int i = 0; i < 3; i++) {
		f.v[i]->setMark();
		coplanarVertices.push_back(f.v[i]);
	}
	surfRep.deleteTriangle(a, b, c);

	propagateDetectCoplanarFaces(c, Segment(b, a), facetRadianTol, coplanarVertices, boundaries, surfRep, 0);
	propagateDetectCoplanarFaces(a, Segment(c, b), facetRadianTol, coplanarVertices, boundaries, surfRep, 0);
	propagateDetectCoplanarFaces(b, Segment(a, c), facetRadianTol, coplanarVertices, boundaries, surfRep, 0);

	for (auto v : coplanarVertices) v->unSetMark();
}

void DelMesher::propagateDetectCoplanarFaces(Vertex *ref, const Segment& s, REAL facetRadianTol,
	std::vector<Vertex *>& coplanarVertices, std::vector<Segment>& boundaries, TriMeshDataStructure& surfRep, int depth) const {
	Vertex *u = s.v[0], *v = s.v[1], *w = NULL;
	if (!surfRep.Adjacent(s, &w)) {
		if (findSegment(Segment(u->getPointedVertex(), v->getPointedVertex()))) boundaries.push_back(Segment(v, u));
		return;
	}

	if (Geometer::dihedralAngle(u->vert, v->vert, ref->vert, w->vert) > facetRadianTol) {
		surfRep.deleteTriangle(u, v, w);
		if (!w->isMarked()) {
			w->setMark();
			coplanarVertices.push_back(w);
		}

		propagateDetectCoplanarFaces(u, Segment(w, v), facetRadianTol, coplanarVertices, boundaries, surfRep, depth + 1);
		propagateDetectCoplanarFaces(v, Segment(u, w), facetRadianTol, coplanarVertices, boundaries, surfRep, depth + 1);
	}
	else
		boundaries.push_back(Segment(v, u));
}

void DelMesher::propagateClean(const Face &f, int depth){
	if (depth < 100) {
		Vertex *a = NULL;
		if (surfaceRep.Contain(Face(f.v[0], f.v[2], f.v[1])) || !meshRep.Adjacent(f, &a))
			return;

		Vertex *b = f.v[0], *c = f.v[1], *d = f.v[2];

		meshRep.deleteTetrahedron(a, b, c, d);

		propagateClean(Face(a, c, d), depth + 1);
		propagateClean(Face(a, b, c), depth + 1);
		propagateClean(Face(a, d, b), depth + 1);
	}
	else {
		std::stack<Face> workingStacks;
		workingStacks.push(f);

		while (!workingStacks.empty()) {
			Face top = workingStacks.top();
			workingStacks.pop();

			Vertex *a = NULL;
			if (!surfaceRep.Contain(Face(top.v[0], top.v[2], top.v[1])) && meshRep.Adjacent(top, &a)) {
				Vertex *b = top.v[0], *c = top.v[1], *d = top.v[2];

				meshRep.deleteTetrahedron(a, b, c, d);

				workingStacks.push(Face(a, c, d));
				workingStacks.push(Face(a, b, c));
				workingStacks.push(Face(a, d, b));
			}
		}
	}
}

bool DelMesher::Encroached(const Segment& s) const {
	Vertex *a = s.v[0], *b = s.v[1], *c = NULL, *d = NULL;
	Tetrahedron t;
	if (!meshRep.adjacent2SegmentFast(s, &t)) return false;

	c = t.v[2]; d = t.v[3];

	DelVector center;
	REAL r = 0;
	Geometer::Orthocenter(a->vert, a->weight, b->vert, b->weight, &center, &r);
	REAL r2 = r * r;

	if (!d->isGhost() && (d->vert - center).length2() < r2) return true;

	Vertex *anchor = d;
	do {
		if (!c->isGhost() && (c->vert - center).length2() < r2)
			return true;

		bool found = meshRep.Adjacent(Face(a, b, c), &d);
		Assert(found);
		std::swap(c, d);
	} while (c != anchor);

	return false;
}

bool DelMesher::Encroached(const Face& f) const {
	Face f2(f.v[0], f.v[2], f.v[1], true);
	Vertex *x, *y;
	if (!meshRep.Adjacent(f, &x) || !meshRep.Adjacent(f2, &y))
		return true;

	DelVector center;
	REAL r = 0;
	Geometer::Orthocircle(f.v[0]->vert, f.v[0]->weight, f.v[1]->vert, f.v[1]->weight, f.v[2]->vert, f.v[2]->weight, &center, &r);
	REAL r2 = r*r;

	DelVector n;
	bool hasNormal = false;
	bool encroached0 = false, encroached1 = false;
	if (!x->isGhost() && ((x->vert - center).length2() - x->weight) <= r2) {
		//determine whether the encroached vertex is eligible or not
		if (Adjacent(f, x)) {
			n = Normalize(Geometer::triangleNormal(f.v[0]->vert, f.v[1]->vert, f.v[2]->vert));
			hasNormal = true;
			REAL dis = fabs((x->vert - center) * n);
			if (matchVertexFlag(x->getVertexType(), VertexType::Vertex_Segment))
				encroached0 = (dis > SQRTF_2 * x->relaxedInsetionRadius);
			else
				encroached0 = (dis > x->relaxedInsetionRadius);
		}
		else 
			encroached0 = true;
	}
	if (!encroached0 && !y->isGhost() && ((y->vert - center).length2() - y->weight) <= r2) {
		//determine whether the encroached vertex is eligible or not
		if (Adjacent(f, y)) {
			if(!hasNormal) n = Normalize(Geometer::triangleNormal(f.v[0]->vert, f.v[1]->vert, f.v[2]->vert));
			REAL dis = fabs((x->vert - center) * n);
			if (matchVertexFlag(x->getVertexType(), VertexType::Vertex_Segment))
				encroached1 = (dis > SQRTF_2 * y->relaxedInsetionRadius);
			else
				encroached1 = (dis > y->relaxedInsetionRadius);
		}
		else
			encroached1 = true;
	}

	return encroached0 || encroached1;
}

bool DelMesher::Encroached(const Segment &s, Vertex *v) const {
	DelVector center;
	REAL r = 0;
	Geometer::Orthocenter(s.v[0]->vert, s.v[0]->weight, s.v[1]->vert, s.v[1]->weight, &center, &r);

	if ((v->vert - center).length2() - v->weight <= r*r)
		return true;

	return false;
}

bool DelMesher::Encroached(const Face &f, Vertex *v) const {
	DelVector center;
	REAL r;
	Geometer::Orthocircle(f.v[0]->vert, f.v[0]->weight, f.v[1]->vert, f.v[1]->weight, f.v[2]->vert, f.v[2]->weight, &center, &r);
	REAL r2 = r * r;

	if (!v->isGhost() && ((v->vert - center).length2() - v->weight) <= r2) {
		//determine whether the encroached vertex is eligible or not
		if (Adjacent(f, v)) {
			DelVector n = Normalize(Geometer::triangleNormal(f.v[0]->vert, f.v[1]->vert, f.v[2]->vert));
			REAL dis = fabs((v->vert - center) * n);
			if (matchVertexFlag(v->getVertexType(), VertexType::Vertex_Segment))
				return dis > SQRTF_2 * v->relaxedInsetionRadius;
			else
				return dis > v->relaxedInsetionRadius;
		}
		else
			return true;
	}

	return false;
}

bool DelMesher::Adjacent(const Segment &s, Vertex *v) const {
	Segment oriSegment = s;
	if (matchVertexFlag(s.v[0]->getVertexType(), VertexType::Vertex_FreeSegment))
		oriSegment = oriSegments[s.v[0]->getOriSegmentIndex()];
	else if (matchVertexFlag(s.v[1]->getVertexType(), VertexType::Vertex_FreeSegment))
		oriSegment = oriSegments[s.v[1]->getOriSegmentIndex()];
	oriSegment.v[0]->setMark(); oriSegment.v[1]->setMark();

	bool adjacent = false;
	if (matchVertexFlag(v->getVertexType(), VertexType::Vertex_Segment)) {
		if (matchVertexFlag(v->getVertexType(), VertexType::Vertex_Free)) {
			Segment seg = oriSegments[v->getOriSegmentIndex()];
			adjacent = seg.v[0]->isMarked() || seg.v[1]->isMarked();
		}
		else
			adjacent = findSegment(Segment(v, oriSegment.v[0], true)) || findSegment(Segment(v, oriSegment.v[1], true));
	}
	else if (matchVertexFlag(v->getVertexType(), VertexType::Vertex_Facet)) {
		Face f;
		surfaceRep.adjacent2Vertex(v, &f);
		Vertex **verts = NULL;
		size_t size = getPolygonVertices(f.index, &verts);

		for (size_t i = 0; i < size; i++) {
			if (verts[i]->isMarked()) {
				adjacent = true;
				break;
			}
		}
	}
	oriSegment.v[0]->unSetMark(); oriSegment.v[1]->unSetMark();

	return adjacent;
}

bool DelMesher::Adjacent(const Face& f, Vertex *v) const {
	Vertex **verts = NULL;
	size_t size = getPolygonVertices(f.index, &verts);

	bool adjacent = false;
	if (matchVertexFlag(v->getVertexType(), VertexType::Vertex_Segment)) {
		if (matchVertexFlag(v->getVertexType(), VertexType::Vertex_Free)) {
			Segment oriSeg = oriSegments[v->getOriSegmentIndex()];
			oriSeg.v[0]->setMark(); oriSeg.v[1]->setMark();

			for (size_t i = 0; i < size; i++) {
				if (verts[i]->isMarked()) {
					adjacent = true;
					break;
				}
			}

			//clean up
			oriSeg.v[0]->unSetMark(); oriSeg.v[1]->unSetMark();
		}
		else {
			for (size_t i = 0; i < size; i++) {
				if (findSegment(Segment(v, verts[i], true))) {
					adjacent = true;
					break;
				}
			}
		}
	}
	else if (matchVertexFlag(v->getVertexType(), VertexType::Vertex_Facet)) {
		Face f2;
		surfaceRep.adjacent2Vertex(v, &f2);

		if (f.index == f2.index) return false;

		Vertex **verts2 = NULL;
		size_t size2 = getPolygonVertices(f2.index, &verts2);

		if (size2 < size) {
			std::swap(verts, verts2);
			std::swap(size, size2);
		}
		for (size_t i = 0; i < size; i++) verts[i]->setMark();
		for (size_t i = 0; i < size2; i++) {
			if (verts2[i]->isMarked()) {
				adjacent = true;
				break;
			}
		}
		//clean up
		for (size_t i = 0; i < size; i++) verts[i]->unSetMark();
	}

	return adjacent;
}

bool DelMesher::Adjacent(const Segment &s, const Face &f) const {
	Segment oriSegment = s;
	if (matchVertexFlag(s.v[0]->getVertexType(), VertexType::Vertex_FreeSegment)) 
		oriSegment = oriSegments[s.v[0]->getOriSegmentIndex()];
	else if (matchVertexFlag(s.v[1]->getVertexType(), VertexType::Vertex_FreeSegment)) 
		oriSegment = oriSegments[s.v[1]->getOriSegmentIndex()];
	oriSegment.v[0]->setMark(); oriSegment.v[1]->setMark();

	bool adjacent = false;
	Vertex **polyVerts;
	size_t size = getPolygonVertices(f.index, &polyVerts);

	for (size_t i = 0; i < size; i++) {
		if (polyVerts[i]->isMarked()) {
			adjacent = true;
			break;
		}
	}

	oriSegment.v[0]->unSetMark(); oriSegment.v[1]->unSetMark();
	return adjacent;
}


Reference<Mesh> DelMesher::generateMesh(int *vertexLableMap){
	for (auto t : meshRep) {
		//set vertices relaxed insertion radius
		Vertex *a = t.v[0], *b = t.v[1], *c = t.v[2], *d = t.v[3];
		DelVector ab = b->vert - a->vert, ac = c->vert - a->vert, ad = d->vert - a->vert,
			bc = c->vert - b->vert, bd = d->vert - b->vert, cd = d->vert - c->vert;
		REAL abLen = ab.length(), acLen = ac.length(), adLen = ad.length(),
			bcLen = bc.length(), bdLen = bd.length(), cdLen = cd.length();

		a->relaxedInsetionRadius = std::min({ a->relaxedInsetionRadius, abLen, acLen, adLen });
		b->relaxedInsetionRadius = std::min({ b->relaxedInsetionRadius, abLen, bcLen, bdLen });
		c->relaxedInsetionRadius = std::min({ c->relaxedInsetionRadius, acLen, bcLen, cdLen });
		d->relaxedInsetionRadius = std::min({ d->relaxedInsetionRadius, adLen, bdLen, cdLen });

		//push tets to skinny queue
		t.setRationAndRadius();
		if (t.getREration() > maxRatio || t.getRadius() > maxRadius)
			skinnyTets.push(t);
	}

	/************************
	       to be done
     ************************/

	//output to mesh
	auto tets = meshRep.getTetrahedrons(false);

	auto polygons = surfaceRep.getTriangles(false);

	int tetCount = tets.size();
	int vertexIndex = 0;
	constexpr int numTetNode = 4;
	std::map<Vertex*, int> vi;
	std::vector<int> tetNodeIndices(numTetNode * tetCount);

	for (int i = 0; i < tetCount; i++) {
		Tetrahedron t = tets[i];
		for (int j = 0; j < 4; j++) {
			if (!t.v[j]->isMarked()) {
				vi[t.v[j]] = vertexIndex++;
				t.v[j]->setMark();
			}
			tetNodeIndices[i * numTetNode + j] = vi[t.v[j]];
		}
	}

	int vertCount = vi.size();
	Mesh *mesh = new TetMesh(vertCount, tets.size(), polygons.size());

	for (int i = 0; i < tetCount; i++) 
		mesh->setElement(i, &tetNodeIndices[i * numTetNode]);

	int *vertIndices = allocAligned<int>(vertCount);
	MeshRelabeler labler(vertCount);
	labler.getNewLables(vertIndices, *mesh);

	int i = 0;
	for (auto pair : vi){
		Vertex *vertex = pair.first;
		mesh->setVertex(vertIndices[pair.second], Vector{ float(vertex->vert.x), float(vertex->vert.y), float(vertex->vert.z) });
	}

	i = 0;
	int surfVertIndices[3];
	for (auto f : polygons){
		surfVertIndices[0] = vertIndices[vi[f.v[0]]];
		surfVertIndices[1] = vertIndices[vi[f.v[1]]];
		surfVertIndices[2] = vertIndices[vi[f.v[2]]];
		mesh->setFacet(i++, surfVertIndices);
	}

	if (vertexLableMap) {
		for (auto oriv :oriVertices)
			vertexLableMap[oriv->getLabel() - 1] = vertIndices[vi[oriv]];
	}

	freeAligned(vertIndices);

	return mesh;
}
}
