#include "stdafx.h"
#include "delMesher.h"
#include "geometer.h"
#include <memory>
#include <set>
#include <random>
#include <stack>
#include <unordered_map>

namespace ODER{

Predicator<REAL> DelTriangulator::predicator;

void DelTriangulator::generateSubPolygons(Vertex **vertices, int *segments, int vertexCount, int segmentCount, 
	const Triangle& ref, bool boundaryOnly, std::vector<Triangle>& output){
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

	this->segments.reserve(segmentCount);
	for (int i = 0; i < segmentCount; i++)
		this->segments.push_back(Segment(vertices[segments[i * 2 + 0]]->getPointedVertex(), vertices[segments[i * 2 + 1]]->getPointedVertex()));

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

	Triangle f(v0, v1, v2);
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
		meshRep.addSegment(s.v[0], s.v[1]);
	}

	//clean out of region triangles
	if (boundaryOnly) {
		for (auto s : this->segments) {
			std::swap(s.v[0], s.v[1]);
			propagateClean(s, 0);
		}
		for (auto f : meshRep)
			output.push_back(Triangle(f.v[0]->getPointedVertex(), f.v[1]->getPointedVertex(), f.v[2]->getPointedVertex()));
	}
	else {
		for (auto f : meshRep) {
			if (!f.v[0]->isGhost())
				output.push_back(Triangle(f.v[0]->getPointedVertex(), f.v[1]->getPointedVertex(), f.v[2]->getPointedVertex()));
		}
	}
}

void DelTriangulator::generateSubPolygons(Vertex **vertices, Segment *segments, int vertexCount, int segmentCount, 
	const Triangle& ref, bool boundaryOnly, std::vector<Triangle>& output) {
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

	this->segments.reserve(segmentCount);
	for (int i = 0; i < segmentCount; i++) {
		Segment seg = segments[i];
		this->segments.push_back(Segment(seg.v[0]->getPointedVertex(), seg.v[1]->getPointedVertex()));
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

	Triangle f(v0, v1, v2);
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
		meshRep.addSegment(s.v[0], s.v[1]);
	}

	//clean out of region triangles
	if (boundaryOnly) {
		for (auto s : this->segments) {
			std::swap(s.v[0], s.v[1]);
			propagateClean(s, 0);
		}
		for (auto f : meshRep) 
			output.push_back(Triangle(f.v[0]->getPointedVertex(), f.v[1]->getPointedVertex(), f.v[2]->getPointedVertex()));
	}
	else {
		for (auto f : meshRep) {
			if (!f.v[0]->isGhost())
				output.push_back(Triangle(f.v[0]->getPointedVertex(), f.v[1]->getPointedVertex(), f.v[2]->getPointedVertex()));
		}
	}
}

void DelTriangulator::insertSegments(const Triangle *triangles, const Segment *segments, int triangleCount, int segmentCount, std::vector<Triangle>& output) {
	meshRep.Clear();
	meshRep.Reserve(triangleCount);
	this->segments.clear();

	//clean pointers
	std::vector<Vertex *> oldPointers;
	oldPointers.reserve(3 * triangleCount);
	for (int i = 0; i < triangleCount; i++) {
		Triangle f = triangles[i];
		for (int j = 0; j < 3; j++) {
			oldPointers.push_back(f.v[j]->getPointedVertex());
			f.v[j]->setVertexPointer(NULL);
		}
	}

	REAL maxArea2 = 0;
	Triangle largetsFace;
	for (int i = 0; i < triangleCount; i++) {
		Triangle f = triangles[i];
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
		this->segments.push_back(newSeg);
		if (!meshRep.Contain(newSeg))
			insertSegment(newSeg);
	}

	//restore pointers
	for (int i = 0; i < triangleCount; i++) {
		Triangle f = triangles[i];
		for (int j = 0; j < 3; j++) 
			f.v[j]->setVertexPointer(oldPointers[i * 3 + j]);
	}

	for (auto f : meshRep)
		output.push_back(Triangle(f.v[0]->getPointedVertex(), f.v[1]->getPointedVertex(), f.v[2]->getPointedVertex()));
}

void DelTriangulator::propagateClean(const Segment& s, int depth) {
	Vertex *w;
	if (!meshRep.Adjacent(s, &w)) return;

	Vertex *u = s.v[0], *v = s.v[1];

	meshRep.deleteTriangle(w, u, v);

	Segment uw(u, w), wv(w, v);
	if (!meshRep.isSegment(uw)) propagateClean(uw, depth + 1);
	if (!meshRep.isSegment(wv)) propagateClean(wv, depth + 1);
}

Triangle DelTriangulator::findPosition(Vertex *u, const Triangle &f) const{
	size_t maxIter = std::numeric_limits<size_t>::max();
	Vertex *a = f.v[0], *b = f.v[1], *c = f.v[2];
	//ghost triangle
	if (a->isGhost()){
		if (predicator.inHalfSpace2D(u->vert, b->vert, c->vert, ghost->vert))
			return f;
		Vertex *x = ghost;
		meshRep.Adjacent(Segment(c, b), &x);
		//Triangle(x, c, b)
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
					return Triangle(ghost, c, b);
				a = x; std::swap(b, c); //Triangle(x, c, b)
			}
			else
				return Triangle(a, b, c, true);
			break;
		case 1:
			meshRep.Adjacent(Segment(b, a), &x);
			if (x->isGhost())
				return Triangle(ghost, b, a);
			c = a; a = x; //Triangle(x, b, a)
			break;
		case 2:
			meshRep.Adjacent(Segment(a, c), &x);
			if (x->isGhost())
				return Triangle(ghost, a, c);
			b = a; a = x; //Triangle(x, a, c)
			break;
		case 3:
			if (Randomnation<2>()){
				meshRep.Adjacent(Segment(b, a), &x);
				if (x->isGhost())
					return Triangle(ghost, b, a);
				c = a; a = x; //Triangle(x, b, a)
			}
			else{
				meshRep.Adjacent(Segment(a, c), &x);
				if (x->isGhost())
					return Triangle(ghost, a, c);
				b = a; a = x; //Triangle(x, a, c)
			}
			break;
		default:
			Severe("Unexpected Case in DelTriangulator::findPosition");
		}
	}
	Severe("Do not find Triangle in DelTriangulator::findPosition");
	//should never return this value
	return Triangle(NULL, NULL, NULL);
}

void DelTriangulator::insertSegment(const Segment& s) {
	std::vector<Vertex *> positive, negtive;
	findCavity(s, positive, negtive);
	triangulateHalfHole(positive);
	triangulateHalfHole(negtive);
}

void DelTriangulator::findCavity(const Segment& s, std::vector<Vertex *>& positive, std::vector<Vertex *>& negtive) {
	Vertex *a = s.v[0], *b = s.v[1];
	if (Randomnation<2>()) std::swap(a, b);
	DelVector aa = a->vert, bb = b->vert, above = ghost->vert;
	Triangle f;

	positive.push_back(a);
	negtive.push_back(a);

	auto begin = meshRep.getIncidentTriangles(a);
	Assert(begin != NULL);
	auto circulator = begin;
	
	//find intersected facet
	bool found = false;
	do {
		f = *circulator;
		if (!f.v[1]->isGhost() && !f.v[2]->isGhost()) {
			if (predicator.orient2d(bb, f.v[1]->vert, f.v[2]->vert, above) < 0) {
				REAL ori1 = predicator.orient2d(aa, bb, f.v[1]->vert, above);
				REAL ori2 = predicator.orient2d(aa, bb, f.v[2]->vert, above);
				Assert(ori1 != 0 && ori2 != 0);
				if ((ori1 > 0) != (ori2 > 0)) found = true;
			}
		}
		++circulator;
	} while (!found && circulator != begin);
	Assert(found);

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

void DelTriangulator::digCavity(Vertex *u, const Segment &s, Triangle *rf, int depth){
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
			*rf = Triangle(u, verts[1], verts[0], true);
	}
}

DelVector DelTriangulator::calculateAbovePoint(int vertexCount, Vertex** vertices, const Triangle& ref) {
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

	DelVector n = REAL(0.5) * sqrt(longest) * Normalize(vec % (c->vert - a->vert));
	DelVector above = a->vert + n;
	if (predicator.orient3d(above, ref.v[0]->vert, ref.v[1]->vert, ref.v[2]->vert) < 0) {
		above = a->vert - n;
		std::swap(vertices[1], vertices[2]);
	}

	return above;
}

Predicator<REAL> DelMesher::predicator;

DelMesher::DelMesher(Vector *surfvs, int *triangles, int numv, int numtri, DelMesherCriteria cri) {
	criteria = cri;

	REAL maxREration = criteria.maxTetRatio;
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

	std::vector<Triangle> oriTriangles;
	oriTriangles.reserve(numtri);
	for (int i = 0; i < numtri; i++) {
		oriSurface.addTriangle(oriVertices[triangles[3 * i]], oriVertices[triangles[3 * i + 1]], oriVertices[triangles[3 * i + 2]]);
		oriTriangles.push_back(Triangle(oriVertices[triangles[3 * i]], oriVertices[triangles[3 * i + 1]], oriVertices[triangles[3 * i + 2]]));
	}

	REAL facetRadianTol = criteria.facetAngleTol / REAL(180) * REAL(M_PI);
	std::vector<Segment> boundaries;
	std::vector<Vertex *> coplanarVertices;
	DelTriangulator triangulator;
	std::vector<Triangle> faces;

	int faceIndex = 0;
	for (auto f : oriTriangles) {
		if (oriSurface.Contain(f)) {
			detectCoplanarFaces(f, facetRadianTol, coplanarVertices, boundaries, oriSurface);
			for (auto& s : boundaries) {
				for (int i = 0; i < 2; i++) {
					if (s.v[i]->getPointedVertex() == NULL) {
						Vertex *newVert = meshRep.allocVertex(s.v[i]->vert, s.v[i]->weight, 
							VertexType(Vertex_Fixed | Vertex_LowDimension));
						surfaceRep.bindVolumeVertex(newVert);
						s.v[i]->setVertexPointer(newVert);
					}
				}
				meshRep.addSegment(s.v[0]->getPointedVertex(), s.v[1]->getPointedVertex());
				std::swap(s.v[0], s.v[1]); //swap here for re-triangulation step followed
			}

			for (auto v : coplanarVertices) {
				if (v->getPointedVertex() == NULL) {
					Vertex *newVert = meshRep.allocVertex(v->vert, v->weight, Vertex_FixedFacet);
					surfaceRep.bindVolumeVertex(newVert);
					v->setVertexPointer(newVert);
				}
			}

			//re-triangulation
			//all boundary facets orient inward
			std::swap(f.v[0], f.v[1]);
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
	for (auto &v : oriVertices) {
		v = v->getPointedVertex();
		v->setVertexPointer(NULL);
	}

	detectAcuteVertices();

	//allocate space
	tobeDeletedFaces.reserve(16);
	newSegsOfFaces.reserve(32);
	tobeDeletedTets.reserve(16);
	newFacesOfTets.reserve(32);

	constrainedTriangulation();
}

DelMesher::~DelMesher() {
	for (auto p : verticesPerPolygon) delete[] p;
}

void DelMesher::triangulation3D(std::vector<Vertex *>& vertices, TetMeshDataStructure& meshRep){
	Assert(vertices.size() > 2);
	if (vertices.size() == 3) {
		Vertex *ghost = meshRep.getGhostVertex();
		meshRep.addTetrahedron(ghost, vertices[0], vertices[1], vertices[2]);
		meshRep.addTetrahedron(ghost, vertices[0], vertices[2], vertices[1]);
		return;
	}

	Vertex *a = vertices[0], *b = vertices[1];
	Vertex *c = NULL, *d = NULL;
	for (std::vector<Vertex *>::size_type i = 2; i < vertices.size(); i++) {
		if (predicator.orientCoplane(a->vert, b->vert, vertices[i]->vert) != REAL(0)) {
			c = vertices[i];
			std::swap(vertices[2], vertices[i]);
			break;
		}
	}
	Assert(c != NULL);
	for (std::vector<Vertex *>::size_type i = 3; i < vertices.size(); i++){
		REAL ori = predicator.orient3d(a->vert, b->vert, c->vert, vertices[i]->vert);
		if (ori != REAL(0)){
			d = vertices[i];
			std::swap(vertices[3], vertices[i]);
			if (ori < REAL(0)) std::swap(a, b);
			break;
		}
	}
	Assert(d != NULL);

	Vertex *ghost = meshRep.getGhostVertex();

	meshRep.addTetrahedron(a, b, c, d);
	meshRep.addTetrahedron(ghost, b, d, c);
	meshRep.addTetrahedron(ghost, a, c, d);
	meshRep.addTetrahedron(ghost, a, b, c);
	meshRep.addTetrahedron(ghost, a, d, b);

	REAL vvv = predicator.orient3d(a->vert, b->vert, c->vert, d->vert);
	VolumeVertexInsertionFlags vifs;
	Tetrahedron rt(a, b, c, d);
	for (std::vector<Vertex *>::size_type i = 4; i < vertices.size(); i++) {
		Vertex *v = vertices[i];
		Tetrahedron t = findPosition(v, rt, meshRep);
		insertVertex(v, t, meshRep, vifs, &rt);
	}
}

void DelMesher::constrainedTriangulation() {
	Geometer::brioSort3d(oriVertices.begin(), oriVertices.end(), [](Vertex *v)->DelVector { return v->vert; });
	triangulation3D(oriVertices, meshRep);
	for (auto iter = meshRep.segmentBegin(); iter != meshRep.segmentEnd(); ++iter) mayMissingSegs.push_back(*iter);
	std::sort(mayMissingSegs.begin(), mayMissingSegs.end(), 
		[](const Segment& left, const Segment& right) {
		return (left.v[0]->vert - left.v[1]->vert).length2() < (right.v[0]->vert - right.v[1]->vert).length2();
	});

	//recover segments
	while (!mayMissingSegs.empty()) {
		Segment s = mayMissingSegs.front();
		mayMissingSegs.pop_front();
		if (meshRep.isSegment(s) && !meshRep.Contain(s)) {
			Tetrahedron t;
			if (!meshRep.fastVertexQueryCheck(s.v[0])) std::swap(s.v[0], s.v[1]);
			if (findIntersectedTetrahedron(s.v[0], s.v[1]->vert, &t)) {
				Assert(s.v[1] != t.v[1] && s.v[1] != t.v[2] && s.v[1] != t.v[3]);
				splitSubSegment(s, findSegmentEncroachedReference(s.v[1], t), false);
			}
			else
				Assert(false);
		}
	}

	//collect faces to be reovered
	for (auto f : surfaceRep) mayMissingFacets.push_back(f);

	//working vectors for face recovery
	std::vector<Vertex *> regionVertices, positiveVertices, negativeVertices;
	std::vector<Segment> regionBoundaries;
	std::vector<Triangle> regionFaces, positiveFaces, negativeFaces;
	std::vector<Tetrahedron> deleted, inserted;
	//recover faces
	while (!mayMissingFacets.empty()) {
		Triangle f = mayMissingFacets.front();
		if (surfaceRep.Contain(f) && !meshRep.Contain(f)) {
			if (!facetRecovery(f, regionVertices, regionBoundaries, regionFaces,
				positiveVertices, positiveFaces, negativeVertices, negativeFaces,
				deleted, inserted)){
				refineRegion(f, true, false, false);
				//recover missing segs in the process of region refinment
				while (!mayMissingSegs.empty()) {
					Segment s = mayMissingSegs.front();
					mayMissingSegs.pop_front();
					if (meshRep.isSegment(s) && !meshRep.Contain(s)) {
						Tetrahedron t;
						if (!meshRep.fastVertexQueryCheck(s.v[0])) std::swap(s.v[0], s.v[1]);
						if (findIntersectedTetrahedron(s.v[0], s.v[1]->vert, &t)) {
							Assert(s.v[1] != t.v[1] && s.v[1] != t.v[2] && s.v[1] != t.v[3]);
							splitSubSegment(s, findSegmentEncroachedReference(s.v[1], t), true);
						}
						else
							Assert(false);
					}
				}
				continue;
			}
		}
		mayMissingFacets.pop_front();
	}
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
			meshRep.Adjacent(Triangle(b, d, c, true), &x);
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
					meshRep.Adjacent(Triangle(d, c, b), &x);
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
				meshRep.Adjacent(Triangle(a, b, c), &x);
				if (x->isGhost())
					return Tetrahedron(x, a, b, c, true);
				//Tetrahedron(x, a, b, c)
				d = c; c = b; b = a; a = x;
				break;
			case 2:
				meshRep.Adjacent(Triangle(a, c, d), &x);
				if (x->isGhost())
					return Tetrahedron(x, a, c, d, true);
				//Tetrahedron(x, a, c, d)
				b = a; a = x;
				break;
			case 3:
				if (Randomnation<2>() == 0){
					meshRep.Adjacent(Triangle(a, b, c), &x);
					if (x->isGhost())
						return Tetrahedron(x, a, b, c, true);
					//Tetrahedron(x, a, b, c)
					d = c; c = b; b = a; a = x;
				}
				else{
					meshRep.Adjacent(Triangle(a, c, d), &x);
					if (x->isGhost())
						return Tetrahedron(x, a, c, d, true);
					//Tetrahedron(x, a, c, d)
					b = a; a = x;
				}
				break;
			case 4:
				meshRep.Adjacent(Triangle(a, d, b), &x);
				if (x->isGhost())
					return Tetrahedron(x, a, d, b, true);
				//Tetrahedron(x, a, d, b)
				c = d; d = b; b = a; a = x;
				break;
			case 5:
				if (Randomnation<2>() == 0){
					meshRep.Adjacent(Triangle(a, b, c), &x);
					if (x->isGhost())
						return Tetrahedron(x, a, b, c, true);
					//Tetrahedron(x, a, b, c)
					d = c; c = b; b = a; a = x;
				}
				else{
					meshRep.Adjacent(Triangle(a, d, b), &x);
					if (x->isGhost())
						return Tetrahedron(x, a, d, b, true);
					//Tetrahedron(x, a, d, b)
					c = d; d = b; b = a; a = x;
				}
				break;
			case 6:
				if (Randomnation<2>() == 0){
					meshRep.Adjacent(Triangle(a, c, d), &x);
					if (x->isGhost())
						return Tetrahedron(x, a, c, d, true);
					//Tetrahedron(x, a, c, d)
					b = a; a = x;
				}
				else{
					meshRep.Adjacent(Triangle(a, d, b), &x);
					if (x->isGhost())
						return Tetrahedron(x, a, d, b, true);
					//Tetrahedron(x, a, d, b)
					c = d; d = b; b = a; a = x;
				}
				break;
			case 7:
			{
				unsigned int term = Randomnation<3>();
				if (term == 0){
					meshRep.Adjacent(Triangle(a, b, c), &x);
					if (x->isGhost())
						return Tetrahedron(x, a, b, c, true);
					//Tetrahedron(x, a, b, c)
					d = c; c = b; b = a; a = x;
				}
				else if (term == 1){
					meshRep.Adjacent(Triangle(a, c, d), &x);
					if (x->isGhost())
						return Tetrahedron(x, a, c, d, true);
					b = a; a = x;
				}
				else{
					meshRep.Adjacent(Triangle(a, d, b), &x);
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

Tetrahedron DelMesher::findPositionWithOcclusion(Vertex *u, const Tetrahedron& t, 
	const TetMeshDataStructure& meshRep, const TriMeshDataStructure& occluderRep) const {
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
		if (a->isGhost()) {
			if (predicator.inOrthoHalfSpace3D(u->vert, u->weight, b->vert, b->weight, c->vert, c->weight, d->vert, d->weight)) {
				if (occluderRep.Contain(Triangle(b, c, d))) return Tetrahedron(NULL, b, c, d);
				else if (occluderRep.Contain(Triangle(b, d, c))) return Tetrahedron(NULL, b, d, c);
				else return t;
			}
			meshRep.Adjacent(Triangle(b, d, c, true), &x);
			a = x; std::swap(d, c);
		}
		else {
			bool ori0 = predicator.orient3d(u->vert, a->vert, b->vert, c->vert) > 0;
			bool ori1 = predicator.orient3d(u->vert, a->vert, c->vert, d->vert) > 0;
			bool ori2 = predicator.orient3d(u->vert, a->vert, d->vert, b->vert) > 0;

			int condition = ori0 + (ori1 << 1) + (ori2 << 2);

			switch (condition) {
			case 0:
				if (predicator.orient3d(u->vert, d->vert, c->vert, b->vert) > 0) {
					if (occluderRep.Contain(Triangle(d, c, b))) return Tetrahedron(NULL, d, c, b);
					else if (occluderRep.Contain(Triangle(d, b, c))) return Tetrahedron(NULL, d, b, c);

					meshRep.Adjacent(Triangle(d, c, b), &x);
					if (x->isGhost())
						return Tetrahedron(x, d, c, b, true);
					//Tetrahedron(x, d, c, b)
					a = x; std::swap(b, d);
				}
				else {
					return Tetrahedron(a, b, c, d, true);
				}
				break;
			case 1:
				if (occluderRep.Contain(Triangle(a, b, c))) return Tetrahedron(NULL, a, b, c);
				else if (occluderRep.Contain(Triangle(a, c, b))) return Tetrahedron(NULL, a, c, b);

				meshRep.Adjacent(Triangle(a, b, c), &x);
				if (x->isGhost())
					return Tetrahedron(x, a, b, c, true);
				//Tetrahedron(x, a, b, c)
				d = c; c = b; b = a; a = x;
				break;
			case 2:
				if (occluderRep.Contain(Triangle(a, c, d))) return Tetrahedron(NULL, a, c, d);
				else if (occluderRep.Contain(Triangle(a, d, c))) return Tetrahedron(NULL, a, d, c);

				meshRep.Adjacent(Triangle(a, c, d), &x);
				if (x->isGhost())
					return Tetrahedron(x, a, c, d, true);
				//Tetrahedron(x, a, c, d)
				b = a; a = x;
				break;
			case 3:
				if (Randomnation<2>() == 0) {
					if (occluderRep.Contain(Triangle(a, b, c))) return Tetrahedron(NULL, a, b, c);
					else if (occluderRep.Contain(Triangle(a, c, b))) return Tetrahedron(NULL, a, c, b);

					meshRep.Adjacent(Triangle(a, b, c), &x);
					if (x->isGhost())
						return Tetrahedron(x, a, b, c, true);
					//Tetrahedron(x, a, b, c)
					d = c; c = b; b = a; a = x;
				}
				else {
					if (occluderRep.Contain(Triangle(a, c, d))) return Tetrahedron(NULL, a, c, d);
					else if (occluderRep.Contain(Triangle(a, d, c))) return Tetrahedron(NULL, a, d, c);

					meshRep.Adjacent(Triangle(a, c, d), &x);
					if (x->isGhost())
						return Tetrahedron(x, a, c, d, true);
					//Tetrahedron(x, a, c, d)
					b = a; a = x;
				}
				break;
			case 4:
				if (occluderRep.Contain(Triangle(a, d, b))) return Tetrahedron(NULL, a, d, b);
				else if (occluderRep.Contain(Triangle(a, b, d))) return Tetrahedron(NULL, a, b, d);

				meshRep.Adjacent(Triangle(a, d, b), &x);
				if (x->isGhost())
					return Tetrahedron(x, a, d, b, true);
				//Tetrahedron(x, a, d, b)
				c = d; d = b; b = a; a = x;
				break;
			case 5:
				if (Randomnation<2>() == 0) {
					if (occluderRep.Contain(Triangle(a, b, c))) return Tetrahedron(NULL, a, b, c);
					else if (occluderRep.Contain(Triangle(a, c, b))) return Tetrahedron(NULL, a, c, b);

					meshRep.Adjacent(Triangle(a, b, c), &x);
					if (x->isGhost())
						return Tetrahedron(x, a, b, c, true);
					//Tetrahedron(x, a, b, c)
					d = c; c = b; b = a; a = x;
				}
				else {
					if (occluderRep.Contain(Triangle(a, d, b))) return Tetrahedron(NULL, a, d, b);
					else if (occluderRep.Contain(Triangle(a, b, d))) return Tetrahedron(NULL, a, b, d);

					meshRep.Adjacent(Triangle(a, d, b), &x);
					if (x->isGhost())
						return Tetrahedron(x, a, d, b, true);
					//Tetrahedron(x, a, d, b)
					c = d; d = b; b = a; a = x;
				}
				break;
			case 6:
				if (Randomnation<2>() == 0) {
					if (occluderRep.Contain(Triangle(a, c, d))) return Tetrahedron(NULL, a, c, d);
					else if (occluderRep.Contain(Triangle(a, d, c))) return Tetrahedron(NULL, a, d, c);

					meshRep.Adjacent(Triangle(a, c, d), &x);
					if (x->isGhost())
						return Tetrahedron(x, a, c, d, true);
					//Tetrahedron(x, a, c, d)
					b = a; a = x;
				}
				else {
					if (occluderRep.Contain(Triangle(a, d, b))) return Tetrahedron(NULL, a, d, b);
					else if (occluderRep.Contain(Triangle(a, b, d))) return Tetrahedron(NULL, a, b, d);

					meshRep.Adjacent(Triangle(a, d, b), &x);
					if (x->isGhost())
						return Tetrahedron(x, a, d, b, true);
					//Tetrahedron(x, a, d, b)
					c = d; d = b; b = a; a = x;
				}
				break;
			case 7:
			{
				unsigned int term = Randomnation<3>();
				if (term == 0) {
					if (occluderRep.Contain(Triangle(a, b, c))) return Tetrahedron(NULL, a, b, c);
					else if (occluderRep.Contain(Triangle(a, c, b))) return Tetrahedron(NULL, a, c, b);

					meshRep.Adjacent(Triangle(a, b, c), &x);
					if (x->isGhost())
						return Tetrahedron(x, a, b, c, true);
					//Tetrahedron(x, a, b, c)
					d = c; c = b; b = a; a = x;
				}
				else if (term == 1) {
					if (occluderRep.Contain(Triangle(a, c, d))) return Tetrahedron(NULL, a, c, d);
					else if (occluderRep.Contain(Triangle(a, d, c))) return Tetrahedron(NULL, a, d, c);

					meshRep.Adjacent(Triangle(a, c, d), &x);
					if (x->isGhost())
						return Tetrahedron(x, a, c, d, true);
					b = a; a = x;
				}
				else {
					if (occluderRep.Contain(Triangle(a, d, b))) return Tetrahedron(NULL, a, d, b);
					else if (occluderRep.Contain(Triangle(a, b, d))) return Tetrahedron(NULL, a, b, d);

					meshRep.Adjacent(Triangle(a, d, b), &x);
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

Triangle DelMesher::findPosition(Vertex *u, const Triangle& f) const{
	size_t maxIter = std::numeric_limits<size_t>::max();
	Vertex *a = f.v[0], *b = f.v[1], *c = f.v[2];

	DelVector uProjected = Geometer::getProjectedPoint(a->vert, b->vert, c->vert, u->vert);
	auto plane = predicator.getProjectionPlane(a->vert, b->vert, c->vert);
	bool faceOrient = predicator.orientCoplane(a->vert, b->vert, c->vert, plane) > 0;

	auto orientTest = [plane, faceOrient](const DelVector& a, const DelVector& b, const DelVector& c) {
		constexpr Predicator<REAL> predicator;
		REAL ori = predicator.orientCoplane(a, b, c, plane);
		return faceOrient ? ori > 0 : ori < 0;
	};

	while (maxIter--){
		Vertex *x = NULL;
		bool ori0 = orientTest(b->vert, a->vert, uProjected);
		bool ori1 = orientTest(a->vert, c->vert, uProjected);

		int condition = ori0 + (ori1 << 1);
		switch (condition){
		case 0:
			if (orientTest(c->vert, b->vert, uProjected)){
				surfaceRep.Adjacent(Segment(c, b), &x);
				a = x; std::swap(b, c); //Triangle(x, c, b)
			}
			else {
				return Triangle(a, b, c, true);
			}
			break;
		case 1:
			surfaceRep.Adjacent(Segment(b, a), &x);
			c = a; a = x; //Triangle(x, b, a)
			break;
		case 2:
			surfaceRep.Adjacent(Segment(a, c), &x);
			b = a; a = x; //Triangle(x, a, c)
			break;
		case 3:
			if (Randomnation<2>() == 0){
				surfaceRep.Adjacent(Segment(b, a), &x);
				c = a; a = x; //Triangle(x, b, a)
			}
			else{
				surfaceRep.Adjacent(Segment(a, c), &x);
				b = a; a = x; //Triangle(x, a, c)
			}
			break;
		default:
			Severe("Unexpected Case in 2d overloaded DelMesher::findPosition");
		}
	}
	Severe("Do not find Triangle in 2d overloaded DelMesher::findPosition");
	//should never return this value
	return Triangle(NULL, NULL, NULL);
}

Triangle DelMesher::findPositionWithOcclusion(Vertex *u, const Triangle& f) const {
	size_t maxIter = std::numeric_limits<size_t>::max();
	Vertex *a = f.v[0], *b = f.v[1], *c = f.v[2];

	DelVector uProjected = Geometer::getProjectedPoint(a->vert, b->vert, c->vert, u->vert);
	auto plane = predicator.getProjectionPlane(a->vert, b->vert, c->vert);
	bool facetOrient = predicator.orientCoplane(a->vert, b->vert, c->vert, plane) > 0;
	
	auto orientTest = [plane, facetOrient](const DelVector& a, const DelVector& b, const DelVector& c) {
		constexpr Predicator<REAL> predicator;
		REAL ori = predicator.orientCoplane(a, b, c, plane);
		return facetOrient ? ori > 0 : ori < 0;
	};

	while (maxIter--) {
		Vertex *x = NULL;
		bool ori0 = orientTest(b->vert, a->vert, uProjected);
		bool ori1 = orientTest(a->vert, c->vert, uProjected);

		int condition = ori0 + (ori1 << 1);
		switch (condition) {
		case 0:
			if (orientTest(c->vert, b->vert, uProjected)) {
				if (meshRep.Contain(Segment(c, b))) return Triangle(NULL, c, b);
				surfaceRep.Adjacent(Segment(c, b), &x);
				a = x; std::swap(b, c); //Triangle(x, c, b)
			}
			else {
				return Triangle(a, b, c, true);
			}
			break;
		case 1:
			if (meshRep.Contain(Segment(b, a))) return Triangle(NULL, b, a);
			surfaceRep.Adjacent(Segment(b, a), &x);
			c = a; a = x; //Triangle(x, b, a)
			break;
		case 2:
			if (meshRep.Contain(Segment(a, c))) return Triangle(NULL, a, c);
			surfaceRep.Adjacent(Segment(a, c), &x);
			b = a; a = x; //Triangle(x, a, c)
			break;
		case 3:
			if (Randomnation<2>() == 0) {
				if (meshRep.Contain(Segment(b, a))) return Triangle(NULL, b, a);
				surfaceRep.Adjacent(Segment(b, a), &x);
				c = a; a = x; //Triangle(x, b, a)
			}
			else {
				if (meshRep.Contain(Segment(a, c))) return Triangle(NULL, a, c);
				surfaceRep.Adjacent(Segment(a, c), &x);
				b = a; a = x; //Triangle(x, a, c)
			}
			break;
		default:
			Severe("Unexpected Case in 2d overloaded DelMesher::findPosition");
		}
	}
	Severe("Do not find Triangle in 2d overloaded DelMesher::findPosition");
	//should never return this value
	return Triangle(NULL, NULL, NULL);
}

void DelMesher::digCavity(Vertex *u, const Triangle& f, TetMeshDataStructure& meshRep,
	const VolumeVertexInsertionFlags& vifs, Tetrahedron *rt){
	Vertex *a = NULL;
	if (!meshRep.Adjacent(f, &a)) {
		if (vifs.cdt)
			vifs.rejected |= (surfaceRep.Contain(f) || surfaceRep.Contain(Triangle(f.v[0], f.v[2], f.v[1])));
		return;
	}

	bool success = true;

	Vertex *b = f.v[0];
	Vertex *c = f.v[1];
	Vertex *d = f.v[2];

	bool deleted = false;
	bool constrained = false;
	if (vifs.cdt) constrained = surfaceRep.Contain(f) || surfaceRep.Contain(Triangle(b, d, c));
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

		if (vifs.missingFacetTest || vifs.encroachFacetTest) {
			Triangle front(b, c, d), back(b, d, c);
			if (surfaceRep.Contain(front)) {
				if (vifs.missingFacetTest) mayMissingFacets.push_back(front);
				if (vifs.encroachFacetTest) mayEncroachedFacets.push_back(front);
			}
			else if (surfaceRep.Contain(back)) {
				if (vifs.missingFacetTest) mayMissingFacets.push_back(back);
				if (vifs.encroachFacetTest) mayEncroachedFacets.push_back(back);
			}
		}
		if (vifs.encroachSegTest || vifs.missingSegTest) {
			Segment ab(a, b), ac(a, c), ad(a, d);
			if (meshRep.isSegment(ab) && !meshRep.testAndMark(a, b)) markedSegments.push_back(ab);
			if (meshRep.isSegment(ac) && !meshRep.testAndMark(a, c)) markedSegments.push_back(ac);
			if (meshRep.isSegment(ad) && !meshRep.testAndMark(a, d)) markedSegments.push_back(ad);
		}

		digCavity(u, Triangle(a, c, d, true), meshRep, vifs, rt);
		digCavity(u, Triangle(a, b, c, true), meshRep, vifs, rt);
		digCavity(u, Triangle(a, d, b, true), meshRep, vifs, rt);
	}
	else {
		if (vifs.trueInsertion && !vifs.cdt) meshRep.addTetrahedron(u, b, d, c);
		else newFacesOfTets.push_back(Triangle(b, d, c));

		if (vifs.cdt && !b->isGhost()) 
			vifs.rejected |= (predicator.orient3d(u->vert, b->vert, d->vert, c->vert) <= REAL(0));

		if (rt) *rt = Tetrahedron(u, b, d, c, true);
		if (vifs.encroachFacetTest) {
			Triangle front(b, d, c), back(b, c, d);
			if (surfaceRep.Contain(front)) mayEncroachedFacets.push_back(front);
			if (surfaceRep.Contain(back)) mayEncroachedFacets.push_back(back);
		}
		if (vifs.skinnyTetTest){
			Tetrahedron t(u, b, d, c);
			if (skinnyTetTest(t)) skinnyTets.push(t);
		}
		if (vifs.insertRadiusTest) {
			REAL minLen2 = (vifs.parent != NULL) ? (vifs.parent->vert - u->vert).length2() 
				: std::numeric_limits<REAL>::max();
			for (int i = 0; i < 3; i++) {
				if (!f.v[i]->isGhost()) {
					REAL len2 = (f.v[i]->vert - u->vert).length2();
					if (minLen2 > len2) {
						vifs.parent = f.v[i];
						minLen2 = len2;
					}
				}
			}
		}
	}
}

void DelMesher::digCavity(Vertex *u, const DelVector& above,const Segment &s, int index, const SurfaceVertexInsertionFlags& vifs){
	Vertex *w = s.v[0], *v = s.v[1], *x = NULL;
	int newIndex = -1;
	if (!surfaceRep.Adjacent(s, &x, &newIndex)) return;

	if (!findSegment(s) &&
		predicator.inOrthoCirclePerturbed(u->vert, u->weight, v->vert, v->weight, w->vert, w->weight, x->vert, x->weight, above) > 0){
		if (vifs.deletionRecord)
			tobeDeletedFaces.push_back(TriangleWithIndex(w, v, x, index, true));
		surfaceRep.deleteTriangle(w, v, x);
		Assert(index == newIndex);
		digCavity(u, above, Segment(x, v), newIndex, vifs);
		digCavity(u, above, Segment(w, x), newIndex, vifs);
	}
	else {
		if (vifs.trueInsertion)
			surfaceRep.addTriangle(u, v, w, index);
		else 
			newSegsOfFaces.push_back(SegmentWithIndex(v, w, index));
		if (vifs.missingTest) mayMissingFacets.push_back(Triangle(u, v, w));
		if (vifs.encroachmentTest) mayEncroachedFacets.push_back(Triangle(u, v, w));
	}
}

bool DelMesher::findIntersectedTetrahedron(Vertex *a, const DelVector& bb, Tetrahedron *t) const {
	if (!meshRep.adjacent2Vertex(a, t)) return false;

	Vertex *b = t->v[1], *c = t->v[2], *d = t->v[3];
	if (b->isGhost()) {
		bool found = meshRep.Adjacent(Triangle(a, c, d), &b);
		if (!found || b->isGhost()) return false;
		std::swap(b, c);
	}
	if (c->isGhost()) {
		bool found = meshRep.Adjacent(Triangle(a, d, b), &c);
		if (!found || c->isGhost()) return false;
		std::swap(b, c);
	}
	if (d->isGhost()) {
		bool found = meshRep.Adjacent(Triangle(a, b, c), &d);
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
			switch (Randomnation<3>()) {
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
			if (Randomnation<2>()) nextMove = Move::Right;
			else nextMove = Move::Horizon;
			break;
		case 2:
			if (Randomnation<2>()) nextMove = Move::Left;
			else nextMove = Move::Horizon;
			break;
		case 3:
			nextMove = Move::Horizon;
			break;
		case 4:
			if (Randomnation<2>()) nextMove = Move::Left;
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
			bool found = meshRep.Adjacent(Triangle(a, c, d), &b);
			terminate = !found || b->isGhost();
			break;
		}
		case Move::Right:
		{
			bool found = meshRep.Adjacent(Triangle(a, d, b), &c);
			terminate = !found || c->isGhost();
			break;
		}
		case Move::Horizon:
		{
			bool found = meshRep.Adjacent(Triangle(a, b, c), &d);
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
			bool found = meshRep.Adjacent(Triangle(b, d, c), &a);
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
				bool found = meshRep.Adjacent(Triangle(b, d, c), &a);
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
					found = meshRep.Adjacent(Triangle(a, c, b), &d);
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

void DelMesher::splitSubSegment(const Segment& s, Vertex* ref, bool missingFacetTest){
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

	VolumeVertexInsertionFlags vifs;
	vifs.missingSegTest = true; vifs.missingFacetTest = missingFacetTest;
	insertVertexOnSegment(steinerVert, s, toBeDeleted, meshRep, vifs);
}

void DelMesher::constrainedRefinment() {
	//working vectors for facet recovery
	std::vector<Vertex *> regionVertices, positiveVertices, negativeVertices;
	std::vector<Segment> regionBoundaries;
	std::vector<Triangle> regionFaces, positiveFaces, negativeFaces;
	std::vector<Tetrahedron> deleted, inserted;

	bool encorachedFacetTest = false, skinnyTest = false;
	bool facetsCollected = false, tetsCollected = false;

	//collect all subsegments
	for (auto iter = meshRep.segmentBegin(); iter != meshRep.segmentEnd(); ++iter) mayEncroachedSegs.push_back(*iter);
refinmentStart:
	//phase 1: refine subsegments that do not meet criteria
	while (!mayEncroachedSegs.empty()) {
		Segment s = mayEncroachedSegs.front();
		mayEncroachedSegs.pop_front();
		REAL length = (s.v[0]->vert - s.v[1]->vert).length();
		if (meshRep.isSegment(s) && (Encroached(s) || length > criteria.maxSegLength)) 
			refineSubSegment(s, NULL, false, encorachedFacetTest, skinnyTest);
	}
	//recover missing subpolygons
	while (!mayMissingFacets.empty()) {
		Triangle f = mayMissingFacets.front();
		if (surfaceRep.Contain(f) && !meshRep.Contain(f)) {
			if (!facetRecovery(f, regionVertices, regionBoundaries, regionFaces,
				positiveVertices, positiveFaces, negativeVertices, negativeFaces,
				deleted, inserted)) {
				refineRegion(f, false, true, true);
				goto refinmentStart;
			}
		}
		mayMissingFacets.pop_front();
	}

	if (!facetsCollected) { //collect subpolygons for the first round
		for (auto f : surfaceRep) mayEncroachedFacets.push_back(f);
		facetsCollected = true;
	}
	encorachedFacetTest = true;

	//phase 2: refine subpolygons that do not meet criteria
	while (!mayEncroachedFacets.empty()) {
		Triangle f = mayEncroachedFacets.front();
		if (surfaceRep.Contain(f)) {
			Vertex *encroached = NULL;
			REAL radius = REAL(0);
			Geometer::Orthocircle(f.v[0]->vert, f.v[0]->weight, f.v[1]->vert, f.v[1]->weight, f.v[2]->vert, f.v[2]->weight,
				(DelVector *)NULL, &radius);

			if (Encroached(f, &encroached)) {
				Assert(encroached);
				f = findPositionWithOcclusion(encroached, f);
				if (f.v[0]) refineSubPolygon(f, skinnyTest);
				else {
					Segment s(f.v[1], f.v[2]);
					bool entwinement = false;
					if (matchVertexFlag(encroached->getVertexType(), VertexType(VertexType::Vertex_Refined 
						| VertexType::Vertex_Facet | VertexType::Vertex_Volume)))
						entwinement = Adjacent(s, encroached);
					else encroached = NULL;
					refineSubSegment(s, encroached, entwinement, encorachedFacetTest, skinnyTest);
				}
			}
			else if (radius > criteria.maxFacetRadius) {
				refineSubPolygon(f, skinnyTest);
				mayEncroachedFacets.pop_front();
			}

			if (!mayEncroachedSegs.empty()) goto refinmentStart;
			while (!mayMissingFacets.empty()) {
				Triangle f = mayMissingFacets.front();
				if (!meshRep.Contain(f)) goto refinmentStart;
				mayMissingFacets.pop_front();
			}
		}
		else mayEncroachedFacets.pop_front();
	}

	if (!tetsCollected) { // collect skinny tets for the first round
		for (auto t : meshRep) {
			if (skinnyTetTest(t)) skinnyTets.push(t);
		}
		tetsCollected = true;
	}
	skinnyTest = true;
	
	// phase 3: refine tetrahedrons that do not meet criteria
	while (!skinnyTets.empty()) {
		Tetrahedron t = skinnyTets.top();
		if (meshRep.Contain(t)) {
			bool processed = false;
			refineTetrahedron(t, processed);
			if (processed) skinnyTets.pop();
			if (!mayEncroachedSegs.empty() || !mayEncroachedFacets.empty() || !mayMissingFacets.empty()) 
				goto refinmentStart;
		}
		else skinnyTets.pop();
	}
}

void DelMesher::refineSubSegment(const Segment &s, Vertex *parent, bool entwinement, bool encorachedFaceTestEnable, bool skinnyTetTestEnable){
	Vertex *a = s.v[0], *b = s.v[1];
	DelVector steinerPoint, steinerPointOther;
	bool aIsAcute = matchVertexFlag(a->getVertexType(), VertexType::Vertex_Acute);
	bool bIsAcute = matchVertexFlag(b->getVertexType(), VertexType::Vertex_Acute);

	bool bothAcute = false;
	bool midSplit = true;
	if (aIsAcute || bIsAcute) {
		midSplit = false;
		bothAcute = (aIsAcute == bIsAcute);
		if (!bothAcute)
			midSplit = aIsAcute ? matchVertexFlag(b->getVertexType(), VertexType::Vertex_Refined) 
			: matchVertexFlag(a->getVertexType(), VertexType::Vertex_Refined);
	}

	if (midSplit)
		steinerPoint = (a->vert + b->vert) * REAL(0.5);
	else {
		//calulate steiner point for acute vertex
		DelVector ab = b->vert - a->vert;
		REAL abLen = ab.length();
		REAL r = nearestPowerOfTwo(REAL(0.5) * abLen);
		if (!bothAcute) {
			if (aIsAcute) steinerPoint = a->vert + (r / abLen) * ab;
			else steinerPoint = b->vert - (r / abLen) * ab;
		}
		else {
			if (r > abLen * REAL(0.5)) r *= REAL(0.5);
			steinerPoint = a->vert + (r / abLen) * ab;

			DelVector bsp = steinerPoint - b->vert;
			REAL bspLen = bsp.length();
			REAL r2 = nearestPowerOfTwo(REAL(0.5) * bspLen);
			steinerPointOther = b->vert + (r2 / bspLen) * bsp;
		}
	}

	VertexType type = VertexType::Vertex_RefinedFreeSegment;
	if (surfaceRep.Contain(s) || surfaceRep.Contain(Segment(b, a))) 
		type = VertexType(type | VertexType::Vertex_Facet);
	Vertex *steinerVert = allocVertex(steinerPoint, REAL(0), type);

	Tetrahedron toBeDeleted;
	if (!meshRep.adjacent2SegmentFast(s, &toBeDeleted)) {
		meshRep.adjacent2Vertex(a, &toBeDeleted);
		toBeDeleted = findPosition(steinerVert, toBeDeleted, meshRep);
	}

	meshRep.deleteSegment(a, b);
	meshRep.addSegment(a, steinerVert);
	meshRep.addSegment(b, steinerVert);
	mayEncroachedSegs.push_back(Segment(a, steinerVert));
	if (!bothAcute) mayEncroachedSegs.push_back(Segment(b, steinerVert));

	VolumeVertexInsertionFlags vifs;
	vifs.cdt = true;
	vifs.missingFacetTest = true;
	vifs.encroachSegTest = true; vifs.encroachFacetTest = encorachedFaceTestEnable;
	vifs.skinnyTetTest = skinnyTetTestEnable;
	vifs.insertRadiusTest = (parent == NULL) || (entwinement && midSplit);
	insertVertexOnSegment(steinerVert, s, toBeDeleted, meshRep, vifs);

	if (vifs.insertRadiusTest) {
		steinerVert->relaxedInsetionRadius = (vifs.parent->vert - steinerPoint).length();
		if (midSplit) {
			if (Adjacent(s, vifs.parent) || vifs.parent == a || vifs.parent == b)
				steinerVert->relaxedInsetionRadius = std::max(steinerVert->relaxedInsetionRadius,
					vifs.parent->relaxedInsetionRadius);
		}
		else {
			if (Adjacent(s, vifs.parent))
				steinerVert->relaxedInsetionRadius = std::max(steinerVert->relaxedInsetionRadius,
					vifs.parent->relaxedInsetionRadius);
		}

		if (entwinement) {
			steinerVert->relaxedInsetionRadius = std::min(parent->relaxedInsetionRadius / SQRTF_2,
				steinerVert->relaxedInsetionRadius);
			REAL midLen = (steinerPoint - a->vert).length();
			steinerVert->relaxedInsetionRadius = std::max(midLen, steinerVert->relaxedInsetionRadius);
		}
	}
	else
		steinerVert->relaxedInsetionRadius = std::min((steinerPoint - a->vert).length(), (steinerPoint - b->vert).length());

	if (bothAcute) {
		Segment ns(b, steinerVert);
		Vertex *steinerVertOther = allocVertex(steinerPointOther, REAL(0), type);

		Tetrahedron toBeDeleted2;
		if (!meshRep.adjacent2SegmentFast(ns, &toBeDeleted2)) {
			meshRep.adjacent2Vertex(a, &toBeDeleted2);
			toBeDeleted2 = findPosition(steinerVertOther, toBeDeleted2, meshRep);
		}

		meshRep.deleteSegment(b, steinerVert);
		meshRep.addSegment(b, steinerVertOther);
		meshRep.addSegment(steinerVert, steinerVertOther);
		mayEncroachedSegs.push_back(Segment(steinerVert, steinerVertOther));
		mayEncroachedSegs.push_back(Segment(b, steinerVertOther));

		vifs.parent = NULL; vifs.rejected = false;
		insertVertexOnSegment(steinerVertOther, ns, toBeDeleted2, meshRep, vifs);

		if (vifs.insertRadiusTest) {
			steinerVertOther->relaxedInsetionRadius = (vifs.parent->vert - steinerPointOther).length();
			if (Adjacent(ns, vifs.parent))
				steinerVertOther->relaxedInsetionRadius = std::max(steinerVertOther->relaxedInsetionRadius,
					vifs.parent->relaxedInsetionRadius);
		}
		else
			steinerVertOther->relaxedInsetionRadius = std::min((steinerPointOther - b->vert).length(), 
			(steinerPointOther - steinerPoint).length());
	}

}

void DelMesher::refineSubPolygon(const Triangle &f, bool skinnyTetTestEnable) {
	Vertex *a = f.v[0], *b = f.v[1], *c = f.v[2], *d = NULL;
	int facetIndex = surfaceRep.getTriangleIndex(a, b, c);
	DelVector centerPoint;
	Geometer::Orthocircle(a->vert, a->weight, b->vert, b->weight, c->vert, c->weight, &centerPoint);
	Vertex *steinerVert = allocVertex(centerPoint, REAL(0), VertexType::Vertex_RefinedFreeFacet);

	VolumeVertexInsertionFlags vifs;
	vifs.cdt = true; vifs.encroachSegTest = true; vifs.encroachFacetTest = true;
	vifs.trueInsertion = false; vifs.insertRadiusTest = true;
	meshRep.Adjacent(f, &d);
	insertVertexOnSurface(steinerVert, f, Tetrahedron(d, a, b, c), meshRep, vifs);

	//detect encorachment
	bool encorach = false;
	Segment encorachedSeg;
	for (auto s : mayEncroachedSegs) {
		if (Encroached(s, steinerVert)) {
			encorach = true;
			encorachedSeg = s;
			break;
		}
	}
	mayEncroachedSegs.clear();

	if (!encorach){
		Assert(!vifs.rejected);
		//collect encorached volume vertices 
		std::vector<Vertex *> volumeVertices;
		volumeVertices.reserve(newFacesOfTets.size());
		for (auto f : newFacesOfTets) {
			for (int i = 0; i < 3; i++) {
				if (!matchOneOfVertexFlags(f.v[i]->getVertexType(), VertexType::Vertex_LowDimension)
					&& !f.v[i]->isMarked()) {
					f.v[i]->setMark();
					volumeVertices.push_back(f.v[i]);
				}
			}
		}
		for (auto v : volumeVertices) v->unSetMark();

		if (volumeVertices.size() == 0) {
			surfaceRep.bindVolumeVertex(steinerVert);
			for (auto s : newSegsOfFaces) surfaceRep.addTriangle(steinerVert, s.v[0], s.v[1], s.index);
			for (auto f : newFacesOfTets) {
				meshRep.addTetrahedron(steinerVert, f.v[0], f.v[1], f.v[2]);
				if (skinnyTetTestEnable) {
					Tetrahedron t(steinerVert, f.v[0], f.v[1], f.v[2]);
					if (skinnyTetTest(t)) skinnyTets.push(t);
				}
			}
			tobeDeletedFaces.clear(); newSegsOfFaces.clear();
			tobeDeletedTets.clear(); newFacesOfTets.clear();
		}
		else {
			for (auto f : tobeDeletedFaces) surfaceRep.addTriangle(f.v[0], f.v[1], f.v[2], f.index);
			for (auto t : tobeDeletedTets) meshRep.addTetrahedron(t.v[0], t.v[1], t.v[2], t.v[3]);
			tobeDeletedFaces.clear(); newSegsOfFaces.clear();
			tobeDeletedTets.clear(); newFacesOfTets.clear();

			VolumeVertexInsertionFlags dvifs;
			dvifs.skinnyTetTest = skinnyTetTestEnable;
			for (auto v : volumeVertices) deleteVertex(v, meshRep, dvifs);

			meshRep.Adjacent(f, &d);
			vifs.trueInsertion = true; vifs.parent = NULL; vifs.skinnyTetTest = skinnyTetTestEnable;
			insertVertexOnSurface(steinerVert, f, Tetrahedron(d, a, b, c), meshRep, vifs);

			Assert(!vifs.rejected);
		}
		//set relaxed insertion radius
		steinerVert->relaxedInsetionRadius = (centerPoint - vifs.parent->vert).length();
		if (Adjacent(facetIndex, vifs.parent)) {
			if (matchVertexFlag(vifs.parent->getVertexType(), VertexType::Vertex_Segment) &&
				steinerVert->relaxedInsetionRadius < SQRTF_2 * vifs.parent->relaxedInsetionRadius)
				steinerVert->relaxedInsetionRadius = SQRTF_2 * vifs.parent->relaxedInsetionRadius;
			else if (matchVertexFlag(vifs.parent->getVertexType(), VertexType::Vertex_Facet) &&
				steinerVert->relaxedInsetionRadius < vifs.parent->relaxedInsetionRadius)
				steinerVert->relaxedInsetionRadius = vifs.parent->relaxedInsetionRadius;
		}
	}
	else{
		//restore
		for (auto f : tobeDeletedFaces) surfaceRep.addTriangle(f.v[0], f.v[1], f.v[2], f.index);
		for (auto t : tobeDeletedTets) meshRep.addTetrahedron(t.v[0], t.v[1], t.v[2], t.v[3]);
		tobeDeletedFaces.clear(); newSegsOfFaces.clear();
		tobeDeletedTets.clear(); newFacesOfTets.clear();
		mayEncroachedFacets.clear();
		//set relaxed insertion radius
		steinerVert->relaxedInsetionRadius = (centerPoint - vifs.parent->vert).length();
		if (Adjacent(facetIndex, vifs.parent)) {
			if (matchVertexFlag(vifs.parent->getVertexType(), VertexType::Vertex_Segment) &&
				steinerVert->relaxedInsetionRadius < SQRTF_2 * vifs.parent->relaxedInsetionRadius)
				steinerVert->relaxedInsetionRadius = SQRTF_2 * vifs.parent->relaxedInsetionRadius;
			else if (matchVertexFlag(vifs.parent->getVertexType(), VertexType::Vertex_Facet) &&
				steinerVert->relaxedInsetionRadius < vifs.parent->relaxedInsetionRadius)
				steinerVert->relaxedInsetionRadius = vifs.parent->relaxedInsetionRadius;
		}
		//split encorached subsegment
		refineSubSegment(encorachedSeg, steinerVert, Adjacent(encorachedSeg, facetIndex), true, skinnyTetTestEnable);
		deallocVertex(steinerVert);
	}
}

void DelMesher::refineTetrahedron(const Tetrahedron& tet, bool& processed) {
	Vertex *a = tet.v[0], *b = tet.v[1], *c = tet.v[2], *d = tet.v[3];
	DelVector center;
	REAL radius = 0;
	Geometer::Orthosphere(a->vert, a->weight, b->vert, b->weight,
		c->vert, c->weight, d->vert, d->weight, &center, &radius);

	Vertex *steinerVert = meshRep.allocVertex(center, REAL(0), VertexType(Vertex_Free | Vertex_Refined));
	Tetrahedron found = findPositionWithOcclusion(steinerVert, tet, meshRep, surfaceRep);

	if (found.v[0]) {
		VolumeVertexInsertionFlags vifs;
		vifs.cdt = true; vifs.encroachSegTest = true; vifs.encroachFacetTest = true;
		vifs.trueInsertion = false;
		insertVertex(steinerVert, found, meshRep, vifs);

		//detect subsegment encorachment first
		bool encorach = false;
		Segment encorachedSeg;
		for (auto s : mayEncroachedSegs) {
			if (Encroached(s, steinerVert)) {
				encorach = true;
				encorachedSeg = s;
				break;
			}
		}
		mayEncroachedSegs.clear();

		if (encorach) {
			for (auto t : tobeDeletedTets) meshRep.addTetrahedron(t.v[0], t.v[1], t.v[2], t.v[3]);
			tobeDeletedTets.clear(); newFacesOfTets.clear();
			mayEncroachedFacets.clear();
			refineSubSegment(encorachedSeg, steinerVert, false, true, true);
			deallocVertex(steinerVert);
			processed = false;
		}
		else {
			//detect subpolygon encorachment next
			bool direct = false; Triangle encorachedFace;

			auto directEncroachmentTest = [](const DelVector& u, const DelVector& a, const DelVector& b, const DelVector& c) {
				constexpr Predicator<REAL> predicator;
				DelVector uProjected = Geometer::getProjectedPoint(u, a, b, c);
				auto plane = predicator.getProjectionPlane(a, b, c);
				bool ori0 = predicator.orientCoplane(uProjected, a, b, plane) >= 0;
				bool ori1 = predicator.orientCoplane(uProjected, b, c, plane) >= 0;
				bool ori2 = predicator.orientCoplane(uProjected, c, a, plane) >= 0;
				return ori0 == ori1 && ori0 == ori2;
			};

			for (auto f : mayEncroachedFacets) {
				if (Encroached(f, steinerVert)) {
					encorach = true;
					encorachedFace = f;
					//detect direct encorachment
					if (directEncroachmentTest(center, f.v[0]->vert, f.v[1]->vert, f.v[2]->vert)) {
						direct = true;
						break;
					}
				}
			}
			mayEncroachedFacets.clear();
			if (encorach) {
				for (auto t : tobeDeletedTets) meshRep.addTetrahedron(t.v[0], t.v[1], t.v[2], t.v[3]);
				tobeDeletedTets.clear(); newFacesOfTets.clear();
				if (direct) refineSubPolygon(encorachedFace, true);
				else {
					encorachedFace = findPositionWithOcclusion(steinerVert, encorachedFace);
					if (encorachedFace.v[0]) refineSubPolygon(encorachedFace, true);
					else refineSubSegment(Segment(encorachedFace.v[1], encorachedFace.v[2]), steinerVert, false, true, true);
				}
				deallocVertex(steinerVert);
				processed = false;
			}
			else {
				Assert(!vifs.rejected);
				steinerVert->relaxedInsetionRadius = radius;

				//reject vertex out of the convex hull
				bool hullOut = false;
				for (auto f : newFacesOfTets) {
					if (f.v[0]->isGhost() || f.v[1]->isGhost() || f.v[2]->isGhost()) {
						hullOut = true;
						break;
					}
				}

				if (!hullOut) {
					for (auto f : newFacesOfTets) {
						meshRep.addTetrahedron(steinerVert, f.v[0], f.v[1], f.v[2]);
						Tetrahedron t(steinerVert, f.v[0], f.v[1], f.v[2]);
						t.setRationAndRadius();
						if (t.getREration() > maxRatio || t.getRadius() > maxRadius)
							skinnyTets.push(t);
					}
				}
				else {
					meshRep.deallocVertex(steinerVert);
					for (auto t : tobeDeletedTets) meshRep.addTetrahedron(t.v[0], t.v[1], t.v[2], t.v[3]);
				}

				tobeDeletedTets.clear(); newFacesOfTets.clear();
				processed = true;
			}
		}
	}
	else {
		Triangle occluder(found.v[1], found.v[2], found.v[3]);
		bool encorach = false;
		Vertex *encorachedVert = NULL;
		for (int i = 0; i < 4; i++) {
			if (Encroached(occluder, tet.v[i])) {
				encorachedVert = tet.v[i];
				encorach = true;
				break;
			}
		}

		processed = false;
		if (encorach) {
			deallocVertex(steinerVert);
			occluder = findPositionWithOcclusion(encorachedVert, occluder);
			if (occluder.v[0]) refineSubPolygon(occluder, true);
			else {
				Segment s(occluder.v[1], occluder.v[2]);
				bool entwinement = false;
				if (matchVertexFlag(encorachedVert->getVertexType(), VertexType(VertexType::Vertex_Refined
					| VertexType::Vertex_Facet | VertexType::Vertex_Volume)))
					entwinement = Adjacent(s, encorachedVert);
				else encorachedVert = NULL;
				refineSubSegment(s, encorachedVert, entwinement, true, true);
			}
		}
		else {
			occluder = findPositionWithOcclusion(steinerVert, occluder);
			if (occluder.v[0] == NULL) {
				refineSubSegment(Segment(occluder.v[1], occluder.v[2]), steinerVert, false, true, true);
				deallocVertex(steinerVert);
			}
			else if (Encroached(occluder, steinerVert)) {
				deallocVertex(steinerVert);
				refineSubPolygon(occluder, true);
			}
			else { //fenced in
				deallocVertex(steinerVert);
				processed = true;
			}
		}
	}
}

//Ensure: also marking region vertices
void DelMesher::findMissingRegion(const Triangle& missed, std::vector<Vertex *> &regionVertices,
	std::vector<Segment>& regionBoundary, std::vector<Triangle>& regionFaces) {
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
	std::vector<Segment>& regionBoundary, std::vector<Triangle>& regionFaces, int depth) {
	if (findSegment(Segment(edge.v[0], edge.v[1], true))) {
		regionBoundary.push_back(Segment(edge.v[1], edge.v[0]));
		return;
	}

	Vertex *w = NULL;
	if (!surfaceRep.Adjacent(edge, &w)) return;

	Vertex *u = edge.v[0], *v = edge.v[1];
	Triangle f(u, v, w);
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

bool DelMesher::findCrossEdge(const Segment& boundary, const std::vector<Triangle>& regionFaces, Segment& cross) const {
	Vertex *a = boundary.v[0], *b = boundary.v[1], *c = NULL, *d = NULL;

	if (meshRep.isSegment(boundary)) {
		Tetrahedron t;
		meshRep.adjacent2SegmentFast(boundary, &t);
		c = t.v[2]; d = t.v[3];
		Assert(c != NULL && d != NULL);
	}
	else {
	    surfaceRep.Adjacent(Segment(b, a), &c);
		meshRep.Adjacent(Triangle(a, b, c), &d);
		std::swap(c, d);
		Assert(c != NULL && d != NULL);
	}

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
		bool found = meshRep.Adjacent(Triangle(a, b, c), &d);
		Assert(found);
		std::swap(c, d);
	} while (c != anchor);

	return hasIntersection;
}

//Require: all region vertices are marked
bool DelMesher::findCavity(const std::vector<Segment>& regionBoundaries, std::vector<Triangle>& regionFaces,
	std::vector<Vertex *>& positiveVertices, std::vector<Triangle>& positiveFaces,
	std::vector<Vertex *>& negativeVertices, std::vector<Triangle>& negativeFaces,
	std::vector<Tetrahedron>& deleted) {
	enum CavityTriangleType { TwoPositive, TwoNegative, TwoCoplanar, PositiveNegativeCoplanar, NegativePositiveCoplanar };
	std::stack<std::pair<Triangle, CavityTriangleType>> triangleStack;
	std::vector<Segment> newSegments;
	std::vector<Triangle> toBeErased;
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
		positiveFaces.push_back(Triangle(seg.v[1], seg.v[0], cross.v[0]));
		negativeVertices.push_back(cross.v[1]);
		negativeFaces.push_back(Triangle(seg.v[0], seg.v[1], cross.v[1]));

		triangleStack.push(std::make_pair(Triangle(cross.v[0], cross.v[1], seg.v[0]),
			CavityTriangleType::PositiveNegativeCoplanar));
		triangleStack.push(std::make_pair(Triangle(cross.v[1], cross.v[0], seg.v[1]),
			CavityTriangleType::NegativePositiveCoplanar));

		cross.v[0]->setMark(); cross.v[1]->setMark();
	}
	else {
		Vertex *posi = cross.v[0], *onFace = cross.v[1];
		toBeErased.push_back(Triangle(seg.v[0], seg.v[1], onFace));

		if (!surfaceRep.isMarked(onFace, seg.v[0])) {
			newSegments.push_back(Segment(onFace, seg.v[0]));
			triangleStack.push(std::make_pair(Triangle(onFace, seg.v[0], posi), CavityTriangleType::TwoCoplanar));
		}
		if (!surfaceRep.isMarked(seg.v[1], onFace)) {
			newSegments.push_back(Segment(seg.v[1], onFace));
			triangleStack.push(std::make_pair(Triangle(seg.v[1], onFace, posi), CavityTriangleType::TwoCoplanar));
		}
	}

	//find other cross tets
	while (!triangleStack.empty()) {
		auto pair = triangleStack.top();
		triangleStack.pop();

		Triangle intersect = pair.first;
		CavityTriangleType faceType = pair.second;

		Vertex *a = intersect.v[0], *b = intersect.v[1], *c = intersect.v[2], *d = NULL;
		if (!meshRep.Adjacent(intersect, &d)) continue;

		Assert(!d->isGhost());

		switch (faceType) {
		case CavityTriangleType::TwoPositive:
		{
			deleted.push_back(Tetrahedron(d, a, b, c));
			meshRep.deleteTetrahedron(d, a, b, c);
			Triangle intersected;
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
					triangleStack.push(std::make_pair(Triangle(a, b, d), CavityTriangleType::TwoPositive));
					triangleStack.push(std::make_pair(Triangle(c, d, b), CavityTriangleType::TwoNegative));
					triangleStack.push(std::make_pair(Triangle(d, c, a), CavityTriangleType::TwoNegative));
				}
				else {
					if (!d->isMarked()) {
						positiveVertices.push_back(d);
						d->setMark();
					}
					positiveFaces.push_back(Triangle(b, a, d));
					triangleStack.push(std::make_pair(Triangle(a, d, c), CavityTriangleType::TwoPositive));
					triangleStack.push(std::make_pair(Triangle(d, b, c), CavityTriangleType::TwoPositive));
				}
			}
			else {
				Assert(hasIntersection);
				positiveFaces.push_back(Triangle(b, a, d));
				triangleStack.push(std::make_pair(Triangle(b, c, d), CavityTriangleType::PositiveNegativeCoplanar));
				triangleStack.push(std::make_pair(Triangle(c, a, d), CavityTriangleType::NegativePositiveCoplanar));
			}
			break;
		}
		case CavityTriangleType::TwoNegative:
		{
			deleted.push_back(Tetrahedron(d, a, b, c));
			meshRep.deleteTetrahedron(d, a, b, c);
			Triangle intersected;
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
					triangleStack.push(std::make_pair(Triangle(a, b, d), CavityTriangleType::TwoNegative));
					triangleStack.push(std::make_pair(Triangle(c, d, b), CavityTriangleType::TwoPositive));
					triangleStack.push(std::make_pair(Triangle(d, c, a), CavityTriangleType::TwoPositive));
				}
				else {
					if (!d->isMarked()) {
						negativeVertices.push_back(d);
						d->setMark();
					}
					negativeFaces.push_back(Triangle(b, a, d));
					triangleStack.push(std::make_pair(Triangle(a, d, c), CavityTriangleType::TwoNegative));
					triangleStack.push(std::make_pair(Triangle(d, b, c), CavityTriangleType::TwoNegative));
				}
			}
			else {
				Assert(hasIntersection);
				negativeFaces.push_back(Triangle(b, a, d));
				triangleStack.push(std::make_pair(Triangle(b, c, d), CavityTriangleType::NegativePositiveCoplanar));
				triangleStack.push(std::make_pair(Triangle(c, a, d), CavityTriangleType::PositiveNegativeCoplanar));
			}
			break;
		}
		case CavityTriangleType::TwoCoplanar:
		{
			Triangle intersected;
			bool hasIntersection = false;
			do {
				for (auto f : regionFaces) {
					if (predicator.Intersection(f.v[0]->vert, f.v[1]->vert, f.v[2]->vert, c->vert, d->vert)) {
						intersected = f;
						hasIntersection = true;
						break;
					}
				}
				if (hasIntersection || !meshRep.Adjacent(Triangle(a, b, d), &c)) break;
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
					negativeFaces.push_back(Triangle(d, b, a));
					positiveFaces.push_back(Triangle(c, a, b));
					triangleStack.push(std::make_pair(Triangle(d, c, a), CavityTriangleType::NegativePositiveCoplanar));
					triangleStack.push(std::make_pair(Triangle(c, d, b), CavityTriangleType::PositiveNegativeCoplanar));
				}
				else {
					toBeErased.push_back(Triangle(b, a, d));
					if (!surfaceRep.isMarked(a, d)) {
						newSegments.push_back(Segment(a, d));
						triangleStack.push(std::make_pair(Triangle(a, d, c), CavityTriangleType::TwoCoplanar));
					}
					if (!surfaceRep.isMarked(d, b)) {
						newSegments.push_back(Segment(d, b));
						triangleStack.push(std::make_pair(Triangle(d, b, c), CavityTriangleType::TwoCoplanar));
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
			Triangle intersected;
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
				positiveFaces.push_back(Triangle(a, c, d));
				triangleStack.push(std::make_pair(Triangle(d, a, b), CavityTriangleType::TwoPositive));
				triangleStack.push(std::make_pair(Triangle(d, b, c), CavityTriangleType::PositiveNegativeCoplanar));
			}
			else if (ori < 0) {
				if (!d->isMarked()) {
					negativeVertices.push_back(d);
					d->setMark();
				}
				negativeFaces.push_back(Triangle(c, b, d));
				triangleStack.push(std::make_pair(Triangle(b, d, a), CavityTriangleType::TwoNegative));
				triangleStack.push(std::make_pair(Triangle(a, d, c), CavityTriangleType::PositiveNegativeCoplanar));
			}
			else {
				Assert(d->isMarked());
				//two vertices on face
				triangleStack.push(std::make_pair(Triangle(a, b, d), CavityTriangleType::PositiveNegativeCoplanar));
				if (!surfaceRep.isMarked(d, c)) {
					if (!surfaceRep.Contain(Segment(d, c))) newSegments.push_back(Segment(d, c));
					triangleStack.push(std::make_pair(Triangle(d, c, a), CavityTriangleType::TwoCoplanar));
				}
				positiveFaces.push_back(Triangle(a, c, d));
				negativeFaces.push_back(Triangle(c, b, d));
			}
			break;
		}
		case CavityTriangleType::NegativePositiveCoplanar:
		{
			deleted.push_back(Tetrahedron(d, a, b, c));
			meshRep.deleteTetrahedron(d, a, b, c);
			DelVector pp = a->vert, qq = b->vert, rr = d->vert;
			Triangle intersected;
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
				negativeFaces.push_back(Triangle(a, c, d));
				triangleStack.push(std::make_pair(Triangle(d, a, b), CavityTriangleType::TwoNegative));
				triangleStack.push(std::make_pair(Triangle(d, b, c), CavityTriangleType::NegativePositiveCoplanar));
			}
			else if (ori > 0) {
				if (!d->isMarked()) {
					positiveVertices.push_back(d);
					d->setMark();
				}
				positiveFaces.push_back(Triangle(c, b, d));
				triangleStack.push(std::make_pair(Triangle(b, d, a), CavityTriangleType::TwoPositive));
				triangleStack.push(std::make_pair(Triangle(a, d, c), CavityTriangleType::NegativePositiveCoplanar));
			}
			else {
				Assert(d->isMarked());
				//two vertices on face
				triangleStack.push(std::make_pair(Triangle(a, b, d), CavityTriangleType::NegativePositiveCoplanar));
				if (!surfaceRep.isMarked(c, d)) {
					if (!surfaceRep.Contain(Segment(c, d))) newSegments.push_back(Segment(c, d));
					triangleStack.push(std::make_pair(Triangle(c, d, b), CavityTriangleType::TwoCoplanar));
				}
				negativeFaces.push_back(Triangle(a, c, d));
				positiveFaces.push_back(Triangle(c, b, d));
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
		int faceIndex = surfaceRep.getTriangleIndex(regionFaces[0].v[0], regionFaces[0].v[1], regionFaces[0].v[2]);
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
					mayMissingFacets.push_back(f);
					regionFaces.push_back(f);
				}
			}
		}
		else {
			for (auto f : toBeErased) 
				surfaceRep.addTriangle(f.v[0], f.v[1], f.v[2], faceIndex);

			for (auto v : positiveVertices) v->unSetMark();
			for (auto v : negativeVertices) v->unSetMark();

			return false;
		}
	}

	for (auto v : positiveVertices) v->unSetMark();
	for (auto v : negativeVertices) v->unSetMark();

	return true;
}

bool DelMesher::facetRecovery(Triangle& f, std::vector<Vertex *>& regionVertices,
	std::vector<Segment>& regionBoundaries, std::vector<Triangle>& regionFaces,
	std::vector<Vertex *>& positiveVertices, std::vector<Triangle>& positiveFaces,
	std::vector<Vertex *>& negativeVertices, std::vector<Triangle>& negativeFaces,
	std::vector<Tetrahedron>& deleted, std::vector<Tetrahedron>& inserted) {

	Triangle seed = f;
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
	return success;
}


bool DelMesher::triangulateCavity(const std::vector<Triangle>& regionFaces, std::vector<Triangle>& boundaryFaces,
	std::vector<Vertex *>& cavityVertices, std::vector<Tetrahedron>& deleted, std::vector<Tetrahedron>& inserted, Triangle& encroached) {
	size_t cavityVerticesCount = cavityVertices.size();
	size_t boundaryFacesSize = boundaryFaces.size();

	std::vector<Triangle> newBoundaryFaces;
	std::vector<Triangle> shellFaces;
	shellFaces.reserve(boundaryFacesSize);
	newBoundaryFaces.reserve(boundaryFacesSize);

	TetMeshDataStructure cavityRep;
	cavityRep.Reserve(cavityVerticesCount);

	//counstruct vertices map
	for (auto& oldVert : cavityVertices) {
		Vertex *newVert = cavityRep.allocVertex(oldVert->vert, oldVert->weight);
		oldVert->setVertexPointer(newVert);
		newVert->setVertexPointer(oldVert);

		oldVert = newVert;
	}

	//construct triangulation
	triangulation3D(cavityVertices, cavityRep);

	int swapCount = 0;
	do {
		shellFaces.clear();
		bool enlarged = false;
		for (auto f : boundaryFaces) {
			Triangle face(f.v[0]->getPointedVertex(), f.v[1]->getPointedVertex(), f.v[2]->getPointedVertex());
			if (!cavityRep.Contain(face)) {
				Vertex *oppo = NULL;
				std::swap(f.v[1], f.v[2]);
				if (meshRep.Adjacent(f, &oppo)) {
					enlarged = true;
					meshRep.deleteTetrahedron(oppo, f.v[0], f.v[1], f.v[2]);
					deleted.push_back(Tetrahedron(oppo, f.v[0], f.v[1], f.v[2]));
					if (oppo->getPointedVertex() == NULL) {
						//alloc vertex
						Vertex *newOppo = cavityRep.allocVertex(oppo->vert, oppo->weight);
						oppo->setVertexPointer(newOppo);
						newOppo->setVertexPointer(oppo);

						//insert it
						Tetrahedron hint;
						bool found = cavityRep.adjacent2Vertex(face.v[0], &hint);
						Assert(found);
						hint = findPosition(newOppo, hint, cavityRep);
						insertVertex(newOppo, hint, cavityRep);

						cavityVertices.push_back(newOppo);
					}
					shellFaces.push_back(Triangle(oppo, f.v[0], f.v[2]));
					shellFaces.push_back(Triangle(oppo, f.v[1], f.v[0]));
					shellFaces.push_back(Triangle(oppo, f.v[2], f.v[1]));
				}

				Triangle back(f.v[0], f.v[2], f.v[1]);
				if (surfaceRep.Contain(f)) mayMissingFacets.push_front(f);
				if (surfaceRep.Contain(back)) mayMissingFacets.push_front(back);
			}
			else
				shellFaces.push_back(f);
		}

		if (!enlarged) break;
		std::swap(shellFaces, boundaryFaces);
		swapCount += 1;
	} while (true);

	bool success = true;
	for (auto f : regionFaces) {
		Triangle face(f.v[0]->getPointedVertex(), f.v[1]->getPointedVertex(), f.v[2]->getPointedVertex());
		//check region faces in cavity and mark it
		if (cavityRep.Contain(face))
			cavityRep.setMark(face.v[0], face.v[1], face.v[2]);
		else {
			//if not, find missing triangle edge put to v[1]-v[2]
			if (!surfaceRep.Contain(f)) std::swap(f.v[1], f.v[2]);
			for (int i = 0; i < 3; i++) {
				Vertex *a = f.v[i], *b = f.v[NEXT_F(i)], *c = f.v[NEXT_F(NEXT_F(i))];
				if (!findSegment(Segment(b, c))) {
					encroached = Triangle(a, b, c);
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
			propagateCleanCavity(Triangle(f.v[0]->getPointedVertex(), f.v[2]->getPointedVertex(), f.v[1]->getPointedVertex()), cavityRep, 0);
		//collect tets
		for (auto t : cavityRep)
			inserted.push_back(Tetrahedron(t.v[0]->getPointedVertex(), t.v[1]->getPointedVertex(), 
				t.v[2]->getPointedVertex(), t.v[3]->getPointedVertex()));
	}

	if (swapCount & 1) std::swap(shellFaces, boundaryFaces);

	for (auto v : cavityVertices) (v->getPointedVertex())->setVertexPointer(NULL);
	return success;
}

void DelMesher::propagateCleanCavity(const Triangle& f, TetMeshDataStructure& cavityRep, int depth) {
	Vertex *a = NULL;
	if (!cavityRep.Adjacent(f, &a) || cavityRep.isMarked(f.v[0], f.v[1], f.v[2]))
		return;

	Vertex *b = f.v[0], *c = f.v[1], *d = f.v[2];
	cavityRep.deleteTetrahedron(a, b, c, d);

	propagateCleanCavity(Triangle(a, c, d), cavityRep, depth + 1);
	propagateCleanCavity(Triangle(a, d, b), cavityRep, depth + 1);
	propagateCleanCavity(Triangle(a, b, c), cavityRep, depth + 1);
}

void DelMesher::refineRegion(const Triangle& regionFacet, bool missingSegTest, bool encroachedSegTest, bool setInsertionRadius) {
	Vertex *steinerVert = allocVertex((regionFacet.v[1]->vert + regionFacet.v[2]->vert) * REAL(0.5), 
		REAL(0), VertexType::Vertex_FreeFacet);

	SurfaceVertexInsertionFlags svifs;
	svifs.trueInsertion = false; svifs.deletionRecord = true;
	insertSurfaceVertex(steinerVert, regionFacet, svifs);

	//check encroachment
	bool encroach = false;
	Segment encrochedSeg;
	for (auto s : newSegsOfFaces) {
		if (meshRep.isSegment(s) && Encroached(s, steinerVert)) {
			encrochedSeg = s;
			encroach = true;
			break;
		}
	}

	if (!encroach) {
		Tetrahedron toBeDeleted;
		meshRep.adjacent2Vertex(regionFacet.v[1], &toBeDeleted);
		toBeDeleted = findPosition(regionFacet.v[1], toBeDeleted, meshRep);

		VolumeVertexInsertionFlags vifs;
		vifs.cdt = true; vifs.trueInsertion = false;
		vifs.missingSegTest = missingSegTest; vifs.missingFacetTest = true;
		vifs.encroachSegTest = encroachedSegTest; vifs.insertRadiusTest = setInsertionRadius;
		insertVertex(steinerVert, toBeDeleted, meshRep, vifs);

		for (auto s : newSegsOfFaces) {
			surfaceRep.addTriangle(steinerVert, s.v[0], s.v[1], s.index);
			mayMissingFacets.push_back(Triangle(steinerVert, s.v[0], s.v[1]));
		}

		if (setInsertionRadius) {
			steinerVert->relaxedInsetionRadius = (vifs.parent->vert - steinerVert->vert).length();
			if (Adjacent(newSegsOfFaces[0].index, steinerVert)) {
				if (matchVertexFlag(vifs.parent->getVertexType(), VertexType::Vertex_Segment) &&
					steinerVert->relaxedInsetionRadius < SQRTF_2 * vifs.parent->relaxedInsetionRadius)
					steinerVert->relaxedInsetionRadius = SQRTF_2 * vifs.parent->relaxedInsetionRadius;
				else if (matchVertexFlag(vifs.parent->getVertexType(), VertexType::Vertex_Facet) &&
					steinerVert->relaxedInsetionRadius < vifs.parent->relaxedInsetionRadius)
					steinerVert->relaxedInsetionRadius = vifs.parent->relaxedInsetionRadius;
			}
		}
		newSegsOfFaces.clear(); tobeDeletedFaces.clear();
	}
	else {
		//fall back
		for (auto f : tobeDeletedFaces) surfaceRep.addTriangle(f.v[0], f.v[1], f.v[2], f.index);
		newSegsOfFaces.clear(); tobeDeletedFaces.clear();

		deallocVertex(steinerVert);
		bool onSurface = surfaceRep.Contain(encrochedSeg) || surfaceRep.Contain(Segment(encrochedSeg.v[1], encrochedSeg.v[0]));
		VertexType type = VertexType::Vertex_FreeSegment;
		if (onSurface) type = VertexType(type | VertexType::Vertex_Facet);

		Vertex *a = encrochedSeg.v[0], *b = encrochedSeg.v[1];
		bool aIsAcute = matchVertexFlag(a->getVertexType(), VertexType::Vertex_Acute);
		bool bIsAcute = matchVertexFlag(b->getVertexType(), VertexType::Vertex_Acute);
		bool midSplit = true;
		if (aIsAcute != bIsAcute)
			midSplit = aIsAcute ? matchVertexFlag(b->getVertexType(), VertexType::Vertex_Free) :
			matchVertexFlag(a->getVertexType(), VertexType::Vertex_Free);

		//if both vertices are acute, just split in half
		if (midSplit)
			steinerVert = allocVertex((a->vert + b->vert) * REAL(0.5), REAL(0), type);
		else {
			//calulate steiner point for acute vertex
			DelVector ab = b->vert - a->vert;
			REAL r = nearestPowerOfTwo(REAL(0.5) * ab.length());
			if (aIsAcute) steinerVert = allocVertex(a->vert + r * ab, REAL(0), type);
			else steinerVert = allocVertex(b->vert - r * ab, REAL(0), type);
		}

		//split the encroched segment and update related information
		Segment s0(a, steinerVert), s1(b, steinerVert);
		meshRep.deleteSegment(a, b);
		meshRep.addSegment(s0.v[0], s0.v[1]);
		meshRep.addSegment(s1.v[0], s1.v[1]);
		if (missingSegTest) {
			mayMissingSegs.push_back(s0);
			mayMissingSegs.push_back(s1);
		}
		if (encroachedSegTest) {
			mayEncroachedSegs.push_back(s0);
			mayEncroachedSegs.push_back(s1);
		}

		if (matchVertexFlag(a->getVertexType(), VertexType::Vertex_FreeSegment))
			steinerVert->setOriSegmentIndex(a->getOriSegmentIndex());
		else if (matchVertexFlag(b->getVertexType(), VertexType::Vertex_FreeSegment))
			steinerVert->setOriSegmentIndex(b->getOriSegmentIndex());
		else {
			steinerVert->setOriSegmentIndex(oriSegments.size());
			oriSegments.push_back(encrochedSeg);
		}

		Tetrahedron toBeDeleted;
		bool found = meshRep.adjacent2SegmentFast(encrochedSeg, &toBeDeleted);
		Assert(found);

		VolumeVertexInsertionFlags newVifs;
		newVifs.cdt = true; newVifs.missingSegTest = true; newVifs.missingFacetTest = true;
		newVifs.insertRadiusTest = setInsertionRadius;
		insertVertexOnSegment(steinerVert, encrochedSeg, toBeDeleted, meshRep, newVifs);

		if (setInsertionRadius) {
			steinerVert->relaxedInsetionRadius = (newVifs.parent->vert - steinerVert->vert).length();
			if (midSplit) {
				if (Adjacent(encrochedSeg, newVifs.parent) || newVifs.parent == a || newVifs.parent == b)
					steinerVert->relaxedInsetionRadius = std::max(steinerVert->relaxedInsetionRadius,
						newVifs.parent->relaxedInsetionRadius);
			}
			else {
				if (Adjacent(encrochedSeg, newVifs.parent))
					steinerVert->relaxedInsetionRadius = std::max(steinerVert->relaxedInsetionRadius,
						newVifs.parent->relaxedInsetionRadius);
			}
		}
	}
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
		digCavity(u, Triangle(b, d, c, true), meshRep, vifs, rt);
		digCavity(u, Triangle(a, c, d, true), meshRep, vifs);
		digCavity(u, Triangle(a, d, b, true), meshRep, vifs);
		digCavity(u, Triangle(a, b, c, true), meshRep, vifs);
	}
	else {
		digCavity(u, Triangle(b, d, c, true), meshRep, vifs, rt);
		digCavity(u, Triangle(a, c, d, true), meshRep, vifs);
		digCavity(u, Triangle(a, d, b, true), meshRep, vifs);
		digCavity(u, Triangle(a, b, c, true), meshRep, vifs);

		if (!vifs.rejected) {
			if (vifs.trueInsertion) {
				for (auto f : newFacesOfTets) meshRep.addTetrahedron(u, f.v[0], f.v[1], f.v[2]);
				newFacesOfTets.clear();
			}
		}
		else {
			if (vifs.trueInsertion) {
				triangulateCavity(u, newFacesOfTets, meshRep, vifs);
				if (rt) meshRep.adjacent2Vertex(u, rt);
				newFacesOfTets.clear();
			}
		}
	}

	if (vifs.encroachSegTest || vifs.missingSegTest)
		for (auto s : markedSegments) meshRep.unSetMark(s.v[0], s.v[1]);

	if (vifs.encroachSegTest)
		for (auto s : markedSegments) mayEncroachedSegs.push_back(s);
	if (vifs.missingSegTest)
		for (auto s : markedSegments) mayMissingSegs.push_back(s);
}

void DelMesher::insertVertexOnSegment(Vertex *u, const Segment& s, const Tetrahedron& tet, TetMeshDataStructure& meshRep, const VolumeVertexInsertionFlags& vifs){
	bool onFacet = matchVertexFlag(u->getVertexType(), VertexType::Vertex_Facet);
	if (onFacet) {
		SurfaceVertexInsertionFlags svifs;
		svifs.missingTest = vifs.missingFacetTest;
		svifs.encroachmentTest = vifs.encroachFacetTest;
		svifs.trueInsertion = vifs.trueInsertion;
		if (vifs.cdt) svifs.trueInsertion = false;
		if (!vifs.trueInsertion) svifs.deletionRecord = true;
		insertSurfaceSegmentVertex(u, s, svifs);
	}
	insertVertex(u, tet, meshRep, vifs);
	if (onFacet && vifs.cdt) {
		if (vifs.trueInsertion) {
			for (auto seg : newSegsOfFaces) surfaceRep.addTriangle(u, seg.v[0], seg.v[1], seg.index);
			newSegsOfFaces.clear();
		}
	}
}

void DelMesher::insertVertexOnSurface(Vertex *u, const Triangle& f, const Tetrahedron& tet,
	TetMeshDataStructure& meshRep, const VolumeVertexInsertionFlags& vifs) {
	SurfaceVertexInsertionFlags svifs;
	svifs.missingTest = vifs.missingFacetTest;
	svifs.encroachmentTest = vifs.encroachFacetTest;
	if (vifs.cdt) svifs.trueInsertion = false;
	if (!vifs.trueInsertion) svifs.deletionRecord = true;
	insertSurfaceVertex(u, f, svifs);

	insertVertex(u, tet, meshRep, vifs);

	if (vifs.cdt) {
		if (vifs.trueInsertion) {
			for (auto seg : newSegsOfFaces) surfaceRep.addTriangle(u, seg.v[0], seg.v[1], seg.index);
			newSegsOfFaces.clear();
		}
	}
}

void DelMesher::triangulateCavity(Vertex *u, std::vector<Triangle>& boundaries, TetMeshDataStructure& meshRep, const VolumeVertexInsertionFlags& vifs) {
	std::vector<Vertex *> cavityVertices;
	cavityVertices.reserve(boundaries.size());

	TetMeshDataStructure cavityRep;
	cavityRep.Reserve(boundaries.size());
	Vertex *nu = cavityRep.allocVertex(u->vert, u->weight);
	u->setVertexPointer(nu); nu->setVertexPointer(u);
	cavityVertices.push_back(nu);

	std::vector<Triangle> cavityShells;
	cavityShells.reserve(boundaries.size());

	std::vector<Triangle> ghostShells;
	ghostShells.reserve(64);

	for (auto f : boundaries) {
		bool hasGhost = false;
		for (int i = 0; i < 3; i++) {
			if (f.v[i]->getPointedVertex() == NULL) {
				Vertex *newVert = NULL;
				if (!f.v[i]->isGhost()) {
					newVert = cavityRep.allocVertex(f.v[i]->vert, f.v[i]->weight);
					cavityVertices.push_back(newVert);
				}
				else newVert = cavityRep.getGhostVertex();

				f.v[i]->setVertexPointer(newVert); newVert->setVertexPointer(f.v[i]);
			}
			if (f.v[i]->isGhost()) hasGhost = true;
		}
		if (!hasGhost) cavityShells.push_back(f);
		else ghostShells.push_back(f);
	}

	int swapCount = 0;
	std::swap(cavityShells, boundaries);
	cavityShells.clear();
	swapCount += 1;

	triangulation3D(cavityVertices, cavityRep);

	do {
		bool enlarge = false;
		for (auto f : boundaries) {
			Triangle face(f.v[0]->getPointedVertex(), f.v[1]->getPointedVertex(), f.v[2]->getPointedVertex());
			Vertex *newOppo = NULL;
			if (cavityRep.Adjacent(face, &newOppo)) {
				Vertex *oppo = NULL;
				Triangle back(f.v[0], f.v[2], f.v[1]);
				if (meshRep.Adjacent(back, &oppo)) {
					if (!surfaceRep.Contain(f) && !surfaceRep.Contain(back) && !oppo->isGhost() && 
						predicator.inOrthoSpherePerturbed(oppo->vert, oppo->weight, f.v[0]->vert, f.v[0]->weight, f.v[2]->vert, f.v[2]->weight,
						f.v[1]->vert, f.v[1]->weight, newOppo->vert, newOppo->weight) > REAL(0)) {
						enlarge = true;

						meshRep.deleteTetrahedron(oppo, f.v[0], f.v[2], f.v[1]);
						if (!vifs.trueInsertion) tobeDeletedTets.push_back(Tetrahedron(oppo, f.v[0], f.v[2], f.v[1]));

						if (oppo->getPointedVertex() == NULL) {
							Vertex *newVert = cavityRep.allocVertex(oppo->vert, oppo->weight);
							newVert->setVertexPointer(oppo); oppo->setVertexPointer(newVert);
							cavityVertices.push_back(newVert);
							insertVertex(newVert, findPosition(newVert, Tetrahedron(newOppo, face.v[0], face.v[1], face.v[2]),
								cavityRep), cavityRep);
						}

						cavityShells.push_back(Triangle(oppo, back.v[0], back.v[2]));
						cavityShells.push_back(Triangle(oppo, back.v[1], back.v[0]));
						cavityShells.push_back(Triangle(oppo, back.v[2], back.v[1]));

						if (vifs.missingSegTest || vifs.encroachSegTest) {
							Segment s0(oppo, f.v[0]), s1(oppo, f.v[1]), s2(oppo, f.v[2]);
							if (meshRep.isSegment(s0) && !meshRep.testAndMark(s0.v[0], s0.v[1]))
								markedSegments.push_back(s0);
							if (meshRep.isSegment(s1) && !meshRep.testAndMark(s1.v[0], s1.v[1]))
								markedSegments.push_back(s1);
							if (meshRep.isSegment(s2) && !meshRep.testAndMark(s2.v[0], s2.v[1]))
								markedSegments.push_back(s2);
						}
					}
					else
						cavityShells.push_back(f);
				}
				else {
					if (vifs.missingFacetTest || vifs.encroachFacetTest) {
						if (surfaceRep.Contain(f)) {
							if (vifs.missingFacetTest) mayMissingFacets.push_back(f);
							if (vifs.encroachFacetTest) mayEncroachedFacets.push_back(f);
						}
						else if (surfaceRep.Contain(back)) {
							if (vifs.missingFacetTest) mayMissingFacets.push_back(back);
							if (vifs.encroachFacetTest) mayEncroachedFacets.push_back(back);
						}
					}
				}
			}
			else {
				Vertex *oppo = NULL;
				Triangle back(f.v[0], f.v[2], f.v[1]);
				if (meshRep.Adjacent(back, &oppo)) {
					Assert(!oppo->isGhost());
					enlarge = true;

					meshRep.deleteTetrahedron(oppo, f.v[0], f.v[2], f.v[1]);
					if (!vifs.trueInsertion) tobeDeletedTets.push_back(Tetrahedron(oppo, f.v[0], f.v[2], f.v[1]));

					if (oppo->getPointedVertex() == NULL) {
						Vertex *newVert = cavityRep.allocVertex(oppo->vert, oppo->weight);
						newVert->setVertexPointer(oppo); oppo->setVertexPointer(newVert);
						cavityVertices.push_back(newVert);
						Tetrahedron hint;
						bool found = cavityRep.adjacent2Vertex(face.v[0], &hint);
						Assert(found);
						hint = findPosition(newVert, hint, cavityRep);
						insertVertex(newVert, hint, cavityRep);
					}

					cavityShells.push_back(Triangle(oppo, back.v[0], back.v[2]));
					cavityShells.push_back(Triangle(oppo, back.v[1], back.v[0]));
					cavityShells.push_back(Triangle(oppo, back.v[2], back.v[1]));

					if (vifs.missingSegTest || vifs.encroachSegTest) {
						Segment s0(oppo, f.v[0]), s1(oppo, f.v[1]), s2(oppo, f.v[2]);
						if (meshRep.isSegment(s0) && !meshRep.testAndMark(s0.v[0], s0.v[1]))
							markedSegments.push_back(s0);
						if (meshRep.isSegment(s1) && !meshRep.testAndMark(s1.v[0], s1.v[1]))
							markedSegments.push_back(s1);
						if (meshRep.isSegment(s2) && !meshRep.testAndMark(s2.v[0], s2.v[1]))
							markedSegments.push_back(s2);
					}
				}
				else {
					if (vifs.missingFacetTest || vifs.encroachFacetTest) {
						if (surfaceRep.Contain(f)) {
							if (vifs.missingFacetTest) mayMissingFacets.push_back(f);
							if (vifs.encroachFacetTest) mayEncroachedFacets.push_back(f);
						}
						else if (surfaceRep.Contain(back)) {
							if (vifs.missingFacetTest) mayMissingFacets.push_back(back);
							if (vifs.encroachFacetTest) mayEncroachedFacets.push_back(back);
						}
					}
				}
			}
		}

		if (!enlarge) break;
		std::swap(boundaries, cavityShells);
		cavityShells.clear();
		swapCount += 1;
	} while (true);

	//mark shell
	for (auto f : ghostShells) 
		cavityRep.setMark(f.v[0]->getPointedVertex(), f.v[1]->getPointedVertex(), f.v[2]->getPointedVertex());
	for (auto f : cavityShells)
		cavityRep.setMark(f.v[0]->getPointedVertex(), f.v[1]->getPointedVertex(), f.v[2]->getPointedVertex());
	//clean tets out of cavity
	for (auto f : cavityShells)
		propagateCleanCavity(Triangle(f.v[0]->getPointedVertex(), f.v[2]->getPointedVertex(), f.v[1]->getPointedVertex()), cavityRep, 0);
	//collect tets
	for (auto t : cavityRep) {
		meshRep.addTetrahedron(t.v[0]->getPointedVertex(), t.v[1]->getPointedVertex(),
			t.v[2]->getPointedVertex(), t.v[3]->getPointedVertex());
		if (vifs.skinnyTetTest) {
			Tetrahedron newTet(t.v[0]->getPointedVertex(), t.v[1]->getPointedVertex(),
				t.v[2]->getPointedVertex(), t.v[3]->getPointedVertex());
			if (skinnyTetTest(newTet)) skinnyTets.push(newTet);
		}
	}

	if (vifs.encroachFacetTest) {
		for (auto f : cavityShells) {
			Triangle back(f.v[0], f.v[2], f.v[1]);
			if (surfaceRep.Contain(f)) mayEncroachedFacets.push_back(f);
			if (surfaceRep.Contain(back)) mayEncroachedFacets.push_back(back);
		}
	}

	for (auto v : cavityVertices) (v->getPointedVertex())->setVertexPointer(NULL);
	if (ghostShells.size() > 0) meshRep.getGhostVertex()->setVertexPointer(NULL);

	if (swapCount % 2) std::swap(boundaries, cavityShells);
}

void DelMesher::deleteVertex(Vertex *u, TetMeshDataStructure& meshRep, const VolumeVertexInsertionFlags& vifs) {
	Assert(!matchOneOfVertexFlags(u->getVertexType(), VertexType::Vertex_LowDimension));
	TetMeshDataStructure cavityRep;
	std::vector<Triangle> boundary;
	boundary.reserve(64);

	makeHole(u, boundary, meshRep);

	//collect vertices to be inserted
	std::vector<Vertex *> cavityVertices;
	cavityVertices.reserve(boundary.size());

	bool hasGhost = false;
	for (auto f : boundary) {
		bool hasGhost = true;
		for (int i = 0; i < 3; i++) {
			if (f.v[i]->getPointedVertex() == NULL) {
				Vertex *cavityVert = NULL;
				if (!f.v[i]->isGhost()) {
					cavityVert = cavityRep.allocVertex(f.v[i]->vert, f.v[i]->weight);
					cavityVertices.push_back(cavityVert);
				}
				else {
					cavityVert = cavityRep.getGhostVertex();
					hasGhost = true;
				}

				cavityVert->setVertexPointer(f.v[i]); 
				f.v[i]->setVertexPointer(cavityVert);
			}
		}

	}

	cavityRep.Reserve(cavityVertices.size());
	triangulation3D(cavityVertices, cavityRep);

	for (auto f : boundary) cavityRep.setMark(f.v[0]->getPointedVertex(), f.v[1]->getPointedVertex(), f.v[2]->getPointedVertex());

	//clean tets out of cavity
	Triangle boundaryFaceOut = Triangle(boundary[0].v[0]->getPointedVertex(), boundary[0].v[2]->getPointedVertex(),
		boundary[0].v[1]->getPointedVertex());
	propagateCleanCavity(boundaryFaceOut, cavityRep, 0);

	for (auto t : cavityRep) {
		meshRep.addTetrahedron(t.v[0]->getPointedVertex(), t.v[1]->getPointedVertex(),
			t.v[2]->getPointedVertex(), t.v[3]->getPointedVertex());

		if (vifs.skinnyTetTest) {
			Tetrahedron newTet(t.v[0]->getPointedVertex(), t.v[1]->getPointedVertex(),
				t.v[2]->getPointedVertex(), t.v[3]->getPointedVertex());
			if (skinnyTetTest(newTet)) skinnyTets.push(newTet);
		}
	}
	
	for (auto v : cavityVertices) (v->getPointedVertex())->setVertexPointer(NULL);
	if (hasGhost) meshRep.getGhostVertex()->setVertexPointer(NULL);
}

void DelMesher::insertSurfaceVertex(Vertex *u, const Triangle& f, const SurfaceVertexInsertionFlags& vifs){
	Vertex *a = f.v[0], *b = f.v[1], *c = f.v[2];
	int faceIndex = surfaceRep.getTriangleIndex(a, b, c);
	if (vifs.trueInsertion) surfaceRep.bindVolumeVertex(u);
	if (vifs.deletionRecord) tobeDeletedFaces.push_back(TriangleWithIndex(a, b, c, faceIndex));
	DelVector aboveVert = Geometer::calculateAbovePoint(a->vert, b->vert, c->vert);
	surfaceRep.deleteTriangle(a, b, c);
	digCavity(u, aboveVert, Segment(b, a), faceIndex, vifs);
	digCavity(u, aboveVert, Segment(c, b), faceIndex, vifs);
	digCavity(u, aboveVert, Segment(a, c), faceIndex, vifs);
}

void DelMesher::insertSurfaceSegmentVertex(Vertex *u, const Segment &s, const SurfaceVertexInsertionFlags& vifs) {
	if (vifs.trueInsertion) surfaceRep.bindVolumeVertex(u);
	Vertex *a = s.v[0], *b = s.v[1];
	Vertex *oppo = NULL;
	int index = -1;
	if (surfaceRep.Adjacent(s, &oppo, &index)) {
		surfaceRep.deleteTriangle(a, b, oppo);
		DelVector above = Geometer::calculateAbovePoint(a->vert, b->vert, oppo->vert);
		digCavity(u, above, Segment(a, oppo), index, vifs);
		digCavity(u, above, Segment(oppo, b), index, vifs);
	}
	if (surfaceRep.Adjacent(Segment(b, a), &oppo, &index)) {
		surfaceRep.deleteTriangle(b, a, oppo);
		DelVector above = Geometer::calculateAbovePoint(b->vert, a->vert, oppo->vert);
		digCavity(u, above, Segment(b, oppo), index, vifs);
		digCavity(u, above, Segment(oppo, a), index, vifs);
	}
}

//Also mark all vertices
void DelMesher::detectCoplanarFaces(const Triangle& f, REAL facetRadianTol,
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

void DelMesher::propagateClean(const Triangle &f, int depth){
	if (depth < 100) {
		Vertex *a = NULL;
		if (surfaceRep.Contain(Triangle(f.v[0], f.v[2], f.v[1])) || !meshRep.Adjacent(f, &a))
			return;

		Vertex *b = f.v[0], *c = f.v[1], *d = f.v[2];

		meshRep.deleteTetrahedron(a, b, c, d);

		propagateClean(Triangle(a, c, d), depth + 1);
		propagateClean(Triangle(a, b, c), depth + 1);
		propagateClean(Triangle(a, d, b), depth + 1);
	}
	else {
		std::stack<Triangle> workingStacks;
		workingStacks.push(f);

		while (!workingStacks.empty()) {
			Triangle top = workingStacks.top();
			workingStacks.pop();

			Vertex *a = NULL;
			if (!surfaceRep.Contain(Triangle(top.v[0], top.v[2], top.v[1])) && meshRep.Adjacent(top, &a)) {
				Vertex *b = top.v[0], *c = top.v[1], *d = top.v[2];

				meshRep.deleteTetrahedron(a, b, c, d);

				workingStacks.push(Triangle(a, c, d));
				workingStacks.push(Triangle(a, b, c));
				workingStacks.push(Triangle(a, d, b));
			}
		}
	}
}

void DelMesher::makeHole(Vertex *u, std::vector<Triangle>& boundaries, TetMeshDataStructure& meshRep) const {
	Tetrahedron t;
	bool found = meshRep.adjacent2Vertex(u, &t);
	Assert(found);

	Vertex *a = t.v[1], *b = t.v[2], *c = t.v[3];
	boundaries.push_back(Triangle(a, b, c));
	meshRep.deleteTetrahedron(u, a, b, c);

	propagateMakeHole(Triangle(u, a, b), boundaries, meshRep, 0);
	propagateMakeHole(Triangle(u, b, c), boundaries, meshRep, 0);
	propagateMakeHole(Triangle(u, c, a), boundaries, meshRep, 0);
}

void DelMesher::propagateMakeHole(const Triangle& f, std::vector<Triangle>& boundaries, TetMeshDataStructure& meshRep, int depth) const {
	if (depth < 100) {
		Vertex *u = f.v[0], *b = f.v[1], *c = f.v[2];
		Vertex *a = NULL;

		if (!meshRep.Adjacent(f, &a)) return;
		boundaries.push_back(Triangle(a, c, b));
		meshRep.deleteTetrahedron(a, u, b, c);

		propagateMakeHole(Triangle(u, a, c), boundaries, meshRep, depth + 1);
		propagateMakeHole(Triangle(u, b, a), boundaries, meshRep, depth + 1);
	}
	else {
		std::stack<Triangle> workingStacks;
		workingStacks.push(f);

		while (!workingStacks.empty()) {
			Triangle top = workingStacks.top();
			workingStacks.pop();

			Vertex *u = top.v[0], *b = top.v[1], *c = top.v[2];
			Vertex *a = NULL;
			if (meshRep.Adjacent(top, &a)) {
				boundaries.push_back(Triangle(a, c, b));
				meshRep.deleteTetrahedron(a, u, b, c);

				workingStacks.push(Triangle(u, a, c));
				workingStacks.push(Triangle(u, b, a));
			}
		}
	}

}

bool DelMesher::Encroached(const Segment& s) const {
	Vertex *a = s.v[0], *b = s.v[1], *c = NULL;
	auto begin = meshRep.getIncidentFacets(s);
	if (begin == NULL) return true;

	constexpr REAL epsilon = REAL(1e-16);
	DelVector center;
	REAL r = 0;
	Geometer::Orthocenter(a->vert, a->weight, b->vert, b->weight, &center, &r);
	REAL r2 = r * r;

	Assert(!isnan(r2));
	auto cirulator = begin;
	do {
		c = cirulator->v[2];
		if (!c->isGhost() && (c->vert - center).length2() - c->weight - r2 <= epsilon * r2)
			return true;

		++cirulator;
	} while (cirulator != begin);

	return false;
}

bool DelMesher::Encroached(const Triangle& f, Vertex **encroachedVert) const {
	int faceIndex = surfaceRep.getTriangleIndex(f.v[0], f.v[1], f.v[2]);
	Triangle f2(f.v[0], f.v[2], f.v[1], true);
	Vertex *x, *y;
	if (!meshRep.Adjacent(f, &x) || !meshRep.Adjacent(f2, &y))
		return true;

	constexpr REAL epsilon = REAL(1e-8);
	DelVector center;
	REAL r = 0;
	Geometer::Orthocircle(f.v[0]->vert, f.v[0]->weight, f.v[1]->vert, f.v[1]->weight, f.v[2]->vert, f.v[2]->weight, &center, &r);
	REAL r2 = r*r;

	DelVector n;
	bool hasNormal = false;
	bool encroached0 = false, encroached1 = false;
	if (!x->isGhost() && ((x->vert - center).length2() - x->weight - r2) <= epsilon * epsilon * r2) {
		//determine whether the encroached vertex is eligible or not
		if (Adjacent(faceIndex, x)) {
			n = Normalize(Geometer::triangleNormal(f.v[0]->vert, f.v[1]->vert, f.v[2]->vert));
			hasNormal = true;
			REAL dis = fabs((x->vert - center) * n);
			if (matchVertexFlag(x->getVertexType(), VertexType::Vertex_Segment))
				encroached0 = (dis - SQRTF_2 * x->relaxedInsetionRadius >= epsilon * dis);
			else
				encroached0 = (dis - x->relaxedInsetionRadius >= epsilon * dis);
		}
		else {
			encroached0 = true;
			*encroachedVert = x;
		}
	}
	if (!encroached0 && !y->isGhost() && ((y->vert - center).length2() - y->weight - r2) <= epsilon * epsilon * r2) {
		//determine whether the encroached vertex is eligible or not
		if (Adjacent(faceIndex, y)) {
			if(!hasNormal) n = Normalize(Geometer::triangleNormal(f.v[0]->vert, f.v[1]->vert, f.v[2]->vert));
			REAL dis = fabs((x->vert - center) * n);
			if (matchVertexFlag(x->getVertexType(), VertexType::Vertex_Segment))
				encroached1 = (dis - SQRTF_2 * y->relaxedInsetionRadius >= epsilon * dis);
			else
				encroached1 = (dis -  y->relaxedInsetionRadius >= epsilon * dis);
		}
		else {
			encroached1 = true;
			*encroachedVert = y;
		}
	}

	return encroached0 || encroached1;
}

bool DelMesher::Encroached(const Segment &s, Vertex *v) const {
	DelVector center;
	REAL r = 0;
	Geometer::Orthocenter(s.v[0]->vert, s.v[0]->weight, s.v[1]->vert, s.v[1]->weight, &center, &r);
	REAL r2 = r * r;

	constexpr REAL epsilon = REAL(1e-8);
	if ((v->vert - center).length2() - v->weight - r2 <= epsilon * epsilon * r2)
		return true;

	return false;
}

bool DelMesher::Encroached(const Triangle &f, Vertex *v) const {
	DelVector center;
	REAL r;
	Geometer::Orthocircle(f.v[0]->vert, f.v[0]->weight, f.v[1]->vert, f.v[1]->weight, f.v[2]->vert, f.v[2]->weight, &center, &r);
	REAL r2 = r * r;

	constexpr REAL epsilon = REAL(1e-8);
	if (!v->isGhost() && fabs((v->vert - center).length2() - v->weight - r2) <= epsilon * epsilon * r2) {
		//determine whether the encroached vertex is eligible or not
		if (Adjacent(surfaceRep.getTriangleIndex(f.v[0], f.v[1], f.v[2]), v)) {
			DelVector n = Normalize(Geometer::triangleNormal(f.v[0]->vert, f.v[1]->vert, f.v[2]->vert));
			REAL dis = fabs((v->vert - center) * n);
			if (matchVertexFlag(v->getVertexType(), VertexType::Vertex_Segment))
				return dis - SQRTF_2 * v->relaxedInsetionRadius >= epsilon * dis;
			else
				return dis - v->relaxedInsetionRadius >= epsilon * dis;
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
			adjacent = (seg.v[0]->isMarked() != seg.v[1]->isMarked());
		}
		else
			adjacent = (v != oriSegment.v[0] && v!= oriSegment.v[1] &&
				(findSegment(Segment(v, oriSegment.v[0])) || findSegment(Segment(v, oriSegment.v[1]))));
	}
	else if (matchVertexFlag(v->getVertexType(), VertexType::Vertex_Facet)) {
		Triangle facet;
		int facetIndex = -1;
		surfaceRep.adjacent2Vertex(v, &facet, &facetIndex);
		Vertex **verts = NULL;
		size_t size = getPolygonVertices(facetIndex, &verts);

		int markedCount = 0;
		for (size_t i = 0; i < size; i++) {
			if (verts[i]->isMarked()) 
				markedCount += 1;
		}
		adjacent = (markedCount == 1);
	}
	oriSegment.v[0]->unSetMark(); oriSegment.v[1]->unSetMark();

	return adjacent;
}

bool DelMesher::Adjacent(int fIndex, Vertex *v) const {
	Vertex **verts = NULL;
	size_t size = getPolygonVertices(fIndex, &verts);

	bool adjacent = false;
	if (matchVertexFlag(v->getVertexType(), VertexType::Vertex_Segment)) {
		if (matchVertexFlag(v->getVertexType(), VertexType::Vertex_Free)) {
			Segment oriSeg = oriSegments[v->getOriSegmentIndex()];
			oriSeg.v[0]->setMark(); oriSeg.v[1]->setMark();

			int markedCount = 0;
			for (size_t i = 0; i < size; i++) {
				if (verts[i]->isMarked())
					markedCount += 1;
			}
			adjacent = (markedCount == 1);

			//clean up
			oriSeg.v[0]->unSetMark(); oriSeg.v[1]->unSetMark();
		}
		else {
			for (size_t i = 0; i < size; i++) {
				if (verts[i] == v) {
					adjacent = false;
					break;
				}
				if (!adjacent && findSegment(Segment(v, verts[i])))
					adjacent = true;
			}
		}
	}
	else if (matchVertexFlag(v->getVertexType(), VertexType::Vertex_Facet)) {
		Triangle f2;
		int f2Index = -1;
		surfaceRep.adjacent2Vertex(v, &f2, &f2Index);

		if (fIndex == f2Index) return false;

		Vertex **verts2 = NULL;
		size_t size2 = getPolygonVertices(f2Index, &verts2);

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

bool DelMesher::Adjacent(const Segment &s, int facetIndex) const {
	Segment oriSegment = s;
	if (matchVertexFlag(s.v[0]->getVertexType(), VertexType::Vertex_FreeSegment)) 
		oriSegment = oriSegments[s.v[0]->getOriSegmentIndex()];
	else if (matchVertexFlag(s.v[1]->getVertexType(), VertexType::Vertex_FreeSegment)) 
		oriSegment = oriSegments[s.v[1]->getOriSegmentIndex()];
	oriSegment.v[0]->setMark(); oriSegment.v[1]->setMark();

	bool adjacent = false;
	Vertex **polyVerts;
	size_t size = getPolygonVertices(facetIndex, &polyVerts);

	int markedCount = 0;
	for (size_t i = 0; i < size; i++) {
		if (polyVerts[i]->isMarked())
			markedCount += 1;
	}
	adjacent = (markedCount == 1);

	oriSegment.v[0]->unSetMark(); oriSegment.v[1]->unSetMark();
	return adjacent;
}


Reference<Mesh> DelMesher::generateMesh(int *vertexLableMap){
	for (auto t : meshRep) {
		//set vertices relaxed insertion radius
		Vertex *a = t.v[0], *b = t.v[1], *c = t.v[2], *d = t.v[3];
		if (!a->isGhost() && !b->isGhost() && !c->isGhost() && !d->isGhost()) {
			DelVector ab = b->vert - a->vert, ac = c->vert - a->vert, ad = d->vert - a->vert,
				bc = c->vert - b->vert, bd = d->vert - b->vert, cd = d->vert - c->vert;
			REAL abLen = ab.length(), acLen = ac.length(), adLen = ad.length(),
				bcLen = bc.length(), bdLen = bd.length(), cdLen = cd.length();

			a->relaxedInsetionRadius = std::min({ a->relaxedInsetionRadius, abLen, acLen, adLen });
			b->relaxedInsetionRadius = std::min({ b->relaxedInsetionRadius, abLen, bcLen, bdLen });
			c->relaxedInsetionRadius = std::min({ c->relaxedInsetionRadius, acLen, bcLen, cdLen });
			d->relaxedInsetionRadius = std::min({ d->relaxedInsetionRadius, adLen, bdLen, cdLen });
		}
	}

	constrainedRefinment();

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
		for (std::vector<Vertex *>::size_type i = 0; i < oriVertices.size(); i++)
			vertexLableMap[i] = vertIndices[vi[oriVertices[i]]];
	}

	freeAligned(vertIndices);

	return NULL;
}
}
