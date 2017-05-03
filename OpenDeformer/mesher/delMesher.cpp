#include "stdafx.h"
#include "delMesher.h"
#include "geometer.h"
#include "meshRelabeler.h"
#include "tetmesh.h"
#include <map>
#include <random>
#include <stack>
#include <unordered_map>

namespace ODER{

Predicator<DelReal> DelTriangulator::predicator;

void DelTriangulator::generateSubPolygons(Vertex **vertices, int *segments, int vertexCount, int segmentCount, 
	const DelVector& above, bool boundaryOnly, std::vector<Triangle>& output){
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
		Vertex *newVert = meshRep.allocVertex(oriVert->point, oriVert->weight);
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

	ghost->point = above;
	v0 = vertices[0]->getPointedVertex();
	v1 = vertices[1]->getPointedVertex();
	v2 = NULL;

	for (int i = 2; i < vertexCount; i++) {
		v2 = vertices[i]->getPointedVertex();
		DelReal ori = predicator.orient2d(v0->point, v1->point, v2->point, above);
		if (ori != 0) {
			std::swap(vertices[2], vertices[i]);
			if (ori < 0) {
				std::swap(vertices[1], vertices[2]);
				std::swap(v1, v2);
			}
			break;
		}
	}
	Assert(v2);

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
	const DelVector& above, bool boundaryOnly, std::vector<Triangle>& output) {
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
		Vertex *newVert = meshRep.allocVertex(oriVert->point, oriVert->weight);
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

	ghost->point = above;
	v0 = vertices[0]->getPointedVertex();
	v1 = vertices[1]->getPointedVertex();
	v2 = NULL;

	for (int i = 2; i < vertexCount; i++) {
		v2 = vertices[i]->getPointedVertex();
		DelReal ori = predicator.orient2d(v0->point, v1->point, v2->point, above);
		if (ori != 0) {
			std::swap(vertices[2], vertices[i]);
			if (ori < 0) {
				std::swap(vertices[1], vertices[2]);
				std::swap(v1, v2);
			}
			break;
		}
	}
	Assert(v2);

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

void DelTriangulator::insertSegments(const Triangle *triangles, const Segment *segments, int triangleCount, int segmentCount, 
	const DelVector& above, std::vector<Triangle>& output) {
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

	for (int i = 0; i < triangleCount; i++) {
		Triangle f = triangles[i];
		for (int j = 0; j < 3; j++) {
			if (f.v[j]->getPointedVertex() == NULL) {
				Vertex *newVert = meshRep.allocVertex(f.v[j]->point, f.v[j]->weight);
				f.v[j]->setVertexPointer(newVert);
				newVert->setVertexPointer(f.v[j]);
			}
		}
		meshRep.addTriangle(f.v[0]->getPointedVertex(), f.v[1]->getPointedVertex(), f.v[2]->getPointedVertex());
	}
	//calculate above point
	ghost = meshRep.getGhostVertex();
	ghost->point = above;

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
		if (predicator.inHalfSpace2D(u->point, b->point, c->point, ghost->point))
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
		bool ori0 = predicator.orient2d(b->point, a->point, u->point, ghost->point) > 0;
		bool ori1 = predicator.orient2d(a->point, c->point, u->point, ghost->point) > 0;

		int condition = ori0 + (ori1 << 1);
		switch (condition){
		case 0:
			if (predicator.orient2d(c->point, b->point, u->point, ghost->point) > 0){
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
	DelVector aa = a->point, bb = b->point, above = ghost->point;
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
			if (predicator.orient2d(bb, f.v[1]->point, f.v[2]->point, above) < 0) {
				DelReal ori1 = predicator.orient2d(aa, bb, f.v[1]->point, above);
				DelReal ori2 = predicator.orient2d(aa, bb, f.v[2]->point, above);
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
		DelReal ori = predicator.orient2d(aa, bb, opposite->point, above);
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
	DelReal *distance = new DelReal[size];
	cavityRep.Clear();

	//costruct vertices map
	std::vector<Vertex *> cavityVertices;
	cavityVertices.reserve(size);
	for (auto vert : vertices) {
		Vertex *cavityVert = cavityRep.allocVertex(vert->point, vert->weight);
		cavityVert->setVertexPointer(vert);
		cavityVertices.push_back(cavityVert);
	}

	Vertex *start = cavityVertices[size - 1], *end = cavityVertices[0];
	prev[0] = size - 1; prev[size - 1] = size - 2;
	next[0] = 1; next[size - 1] = 0;
	distance[0] = distance[size - 1] = 0;

	DelVector ab = end->point - start->point;
	for (int i = 1; i < size - 1; i++) {
		prev[i] = i - 1;
		next[i] = i + 1;
		distance[i] = ((cavityVertices[i]->point - start->point) % ab).length2();
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
		outCircle = predicator.inCirclePerturbed(u->point, v->point, w->point, x->point, ghost->point) > 0;

		if (oriTest) negativeOri = predicator.orient2d(u->point, v->point, w->point, ghost->point) <= 0;
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

	if (exited)
		deleted = x->isMarked() && 
		predicator.inCirclePerturbed(u->point, v->point, w->point, x->point, ghost->point) > 0;

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

	//restore convex poly
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
		deleted = predicator.inHalfSpace2D(u->point, w->point, v->point, ghost->point);
	}
	else{
		deleted = predicator.inCirclePerturbed(u->point, verts[1]->point, 
			                               verts[0]->point, verts[2]->point, ghost->point) > 0;
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

Predicator<DelReal> DelMesher::predicator;

DelMesher::DelMesher(Vector3f *surfvs, int *triangles, int numv, int numtri, DelMesherCriteria cri) {
	criteria = cri;

	oriVertices.reserve(numv);
	TriMeshDataStructure oriSurface;
	for (int i = 0; i < numv; i++) {
		boundBox.Insert(surfvs[i]);
		Vertex *vert = oriSurface.allocVertex(DelVector{ surfvs[i].x, surfvs[i].y, surfvs[i].z }, DelReal(0));
		oriVertices.push_back(vert);
	}
	boundBox.Scale(DelReal(1.05));

	std::vector<Triangle> oriTriangles;
	oriTriangles.reserve(numtri);
	for (int i = 0; i < numtri; i++) {
		oriSurface.addTriangle(oriVertices[triangles[3 * i]], oriVertices[triangles[3 * i + 1]], oriVertices[triangles[3 * i + 2]]);
		oriTriangles.push_back(Triangle(oriVertices[triangles[3 * i]], oriVertices[triangles[3 * i + 1]], oriVertices[triangles[3 * i + 2]]));
	}

	DelReal facetRadianTol = criteria.facetAngleTol / DelReal(180) * DelReal(M_PI);
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
						Vertex *newVert = meshRep.allocVertex(s.v[i]->point, s.v[i]->weight, 
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
					Vertex *newVert = meshRep.allocVertex(v->point, v->weight, Vertex_FixedFacet);
					surfaceRep.bindVolumeVertex(newVert);
					v->setVertexPointer(newVert);
				}
			}

			//re-triangulation
			//all boundary facets orient inward
			std::swap(f.v[0], f.v[1]);
			if (coplanarVertices.size() > 3) {
				DelVector above = calculateAbovePoint(coplanarVertices, f);
				triangulator.generateSubPolygons(&coplanarVertices[0], &boundaries[0], coplanarVertices.size(), boundaries.size(), above, true, faces);
				for (auto face : faces) 
					surfaceRep.addTriangle(face.v[0]->getPointedVertex(), face.v[1]->getPointedVertex(), face.v[2]->getPointedVertex(), faceIndex);
				abovePoints.push_back(above);
			}
			else {
				surfaceRep.addTriangle(f.v[0]->getPointedVertex(), f.v[1]->getPointedVertex(), f.v[2]->getPointedVertex(), faceIndex);
				abovePoints.push_back(Geometer::calculateAbovePoint(f.v[0]->point, f.v[1]->point, f.v[2]->point));
			}

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
		if (predicator.orientCoplane(a->point, b->point, vertices[i]->point) != DelReal(0)) {
			c = vertices[i];
			std::swap(vertices[2], vertices[i]);
			break;
		}
	}
	Assert(c != NULL);
	for (std::vector<Vertex *>::size_type i = 3; i < vertices.size(); i++){
		DelReal ori = predicator.orient3d(a->point, b->point, c->point, vertices[i]->point);
		if (ori != DelReal(0)){
			d = vertices[i];
			std::swap(vertices[3], vertices[i]);
			if (ori < DelReal(0)) std::swap(a, b);
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

	VolumeVertexInsertionFlags vifs;
	Tetrahedron rt(a, b, c, d);
	for (std::vector<Vertex *>::size_type i = 4; i < vertices.size(); i++) {
		Vertex *v = vertices[i];
		Tetrahedron t = findPosition(v, rt, meshRep);
		insertVertex(v, t, meshRep, vifs, &rt);
	}
}

void DelMesher::constrainedTriangulation() {
	for (int i = 0; i < oriVertices.size(); i++) oriVertices[i]->setVertexPointer((Vertex *)i);
	Geometer::brioSort3d(oriVertices.begin(), oriVertices.end(), [](Vertex *v)->DelVector { return v->point; });
	triangulation3D(oriVertices, meshRep);
	for (int i = 0; i < oriVertices.size(); i++) {
		int oriIndex = (int)oriVertices[i]->getPointedVertex();
		while (oriIndex != i) {
			std::swap(oriVertices[i], oriVertices[oriIndex]);
			oriIndex = (int)oriVertices[i]->getPointedVertex();
		}
		oriVertices[i]->setVertexPointer(NULL);
	}

	for (auto iter = meshRep.segmentBegin(); iter != meshRep.segmentEnd(); ++iter) mayMissingSegs.push_back(*iter);
	std::sort(mayMissingSegs.begin(), mayMissingSegs.end(), 
		[](const Segment& left, const Segment& right) {
		return (left.v[0]->point - left.v[1]->point).length2() < (right.v[0]->point - right.v[1]->point).length2();
	});

	//recover segments
	while (!mayMissingSegs.empty()) {
		Segment s = mayMissingSegs.front();
		mayMissingSegs.pop_front();
		if (meshRep.isSegment(s) && !meshRep.Contain(s)) {
			Tetrahedron t;
			if (!meshRep.fastVertexQueryCheck(s.v[0])) std::swap(s.v[0], s.v[1]);
			if (findIntersectedTetrahedron(s.v[0], s.v[1]->point, &t)) {
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
				refineRegion(f, true, false);
				//recover missing segs in the process of region refinment
				while (!mayMissingSegs.empty()) {
					Segment s = mayMissingSegs.front();
					mayMissingSegs.pop_front();
					if (meshRep.isSegment(s) && !meshRep.Contain(s)) {
						Tetrahedron t;
						if (!meshRep.fastVertexQueryCheck(s.v[0])) std::swap(s.v[0], s.v[1]);
						if (findIntersectedTetrahedron(s.v[0], s.v[1]->point, &t)) {
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

DelVector DelMesher::calculateAbovePoint(const std::vector<Vertex *> &vertices, const Triangle& ref) {
	std::vector<Vertex *>::size_type vertexCount = vertices.size();
	const Vertex *a = vertices[0];
	const Vertex *b = NULL, *c = NULL;

	DelReal longest = 0;
	for (std::vector<Vertex *>::size_type i = 1; i < vertexCount; i++) {
		DelReal length = (vertices[i]->point - a->point).length2();
		if (length > longest) {
			b = vertices[i];
			longest = length;
		}
	}

	DelVector vec = b->point - a->point;
	DelReal largest = 0;
	for (std::vector<Vertex *>::size_type i = 2; i < vertexCount; i++) {
		DelReal area = (vec % (vertices[i]->point - a->point)).length2();
		if (area > largest) {
			c = vertices[i];
			largest = area;
		}
	}

	DelVector n = DelReal(0.5) * sqrt(longest) * Normalize(vec % (c->point - a->point));
	DelVector above = a->point + n;
	if (predicator.orient3d(above, ref.v[0]->point, ref.v[1]->point, ref.v[2]->point) < 0)
		above = a->point - n;

	return above;
}

void DelMesher::detectAcuteVertices() const {
	constexpr int maxCount = 5;
	constexpr DelReal acuteRadian = DelReal(89) / DelReal(180) * DelReal(M_PI);
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
						if (Geometer::interiorAngle(a->point, b->point, v->point) < acuteRadian) {
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
			if (predicator.inHalfSpace3D(u->point, b->point, c->point, d->point))
				return t;
			meshRep.Adjacent(Triangle(b, d, c, true), &x);
			a = x; std::swap(d, c);
		}
		else{
			bool ori0 = predicator.orient3d(u->point, a->point, b->point, c->point) > 0;
			bool ori1 = predicator.orient3d(u->point, a->point, c->point, d->point) > 0;
			bool ori2 = predicator.orient3d(u->point, a->point, d->point, b->point) > 0;

			int condition = ori0 + (ori1 << 1) + (ori2 << 2);

			switch (condition){
			case 0:
				if (predicator.orient3d(u->point, d->point, c->point, b->point) > 0){
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
			if (predicator.inHalfSpace3D(u->point, b->point, c->point, d->point)) {
				if (occluderRep.Contain(Triangle(b, c, d))) return Tetrahedron(NULL, b, c, d);
				else if (occluderRep.Contain(Triangle(b, d, c))) return Tetrahedron(NULL, b, d, c);
				else return t;
			}
			meshRep.Adjacent(Triangle(b, d, c, true), &x);
			a = x; std::swap(d, c);
		}
		else {
			bool ori0 = predicator.orient3d(u->point, a->point, b->point, c->point) > 0;
			bool ori1 = predicator.orient3d(u->point, a->point, c->point, d->point) > 0;
			bool ori2 = predicator.orient3d(u->point, a->point, d->point, b->point) > 0;

			int condition = ori0 + (ori1 << 1) + (ori2 << 2);

			switch (condition) {
			case 0:
				if (predicator.orient3d(u->point, d->point, c->point, b->point) > 0) {
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

	DelVector uProjected = Geometer::getProjectedPoint(a->point, b->point, c->point, u->point);
	auto plane = predicator.getProjectionPlane(a->point, b->point, c->point);
	bool faceOrient = predicator.orientCoplane(a->point, b->point, c->point, plane) > 0;

	auto orientTest = [plane, faceOrient](const DelVector& a, const DelVector& b, const DelVector& c) {
		constexpr Predicator<DelReal> predicator;
		DelReal ori = predicator.orientCoplane(a, b, c, plane);
		return faceOrient ? ori > 0 : ori < 0;
	};

	while (maxIter--){
		Vertex *x = NULL;
		bool ori0 = orientTest(b->point, a->point, uProjected);
		bool ori1 = orientTest(a->point, c->point, uProjected);

		int condition = ori0 + (ori1 << 1);
		switch (condition){
		case 0:
			if (orientTest(c->point, b->point, uProjected)){
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

	DelVector uProjected = Geometer::getProjectedPoint(a->point, b->point, c->point, u->point);
	auto plane = predicator.getProjectionPlane(a->point, b->point, c->point);
	bool facetOrient = predicator.orientCoplane(a->point, b->point, c->point, plane) > 0;
	
	auto orientTest = [plane, facetOrient](const DelVector& a, const DelVector& b, const DelVector& c) {
		constexpr Predicator<DelReal> predicator;
		DelReal ori = predicator.orientCoplane(a, b, c, plane);
		return facetOrient ? ori > 0 : ori < 0;
	};

	while (maxIter--) {
		Vertex *x = NULL;
		bool ori0 = orientTest(b->point, a->point, uProjected);
		bool ori1 = orientTest(a->point, c->point, uProjected);

		int condition = ori0 + (ori1 << 1);
		switch (condition) {
		case 0:
			if (orientTest(c->point, b->point, uProjected)) {
				if (meshRep.isSegment(Segment(c, b))) return Triangle(NULL, b, c);
				surfaceRep.Adjacent(Segment(c, b), &x);
				a = x; std::swap(b, c); //Triangle(x, c, b)
			}
			else {
				return Triangle(a, b, c, true);
			}
			break;
		case 1:
			if (meshRep.isSegment(Segment(b, a))) return Triangle(NULL, a, b);
			surfaceRep.Adjacent(Segment(b, a), &x);
			c = a; a = x; //Triangle(x, b, a)
			break;
		case 2:
			if (meshRep.isSegment(Segment(a, c))) return Triangle(NULL, c, a);
			surfaceRep.Adjacent(Segment(a, c), &x);
			b = a; a = x; //Triangle(x, a, c)
			break;
		case 3:
			if (Randomnation<2>() == 0) {
				if (meshRep.isSegment(Segment(b, a))) return Triangle(NULL, a, b);
				surfaceRep.Adjacent(Segment(b, a), &x);
				c = a; a = x; //Triangle(x, b, a)
			}
			else {
				if (meshRep.isSegment(Segment(a, c))) return Triangle(NULL, c, a);
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
	const VolumeVertexInsertionFlags& vifs, int depth, Tetrahedron *rt){
	if (depth < 200) {
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
				deleted = predicator.inHalfSpace3D(u->point, b->point, c->point, d->point);
			else if (b->isGhost()) 
				deleted = predicator.inHalfSpace3D(u->point, a->point, d->point, c->point);
			else 
				deleted = predicator.inSpherePerturbed(u->point, b->point, d->point, c->point, a->point) > 0;
		}

		if (deleted) {
			if (!vifs.trueInsertion) tobeDeletedTets.push_back(Tetrahedron(a, b, c, d));
			meshRep.deleteTetrahedron(a, b, c, d);

			if (vifs.missingFacetTest || vifs.refineFacetTest) {
				Triangle front(b, c, d), back(b, d, c);
				if (surfaceRep.Contain(front)) {
					if (vifs.missingFacetTest) mayMissingFacets.push_back(front);
					if (vifs.refineFacetTest) mayRefineFacets.push_back(front);
				}
				else if (surfaceRep.Contain(back)) {
					if (vifs.missingFacetTest) mayMissingFacets.push_back(back);
					if (vifs.refineFacetTest) mayRefineFacets.push_back(back);
				}
			}
			if (vifs.refineSegTest || vifs.missingSegTest) {
				Segment ab(a, b), ac(a, c), ad(a, d);
				if (meshRep.isSegment(ab) && !meshRep.testAndMark(a, b)) markedSegments.push_back(ab);
				if (meshRep.isSegment(ac) && !meshRep.testAndMark(a, c)) markedSegments.push_back(ac);
				if (meshRep.isSegment(ad) && !meshRep.testAndMark(a, d)) markedSegments.push_back(ad);
			}

			digCavity(u, Triangle(a, c, d, true), meshRep, vifs, depth + 1, rt);
			digCavity(u, Triangle(a, b, c, true), meshRep, vifs, depth + 1, rt);
			digCavity(u, Triangle(a, d, b, true), meshRep, vifs, depth + 1, rt);
		}
		else {
			if (vifs.trueInsertion && !vifs.cdt) meshRep.addTetrahedron(u, b, d, c);
			else newFacesOfTets.push_back(Triangle(b, d, c));

			if (vifs.cdt && !b->isGhost())
				vifs.rejected |= (predicator.orient3d(u->point, b->point, d->point, c->point) <= 0);

			if (rt) *rt = Tetrahedron(u, b, d, c);
			if (vifs.refineFacetTest) {
				Triangle front(b, d, c), back(b, c, d);
				if (surfaceRep.Contain(front)) mayRefineFacets.push_back(front);
				if (surfaceRep.Contain(back)) mayRefineFacets.push_back(back);
			}
			if (!vifs.cdt && vifs.skinnyTetTest) 
				maySkinnyTets.push_back(Tetrahedron(u, b, d, c));

			if (vifs.insertRadiusTest) {
				DelReal minLen2 = (vifs.parent != NULL) ? (vifs.parent->point - u->point).length2()
					: std::numeric_limits<DelReal>::max();
				for (int i = 0; i < 3; i++) {
					if (!f.v[i]->isGhost()) {
						DelReal len2 = (f.v[i]->point - u->point).length2();
						if (minLen2 > len2) {
							vifs.parent = f.v[i];
							minLen2 = len2;
						}
					}
				}
			}
		}
	}
	else {
		std::stack<Triangle> workingStacks;
		workingStacks.push(f);

		while (!workingStacks.empty()) {
			Triangle top = workingStacks.top();
			workingStacks.pop();
			Vertex *a = NULL;
			if (!meshRep.Adjacent(top, &a)) {
				if (vifs.cdt)
					vifs.rejected |= (surfaceRep.Contain(top) || surfaceRep.Contain(Triangle(top.v[0], top.v[2], top.v[1])));
				continue;
			}

			bool success = true;

			Vertex *b = top.v[0];
			Vertex *c = top.v[1];
			Vertex *d = top.v[2];

			bool deleted = false;
			bool constrained = false;
			if (vifs.cdt) constrained = surfaceRep.Contain(top) || surfaceRep.Contain(Triangle(b, d, c));
			if (!constrained) {
				if (a->isGhost())
					deleted = predicator.inHalfSpace3D(u->point, b->point, c->point, d->point);
				else if (b->isGhost())
					deleted = predicator.inHalfSpace3D(u->point, a->point, d->point, c->point);
				else
					deleted = predicator.inSpherePerturbed(u->point, b->point, d->point, c->point, a->point) > 0;
			}

			if (deleted) {
				if (!vifs.trueInsertion) tobeDeletedTets.push_back(Tetrahedron(a, b, c, d));
				meshRep.deleteTetrahedron(a, b, c, d);

				if (vifs.missingFacetTest || vifs.refineFacetTest) {
					Triangle front(b, c, d), back(b, d, c);
					if (surfaceRep.Contain(front)) {
						if (vifs.missingFacetTest) mayMissingFacets.push_back(front);
						if (vifs.refineFacetTest) mayRefineFacets.push_back(front);
					}
					else if (surfaceRep.Contain(back)) {
						if (vifs.missingFacetTest) mayMissingFacets.push_back(back);
						if (vifs.refineFacetTest) mayRefineFacets.push_back(back);
					}
				}
				if (vifs.refineSegTest || vifs.missingSegTest) {
					Segment ab(a, b), ac(a, c), ad(a, d);
					if (meshRep.isSegment(ab) && !meshRep.testAndMark(a, b)) markedSegments.push_back(ab);
					if (meshRep.isSegment(ac) && !meshRep.testAndMark(a, c)) markedSegments.push_back(ac);
					if (meshRep.isSegment(ad) && !meshRep.testAndMark(a, d)) markedSegments.push_back(ad);
				}

				workingStacks.push(Triangle(a, c, d, true));
				workingStacks.push(Triangle(a, b, c, true));
				workingStacks.push(Triangle(a, d, b, true));
			}
			else {
				if (vifs.trueInsertion && !vifs.cdt) meshRep.addTetrahedron(u, b, d, c);
				else newFacesOfTets.push_back(Triangle(b, d, c));

				if (vifs.cdt && !b->isGhost())
					vifs.rejected |= (predicator.orient3d(u->point, b->point, d->point, c->point) <= 0);

				if (rt) *rt = Tetrahedron(u, b, d, c);
				if (vifs.refineFacetTest) {
					Triangle front(b, d, c), back(b, c, d);
					if (surfaceRep.Contain(front)) mayRefineFacets.push_back(front);
					if (surfaceRep.Contain(back)) mayRefineFacets.push_back(back);
				}
				if (!vifs.cdt && vifs.skinnyTetTest)
					maySkinnyTets.push_back(Tetrahedron(u, b, d, c));

				if (vifs.insertRadiusTest) {
					DelReal minLen2 = (vifs.parent != NULL) ? (vifs.parent->point - u->point).length2()
						: std::numeric_limits<DelReal>::max();
					for (int i = 0; i < 3; i++) {
						if (!top.v[i]->isGhost()) {
							DelReal len2 = (top.v[i]->point - u->point).length2();
							if (minLen2 > len2) {
								vifs.parent = top.v[i];
								minLen2 = len2;
							}
						}
					}
				}
			}
		}
	}
}

void DelMesher::digCavity(Vertex *u, const Segment &s, const DelVector& above, int index, const SurfaceVertexInsertionFlags& vifs){
	Vertex *w = s.v[0], *v = s.v[1], *x = NULL;
	int newIndex = -1;

	bool notBoundary = !meshRep.isSegment(s);
	if (!surfaceRep.Adjacent(s, &x, &newIndex)) {
		if (notBoundary) return;
	}

	if (notBoundary && predicator.inCirclePerturbed(u->point, v->point, w->point, x->point, above) > 0){
		if (vifs.deletionRecord)
			tobeDeletedFaces.push_back(TriangleWithIndex(w, v, x, index, true));
		surfaceRep.deleteTriangle(w, v, x);
		Assert(index == newIndex);
		digCavity(u, Segment(x, v), above, newIndex, vifs);
		digCavity(u, Segment(w, x), above, newIndex, vifs);
	}
	else {
		if (vifs.trueInsertion)
			surfaceRep.addTriangle(u, v, w, index);
		else 
			newSegsOfFaces.push_back(SegmentWithIndex(v, w, index));
		if (vifs.missingTest) mayMissingFacets.push_back(Triangle(u, v, w));
		if (vifs.refineTest) mayRefineFacets.push_back(Triangle(u, v, w));
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

	const Predicator<DelReal> predicator;
	enum Move { Left, Right, Horizon };
	bool terminate = false;
	do {
		bool ori0 = predicator.orient3d(bb, a->point, d->point, c->point) >= 0;
		bool ori1 = predicator.orient3d(bb, a->point, b->point, d->point) >= 0;
		bool ori2 = predicator.orient3d(bb, a->point, c->point, b->point) >= 0;

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

	DelVector startVert = a->point, endVert = end->point;
	Vertex *ref = b;
	DelReal maxAngle = Geometer::interiorAngle(ref->point, startVert, endVert);

	DelReal angle = Geometer::interiorAngle(c->point, startVert, endVert);
	if (angle > maxAngle) {
		ref = c;
		maxAngle = angle;
	}
	angle = Geometer::interiorAngle(d->point, startVert, endVert);
	if (angle > maxAngle) {
		ref = d;
		maxAngle = angle;
	}

	enum Intersection { AcrossFace, AcrossEdge };
	Intersection intersect = Intersection::AcrossFace;
	if (predicator.orient3d(endVert, a->point, c->point, b->point) == 0)
		intersect = Intersection::AcrossEdge;
	else if (predicator.orient3d(endVert, a->point, d->point, c->point) == 0) {
		std::swap(d, b);
		std::swap(b, c);
		intersect = Intersection::AcrossEdge;
	}
	else if (predicator.orient3d(endVert, a->point, b->point, d->point) == 0) {
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
				angle = Geometer::interiorAngle(a->point, startVert, endVert);
				if (angle > maxAngle) {
					maxAngle = angle;
					ref = a;
				}
				std::swap(a, d);
				bool dbc = predicator.Intersection(d->point, b->point, c->point, startVert, endVert);
				bool dab = predicator.Intersection(d->point, a->point, b->point, startVert, endVert);
				bool dca = predicator.Intersection(d->point, c->point, a->point, startVert, endVert);
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
					angle = Geometer::interiorAngle(a->point, startVert, endVert);
					if (angle > maxAngle) {
						maxAngle = angle;
						ref = a;
					}
					std::swap(a, d);
					bool dab = predicator.Intersection(d->point, a->point, b->point, startVert, endVert);
					bool dca = predicator.Intersection(d->point, c->point, a->point, startVert, endVert);
					bool coplane = predicator.orient3d(endVert, b->point, c->point, d->point) == DelReal(0);
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
	DelVector aa = a->point, bb = b->point;
	DelVector ab = bb - aa;
	DelReal abLen = ab.length();
	DelReal invabLen = DelReal(1) / abLen;

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
		DelVector acutePoint = acuteVert->point;
		DelVector otherPoint = otherVert->point;
		DelVector refPoint = ref->point;
		DelReal r = (refPoint - acutePoint).length();
		steinerPoint = acutePoint + (r * invabLen) * (otherPoint - acutePoint);

		DelReal steOtherLen = (steinerPoint - otherPoint).length();
		DelReal steRefLen = (steinerPoint - refPoint).length();

		if (steOtherLen < steRefLen) {
			r = nearestPowerOfTwo(DelReal(0.5) * abLen);
			steinerPoint = acutePoint + (r * invabLen) * (otherPoint - acutePoint);
		}
	}
	else {
		DelVector ar = ref->point - aa;
		DelVector br = ref->point - bb;
		DelReal arLen = ar.length();
		DelReal brLen = br.length();
		if (arLen < DelReal(0.5) * abLen)
			steinerPoint = aa + (arLen * invabLen) * ab;
		else if (brLen < DelReal(0.5) * abLen)
			steinerPoint = bb - (brLen * invabLen) * ab;
		else
			steinerPoint = DelReal(0.5) * (aa + bb);
	}

	//check to see steiner point isn't too close to both end point
	DelReal t = (steinerPoint - aa).length() * invabLen;
	if (t < DelReal(0.2) || t > DelReal(0.8)) {
		if (!acuteSplit) 
			steinerPoint = DelReal(0.5) * (aa + bb);
		else {
			DelReal r = nearestPowerOfTwo(DelReal(0.5) * abLen);
			steinerPoint = acuteVert->point + (r * invabLen) * (otherVert->point - acuteVert->point);
		}
	}

	VertexType type = VertexType::Vertex_FreeSegment;
	bool onSurface = surfaceRep.Contain(s) || surfaceRep.Contain(Segment(s.v[1], s.v[0]));
	if (onSurface) type = VertexType(type | VertexType::Vertex_Facet);
	Vertex *steinerVert = allocVertex(steinerPoint, DelReal(0), type);

	mayMissingSegs.push_front(Segment(a, steinerVert));
	mayMissingSegs.push_front(Segment(steinerVert, b));

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
	for (auto iter = meshRep.segmentBegin(); iter != meshRep.segmentEnd(); ++iter) mayRefineSegs.push_back(*iter);
refinmentStart:
	//phase 1: refine subsegments that do not meet criteria
	while (!mayRefineSegs.empty()) {
		Segment s = mayRefineSegs.front();
		mayRefineSegs.pop_front();
		if (meshRep.isSegment(s) && (Encroached(s) || s.getLength() > criteria.maxSegLength))
			refineSubSegment(s, NULL, false, encorachedFacetTest, skinnyTest);
	}
	//recover missing subpolygons
	while (!mayMissingFacets.empty()) {
		Triangle f = mayMissingFacets.front();
		if (surfaceRep.Contain(f) && !meshRep.Contain(f)) {
			if (!facetRecovery(f, regionVertices, regionBoundaries, regionFaces,
				positiveVertices, positiveFaces, negativeVertices, negativeFaces,
				deleted, inserted, encorachedFacetTest, skinnyTest)) {
				refineRegion(f, false, true, true, skinnyTest);
				goto refinmentStart;
			}
		}
		mayMissingFacets.pop_front();
	}

	if (!facetsCollected) { //collect subpolygons for the first round
		for (auto f : surfaceRep) mayRefineFacets.push_back(f);
		facetsCollected = true;
	}
	encorachedFacetTest = true;

	//phase 2: refine subpolygons that do not meet criteria
	while (!mayRefineFacets.empty()) {
		TriangleWithGeometry f = mayRefineFacets.front();
		if (surfaceRep.Contain(f)) {
			Vertex *encroached = NULL;
			if (Encroached(f, &encroached)) {
				TriangleWithGeometry splited = findPositionWithOcclusion(encroached, f);
				if (splited.v[0]) refineSubPolygon(splited, skinnyTest);
				else {
					Segment s(splited.v[1], splited.v[2]);
					bool entwinement = false;
					if (matchVertexFlag(encroached->getVertexType(), VertexType::Vertex_Refined) &&
						!matchVertexFlag(encroached->getVertexType(), VertexType::Vertex_Segment))
						entwinement = Adjacent(s, encroached);
					else encroached = NULL;
					refineSubSegment(s, encroached, entwinement, encorachedFacetTest, skinnyTest);
				}
			}
			else if (f.getRadius() > criteria.maxFacetRadius) {
				refineSubPolygon(f, skinnyTest);
				mayRefineFacets.pop_front();
			}
			else mayRefineFacets.pop_front();

			if (!mayRefineSegs.empty() || !mayMissingFacets.empty()) goto refinmentStart;
		}
		else mayRefineFacets.pop_front();
	}

	if (!tetsCollected) { // collect skinny tets for the first round
		for (auto t : meshRep) {
			if (skinnyTetTest(t)) maySkinnyTets.push_back(t);
		}
		tetsCollected = true;
	}
	skinnyTest = true;
	
	// phase 3: refine tetrahedrons that do not meet criteria
	while (!maySkinnyTets.empty()) {
		Tetrahedron t = maySkinnyTets.front();
		if (meshRep.Contain(t) && skinnyTetTest(t)) {
			bool processed = false;
			refineTetrahedron(t, processed);
			if (processed) maySkinnyTets.pop_front();
			if (!mayRefineSegs.empty() || !mayRefineFacets.empty() || !mayMissingFacets.empty()) 
				goto refinmentStart;
		}
		else maySkinnyTets.pop_front();
	}
}

void DelMesher::refineSubSegment(const Segment &s, Vertex *parent, bool entwinement, bool refineFacetTestEnable, bool skinnyTetTestEnable){
	Vertex *a = s.v[0], *b = s.v[1];
	DelVector ab = b->point - a->point;
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
		steinerPoint = a->point + (ab * DelReal(0.5));
	else {
		//calulate steiner point for acute vertex
		DelReal abLen = ab.length();
		DelReal r = nearestPowerOfTwo(DelReal(0.5) * abLen);
		if (!bothAcute) {
			if (aIsAcute) steinerPoint = a->point + (r / abLen) * ab;
			else steinerPoint = b->point - (r / abLen) * ab;
		}
		else {
			if (r > abLen * DelReal(0.5)) r *= DelReal(0.5);
			steinerPoint = a->point + (r / abLen) * ab;

			DelReal r2 = nearestPowerOfTwo(DelReal(0.5) * (steinerPoint - b->point).length());
			steinerPointOther = b->point - (r2 / abLen) * ab;
		}
	}

	VertexType type = VertexType::Vertex_RefinedFreeSegment;
	if (surfaceRep.Contain(s) || surfaceRep.Contain(Segment(b, a))) 
		type = VertexType(type | VertexType::Vertex_Facet);
	Vertex *steinerVert = allocVertex(steinerPoint, DelReal(0), type);

	Tetrahedron toBeDeleted;
	if (!meshRep.adjacent2SegmentFast(s, &toBeDeleted)) {
		if (meshRep.fastVertexQueryCheck(a)) meshRep.adjacent2Vertex(a, &toBeDeleted);
		else meshRep.adjacent2Vertex(b, &toBeDeleted);
	}
	toBeDeleted = findPosition(steinerVert, toBeDeleted, meshRep);

	mayRefineSegs.push_back(Segment(a, steinerVert));
	if (!bothAcute) mayRefineSegs.push_back(Segment(b, steinerVert));

	if (matchVertexFlag(a->getVertexType(), VertexType::Vertex_FreeSegment))
		steinerVert->setOriSegmentIndex(a->getOriSegmentIndex());
	else if (matchVertexFlag(b->getVertexType(), VertexType::Vertex_FreeSegment))
		steinerVert->setOriSegmentIndex(b->getOriSegmentIndex());
	else {
		steinerVert->setOriSegmentIndex(oriSegments.size());
		oriSegments.push_back(s);
	}

	VolumeVertexInsertionFlags vifs;
	vifs.cdt = true;
	vifs.missingFacetTest = true;
	vifs.refineSegTest = true; vifs.refineFacetTest = refineFacetTestEnable;
	vifs.skinnyTetTest = skinnyTetTestEnable;
	vifs.insertRadiusTest = (parent == NULL) || (entwinement && midSplit);
	insertVertexOnSegment(steinerVert, s, toBeDeleted, meshRep, vifs);

	if (vifs.insertRadiusTest) {
		steinerVert->relaxedInsetionRadius = (vifs.parent->point - steinerPoint).length();
		if (midSplit) {
			if (Adjacent(s, vifs.parent)) {
				DelReal parentRIS = vifs.parent->relaxedInsetionRadius;
				if (!matchVertexFlag(vifs.parent->getVertexType(), VertexType::Vertex_Segment) &&
					matchVertexFlag(vifs.parent->getVertexType(), VertexType::Vertex_Refined)) {
					parentRIS /= DelReal(SQRT_2);
					if (parent == NULL) parentRIS = std::max(ab.length() * DelReal(0.5), parentRIS);
				}
				steinerVert->relaxedInsetionRadius = parentRIS;
			}
			if (vifs.parent == a || vifs.parent == b)
				steinerVert->relaxedInsetionRadius = std::max(steinerVert->relaxedInsetionRadius,
					std::min(a->relaxedInsetionRadius, b->relaxedInsetionRadius));
		}
		else {
			if (Adjacent(s, vifs.parent)) {
				steinerVert->relaxedInsetionRadius = std::max(steinerVert->relaxedInsetionRadius, 
					vifs.parent->relaxedInsetionRadius);
				if (!matchVertexFlag(vifs.parent->getVertexType(), VertexType::Vertex_Segment) &&
					matchVertexFlag(vifs.parent->getVertexType(), VertexType::Vertex_Refined))
					steinerVert->relaxedInsetionRadius = std::min((steinerPoint - b->point).length(),
					(steinerPoint - a->point).length());
			}
		}
		if (!matchVertexFlag(vifs.parent->getVertexType(), VertexType::Vertex_LowDimension))
			steinerVert->relaxedInsetionRadius = std::max(steinerVert->relaxedInsetionRadius,
				vifs.parent->relaxedInsetionRadius / DelReal(SQRT_2));

		if (entwinement) {
			steinerVert->relaxedInsetionRadius = std::min(parent->relaxedInsetionRadius / DelReal(SQRT_2),
				steinerVert->relaxedInsetionRadius);
			DelReal midLen = DelReal(0.5) * ab.length();
			steinerVert->relaxedInsetionRadius = std::max(midLen, steinerVert->relaxedInsetionRadius);
		}
	}
	else
		steinerVert->relaxedInsetionRadius = std::min((steinerPoint - a->point).length(), (steinerPoint - b->point).length());

	if (bothAcute) {
		Segment ns(b, steinerVert);
		Vertex *steinerVertOther = allocVertex(steinerPointOther, DelReal(0), type);

		Tetrahedron toBeDeleted2;
		if (!meshRep.adjacent2SegmentFast(ns, &toBeDeleted2)) {
			if (meshRep.fastVertexQueryCheck(b)) meshRep.adjacent2Vertex(b, &toBeDeleted2);
			else meshRep.adjacent2Vertex(steinerVert, &toBeDeleted2);
		}
		toBeDeleted2 = findPosition(steinerVertOther, toBeDeleted2, meshRep);

		mayRefineSegs.push_back(Segment(steinerVert, steinerVertOther));
		mayRefineSegs.push_back(Segment(b, steinerVertOther));

		steinerVertOther->setOriSegmentIndex(steinerVert->getOriSegmentIndex());

		vifs.parent = NULL; vifs.rejected = false;
		insertVertexOnSegment(steinerVertOther, ns, toBeDeleted2, meshRep, vifs);

		if (vifs.insertRadiusTest) {
			steinerVertOther->relaxedInsetionRadius = (vifs.parent->point - steinerPointOther).length();
			if (Adjacent(ns, vifs.parent)) {
				steinerVertOther->relaxedInsetionRadius = std::max(steinerVertOther->relaxedInsetionRadius,
					vifs.parent->relaxedInsetionRadius);
				if (!matchVertexFlag(vifs.parent->getVertexType(), VertexType::Vertex_Segment) &&
					matchVertexFlag(vifs.parent->getVertexType(), VertexType::Vertex_Refined))
					steinerVertOther->relaxedInsetionRadius = std::min((steinerPointOther - b->point).length(), 
					(steinerPointOther - steinerPoint).length());
			}
			if (!matchVertexFlag(vifs.parent->getVertexType(), VertexType::Vertex_LowDimension))
				steinerVert->relaxedInsetionRadius = std::max(steinerVert->relaxedInsetionRadius,
					vifs.parent->relaxedInsetionRadius / DelReal(SQRT_2));
		}
		else
			steinerVertOther->relaxedInsetionRadius = std::min((steinerPointOther - b->point).length(), 
			(steinerPointOther - steinerPoint).length());
	}
}

void DelMesher::refineSubPolygon(TriangleWithGeometry &f, bool skinnyTetTestEnable) {
	f.setGeometricProperties();
	Vertex *a = f.v[0], *b = f.v[1], *c = f.v[2], *d = NULL;
	int facetIndex = surfaceRep.getTriangleIndex(a, b, c);
	DelVector centerPoint = f.getCircumcenter();
	DelReal radius2 = f.getRadius() * f.getRadius();
	Assert(radius2 > 0);

	Vertex *steinerVert = allocVertex(centerPoint, DelReal(0), VertexType::Vertex_RefinedFreeFacet);

	VolumeVertexInsertionFlags vifs;
	vifs.cdt = true; vifs.refineSegTest = true; vifs.refineFacetTest = true;
	vifs.trueInsertion = false; vifs.insertRadiusTest = true;
	meshRep.Adjacent(f, &d); Assert(d);
	insertVertexOnSurface(steinerVert, f, findPosition(steinerVert, Tetrahedron(d, a, b, c), meshRep), meshRep, vifs);

	//set relaxed insertion radius
	steinerVert->relaxedInsetionRadius = (centerPoint - vifs.parent->point).length();
	if (Adjacent(facetIndex, vifs.parent)) {
		DelReal parentRIS = vifs.parent->relaxedInsetionRadius;
		if (matchVertexFlag(vifs.parent->getVertexType(), VertexType::Vertex_Segment))
			parentRIS *= DelReal(SQRT_2);
		steinerVert->relaxedInsetionRadius = std::max(steinerVert->relaxedInsetionRadius,
			parentRIS);
	}

	//detect encorachment
	bool encorach = false;
	Segment encorachedSeg;
	for (auto s : mayRefineSegs) {
		if (Encroached(s, steinerVert)) {
			encorach = true;
			encorachedSeg = s;
			break;
		}
	}
	mayRefineSegs.clear();

	if (!encorach) {
		std::vector<Vertex *> encorachedVolumeVerts;
		for (auto face : newFacesOfTets) {
			for (int i = 0; i < 3; i++) {
				if (!matchOneOfVertexFlags(face.v[i]->getVertexType(), VertexType::Vertex_LowDimension) &&
					Encroached(centerPoint, radius2, face.v[i]) && !face.v[i]->isMarked()) {
					face.v[i]->setMark();
					encorachedVolumeVerts.push_back(face.v[i]);
				}
			}
		}

		if (encorachedVolumeVerts.empty()) {
			if (!vifs.rejected) {
				for (auto face : newFacesOfTets) {
					meshRep.addTetrahedron(steinerVert, face.v[0], face.v[1], face.v[2]);
					if (skinnyTetTestEnable) {
						Tetrahedron newTet(Tetrahedron(steinerVert, face.v[0], face.v[1], face.v[2]));
						if (skinnyTetTest(newTet)) maySkinnyTets.push_back(newTet);
					}
				}
			}
			else { //SOS routine
				vifs.trueInsertion = true; vifs.missingFacetTest = true;
				vifs.skinnyTetTest = skinnyTetTestEnable;
				triangulateCavity(steinerVert, newFacesOfTets, meshRep, vifs);
			}

			surfaceRep.bindVolumeVertex(steinerVert);
			for (auto s : newSegsOfFaces) {
				if (!meshRep.Contain(Triangle(steinerVert, s.v[0], s.v[1])))
					mayMissingFacets.push_back(Triangle(steinerVert, s.v[0], s.v[1]));
				surfaceRep.addTriangle(steinerVert, s.v[0], s.v[1], s.index);
			}
			tobeDeletedFaces.clear(); newSegsOfFaces.clear();
			tobeDeletedTets.clear(); newFacesOfTets.clear();
		}
		else {
			//restore
			for (auto face : tobeDeletedFaces) surfaceRep.addTriangle(face.v[0], face.v[1], face.v[2], face.index);
			for (auto t : tobeDeletedTets) meshRep.addTetrahedron(t.v[0], t.v[1], t.v[2], t.v[3]);
			tobeDeletedFaces.clear(); newSegsOfFaces.clear();
			tobeDeletedTets.clear(); newFacesOfTets.clear();

			VolumeVertexInsertionFlags dvifs;
			dvifs.skinnyTetTest = skinnyTetTestEnable;
			for (auto v : encorachedVolumeVerts) {
				v->unSetMark();
				deleteVertex(v, meshRep, dvifs);
			}
			deallocVertex(steinerVert);
			refineSubPolygon(f, skinnyTetTestEnable);
		}
	}
	else {
		//restore
		for (auto face : tobeDeletedFaces) surfaceRep.addTriangle(face.v[0], face.v[1], face.v[2], face.index);
		for (auto t : tobeDeletedTets) meshRep.addTetrahedron(t.v[0], t.v[1], t.v[2], t.v[3]);
		tobeDeletedFaces.clear(); newSegsOfFaces.clear();
		tobeDeletedTets.clear(); newFacesOfTets.clear();

		//split encorached subsegment
		refineSubSegment(encorachedSeg, steinerVert, Adjacent(encorachedSeg, facetIndex), true, skinnyTetTestEnable);
		deallocVertex(steinerVert);
	}
}

void DelMesher::refineTetrahedron(const Tetrahedron& tet, bool& processed) {
	Vertex *a = tet.v[0], *b = tet.v[1], *c = tet.v[2], *d = tet.v[3];
	DelVector center = tet.getCircumcenter();
	DelReal radius = tet.getRadius();

	Assert(radius >= 0);
	Vertex *steinerVert = meshRep.allocVertex(center, DelReal(0), VertexType(Vertex_Free | Vertex_Refined));
	Tetrahedron found = findPositionWithOcclusion(steinerVert, tet, meshRep, surfaceRep);

	if (found.v[0]) {
		//reject out of hull vertex
		if (found.v[0]->isGhost() || found.v[1]->isGhost() || found.v[2]->isGhost() || found.v[3]->isGhost()) {
			meshRep.deallocVertex(steinerVert);
			processed = true;
			return;
		}

		VolumeVertexInsertionFlags vifs;
		vifs.cdt = true; vifs.refineSegTest = true; vifs.refineFacetTest = true;
		vifs.trueInsertion = false;
		insertVertex(steinerVert, found, meshRep, vifs);

		//detect subsegment encorachment first
		bool encorach = false;
		Segment encorachedSeg;
		for (auto s : mayRefineSegs) {
			if (Encroached(s, steinerVert)) {
				encorach = true;
				encorachedSeg = s;
				break;
			}
		}
		mayRefineSegs.clear();

		if (encorach) {
			for (auto t : tobeDeletedTets) meshRep.addTetrahedron(t.v[0], t.v[1], t.v[2], t.v[3]);
			tobeDeletedTets.clear(); newFacesOfTets.clear();
			mayRefineFacets.clear();
			refineSubSegment(encorachedSeg, steinerVert, false, true, true);
			deallocVertex(steinerVert);
			processed = false;
		}
		else {
			bool direct = false; TriangleWithGeometry encorachedFace;

			auto directEncroachmentTest = [](const DelVector& u, const DelVector& a, const DelVector& b, const DelVector& c) {
				constexpr Predicator<DelReal> predicator;
				DelVector uProjected = Geometer::getProjectedPoint(u, a, b, c);
				auto plane = predicator.getProjectionPlane(a, b, c);
				bool ori0 = predicator.orientCoplane(uProjected, a, b, plane) >= 0;
				bool ori1 = predicator.orientCoplane(uProjected, b, c, plane) >= 0;
				bool ori2 = predicator.orientCoplane(uProjected, c, a, plane) >= 0;
				return ori0 == ori1 && ori0 == ori2;
			};

			for (auto f : mayRefineFacets) {
				if (Encroached(f, steinerVert)) {
					encorach = true;
					encorachedFace = f;
					//detect direct encorachment
					if (directEncroachmentTest(center, f.v[0]->point, f.v[1]->point, f.v[2]->point)) {
						direct = true;
						break;
					}
				}
			}
			mayRefineFacets.clear();
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

				for (auto f : newFacesOfTets) {
					meshRep.addTetrahedron(steinerVert, f.v[0], f.v[1], f.v[2]);
					Tetrahedron t(steinerVert, f.v[0], f.v[1], f.v[2]);
					if (skinnyTetTest(t)) maySkinnyTets.push_back(t);
				}

				tobeDeletedTets.clear(); newFacesOfTets.clear();
				processed = true;
			}
		}
	}
	else {
		TriangleWithGeometry occluder(found.v[1], found.v[2], found.v[3]);
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
				if (matchVertexFlag(encorachedVert->getVertexType(), VertexType::Vertex_Refined) && 
					!matchVertexFlag(encorachedVert->getVertexType(), VertexType::Vertex_Segment))
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
				if (predicator.Intersection(f.v[0]->point, f.v[1]->point, f.v[2]->point, d->point, c->point) &&
					predicator.orient3d(c->point, f.v[0]->point, f.v[1]->point, f.v[2]->point) > DelReal(0)) {
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
bool DelMesher::findCavity(const std::vector<Segment>& regionBoundaries, const std::vector<Triangle>& regionFaces,
	std::vector<Vertex *>& positiveVertices, std::vector<Triangle>& positiveFaces,
	std::vector<Vertex *>& negativeVertices, std::vector<Triangle>& negativeFaces,
	std::vector<Tetrahedron>& deleted, bool refineFacetTest) {
	enum CavityTriangleType { TwoPositive, TwoNegative, TwoCoplanar, PositiveNegativeCoplanar, NegativePositiveCoplanar };
	std::stack<std::pair<Triangle, CavityTriangleType>> triangleStack;
	std::vector<Segment> newSegments;
	std::vector<Triangle> newFacets;
	newFacets.reserve(regionFaces.size());

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
		newFacets.push_back(Triangle(seg.v[0], seg.v[1], onFace));

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
				if (predicator.Intersection(f.v[0]->point, f.v[1]->point, f.v[2]->point, a->point, b->point, d->point)) {
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
				if (predicator.Intersection(f.v[0]->point, f.v[1]->point, f.v[2]->point, a->point, b->point, d->point)) {
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
					if (predicator.Intersection(f.v[0]->point, f.v[1]->point, f.v[2]->point, c->point, d->point)) {
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
					newFacets.push_back(Triangle(b, a, d));
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
			DelVector pp = a->point, qq = b->point, rr = d->point;
			Triangle intersected;
			bool hasIntersection = false;
			for (auto f : regionFaces) {
				if (predicator.Intersection(f.v[0]->point, f.v[1]->point, f.v[2]->point, pp, rr) ||
					predicator.Intersection(f.v[0]->point, f.v[1]->point, f.v[2]->point, qq, rr)) {
					intersected = f;
					hasIntersection = true;
					break;
				}
			}
			Assert(hasIntersection);

			DelReal ori = predicator.orient3d(rr, intersected.v[0]->point, intersected.v[1]->point, intersected.v[2]->point);
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
			DelVector pp = a->point, qq = b->point, rr = d->point;
			Triangle intersected;
			bool hasIntersection = false;
			for (auto f : regionFaces) {
				if (predicator.Intersection(f.v[0]->point, f.v[1]->point, f.v[2]->point, pp, rr) ||
					predicator.Intersection(f.v[0]->point, f.v[1]->point, f.v[2]->point, qq, rr)) {
					intersected = f;
					hasIntersection = true;
					break;
				}
			}
			Assert(hasIntersection);

			DelReal ori = predicator.orient3d(rr, intersected.v[0]->point, intersected.v[1]->point, intersected.v[2]->point);
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
		if (newFacets.size() != regionFaces.size()) {
			DelTriangulator triangulator;
			for (auto f : newFacets) surfaceRep.addTriangle(f.v[0], f.v[1], f.v[2], faceIndex);
			newFacets.clear();
			triangulator.insertSegments(&regionFaces[0], &newSegments[0], regionFaces.size(), newSegments.size(), abovePoints[faceIndex], newFacets);
			for (auto f : newFacets) {
				if (!surfaceRep.Contain(f)) {
					surfaceRep.addTriangle(f.v[0], f.v[1], f.v[2], faceIndex);
					mayMissingFacets.push_back(f);
				}
			}
		}
		else {
			for (auto f : newFacets) surfaceRep.addTriangle(f.v[0], f.v[1], f.v[2], faceIndex);
		}

		if (refineFacetTest) {
			for (auto f : newFacets) mayRefineFacets.push_back(f);
		}

		for (auto v : positiveVertices) v->unSetMark();
		for (auto v : negativeVertices) v->unSetMark();

		return false;
	}

	for (auto v : positiveVertices) v->unSetMark();
	for (auto v : negativeVertices) v->unSetMark();

	return true;
}

bool DelMesher::facetRecovery(Triangle& f, std::vector<Vertex *>& regionVertices,
	std::vector<Segment>& regionBoundaries, std::vector<Triangle>& regionFaces,
	std::vector<Vertex *>& positiveVertices, std::vector<Triangle>& positiveFaces,
	std::vector<Vertex *>& negativeVertices, std::vector<Triangle>& negativeFaces,
	std::vector<Tetrahedron>& deleted, std::vector<Tetrahedron>& inserted, 
	bool refineFacetTest, bool skinnyTest) {

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
		deleted, refineFacetTest)) {
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
			for (auto t : inserted) {
				meshRep.addTetrahedron(t.v[0], t.v[1], t.v[2], t.v[3]);
				if (skinnyTest && skinnyTetTest(t)) maySkinnyTets.push_back(t);
			}
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
		Vertex *newVert = cavityRep.allocVertex(oldVert->point, oldVert->weight);
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
						Vertex *newOppo = cavityRep.allocVertex(oppo->point, oppo->weight);
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

void DelMesher::refineRegion(const Triangle& regionFacet, bool missingSegTest, bool refineSegTest,
	bool setInsertionRadius, bool skinnyTest) {
	Vertex *steinerVert = allocVertex((regionFacet.v[1]->point + regionFacet.v[2]->point) * DelReal(0.5), 
		DelReal(0), VertexType::Vertex_FreeFacet);

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
		toBeDeleted = findPosition(steinerVert, toBeDeleted, meshRep);

		VolumeVertexInsertionFlags vifs;
		vifs.cdt = true; vifs.trueInsertion = false;
		vifs.missingSegTest = missingSegTest; vifs.missingFacetTest = true;
		vifs.refineSegTest = refineSegTest; vifs.insertRadiusTest = setInsertionRadius;
		vifs.skinnyTetTest = skinnyTest;
		insertVertex(steinerVert, toBeDeleted, meshRep, vifs);

		for (auto s : newSegsOfFaces) {
			surfaceRep.addTriangle(steinerVert, s.v[0], s.v[1], s.index);
			mayMissingFacets.push_back(Triangle(steinerVert, s.v[0], s.v[1]));
		}

		if (setInsertionRadius) {
			steinerVert->relaxedInsetionRadius = (vifs.parent->point - steinerVert->point).length();
			DelReal parentRIS = vifs.parent->relaxedInsetionRadius;
			if (matchVertexFlag(vifs.parent->getVertexType(), VertexType::Vertex_Segment))
				parentRIS *= DelReal(SQRT_2);
			steinerVert->relaxedInsetionRadius = std::max(steinerVert->relaxedInsetionRadius,
				parentRIS);
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
		if (matchVertexFlag(a->getVertexType(), VertexType::Vertex_Refined) ||
			matchVertexFlag(b->getVertexType(), VertexType::Vertex_Refined))
			type = VertexType(type | VertexType::Vertex_Refined);

		steinerVert = allocVertex((a->point + b->point) * DelReal(0.5), DelReal(0), type);

		//split the encroched segment and update related information
		Segment s0(a, steinerVert), s1(b, steinerVert);
		if (missingSegTest) {
			mayMissingSegs.push_back(s0);
			mayMissingSegs.push_back(s1);
		}
		if (refineSegTest) {
			mayRefineSegs.push_back(s0);
			mayRefineSegs.push_back(s1);
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
		toBeDeleted = findPosition(steinerVert, toBeDeleted, meshRep);

		VolumeVertexInsertionFlags newVifs;
		newVifs.cdt = true; newVifs.missingSegTest = true; newVifs.missingFacetTest = true;
		newVifs.insertRadiusTest = setInsertionRadius; newVifs.skinnyTetTest = skinnyTest;
		insertVertexOnSegment(steinerVert, encrochedSeg, toBeDeleted, meshRep, newVifs);

		if (setInsertionRadius) {
			steinerVert->relaxedInsetionRadius = (newVifs.parent->point - steinerVert->point).length();
			if (Adjacent(encrochedSeg, newVifs.parent)) {
				DelReal parentRIS = newVifs.parent->relaxedInsetionRadius;
				if (!matchVertexFlag(newVifs.parent->getVertexType(), VertexType::Vertex_Segment) &&
					matchVertexFlag(newVifs.parent->getVertexType(), VertexType::Vertex_Refined))
					parentRIS = std::max((steinerVert->point - a->point).length(), parentRIS / DelReal(SQRT_2));
				steinerVert->relaxedInsetionRadius = std::max(steinerVert->relaxedInsetionRadius,
					parentRIS);
			}
			if (newVifs.parent == a || newVifs.parent == b)
				steinerVert->relaxedInsetionRadius = std::max(steinerVert->relaxedInsetionRadius,
					std::min(a->relaxedInsetionRadius, b->relaxedInsetionRadius));
		}
	}
}

void DelMesher::insertVertex(Vertex *u, const Tetrahedron& tet, TetMeshDataStructure& meshRep,
	const VolumeVertexInsertionFlags& vifs, Tetrahedron *rt){
	Vertex *a = tet.v[0], *b = tet.v[1], *c = tet.v[2], *d = tet.v[3];

	if (!vifs.trueInsertion) tobeDeletedTets.push_back(tet);
	meshRep.deleteTetrahedron(a, b, c, d);

	if (vifs.refineSegTest || vifs.missingSegTest) {
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

	digCavity(u, Triangle(b, d, c, true), meshRep, vifs, 0, rt);
	digCavity(u, Triangle(a, c, d, true), meshRep, vifs, 0);
	digCavity(u, Triangle(a, d, b, true), meshRep, vifs, 0);
	digCavity(u, Triangle(a, b, c, true), meshRep, vifs, 0);

	if (vifs.refineSegTest || vifs.missingSegTest) {
		for (auto s : markedSegments) {
			meshRep.unSetMark(s.v[0], s.v[1]);
			if (vifs.refineSegTest) mayRefineSegs.push_back(s);
			if (vifs.missingSegTest) mayMissingSegs.push_back(s);
		}
		markedSegments.clear();
	}

	if (vifs.cdt && vifs.trueInsertion) {
		if (!vifs.rejected) {
			for (auto f : newFacesOfTets) {
				meshRep.addTetrahedron(u, f.v[0], f.v[1], f.v[2]);
				if (vifs.skinnyTetTest)
					maySkinnyTets.push_back(Tetrahedron(u, f.v[0], f.v[1], f.v[2]));
			}
			newFacesOfTets.clear();
		}
		else {
			triangulateCavity(u, newFacesOfTets, meshRep, vifs);
			if (rt) meshRep.adjacent2Vertex(u, rt);
			newFacesOfTets.clear();
		}
	}
}

void DelMesher::insertVertexOnSegment(Vertex *u, const Segment& s, const Tetrahedron& tet, TetMeshDataStructure& meshRep, const VolumeVertexInsertionFlags& vifs){
	bool onFacet = matchVertexFlag(u->getVertexType(), VertexType::Vertex_Facet);
	if (onFacet) {
		SurfaceVertexInsertionFlags svifs;
		svifs.missingTest = vifs.missingFacetTest;
		svifs.refineTest = vifs.refineFacetTest;
		svifs.trueInsertion = vifs.trueInsertion;
		if (vifs.cdt) svifs.trueInsertion = false;
		if (!vifs.trueInsertion) svifs.deletionRecord = true;
		insertSurfaceSegmentVertex(u, s, svifs);
	}
	if (vifs.trueInsertion) {
		meshRep.deleteSegment(s.v[0], s.v[1]);
		meshRep.addSegment(u, s.v[0]);
		meshRep.addSegment(u, s.v[1]);
	}
	insertVertex(u, tet, meshRep, vifs);
	if (onFacet && vifs.cdt) {
		if (vifs.trueInsertion) {
			surfaceRep.bindVolumeVertex(u);
			for (auto seg : newSegsOfFaces) surfaceRep.addTriangle(u, seg.v[0], seg.v[1], seg.index);
			newSegsOfFaces.clear();
		}
	}
}

void DelMesher::insertVertexOnSurface(Vertex *u, const Triangle& f, const Tetrahedron& tet,
	TetMeshDataStructure& meshRep, const VolumeVertexInsertionFlags& vifs) {
	SurfaceVertexInsertionFlags svifs;
	svifs.missingTest = vifs.missingFacetTest;
	svifs.refineTest = vifs.refineFacetTest;
	if (vifs.cdt) svifs.trueInsertion = false;
	if (!vifs.trueInsertion) svifs.deletionRecord = true;
	insertSurfaceVertex(u, f, svifs);

	insertVertex(u, tet, meshRep, vifs);

	if (vifs.cdt) {
		if (vifs.trueInsertion) {
			surfaceRep.bindVolumeVertex(u);
			for (auto seg : newSegsOfFaces) surfaceRep.addTriangle(u, seg.v[0], seg.v[1], seg.index);
			newSegsOfFaces.clear();
		}
	}
}

void DelMesher::triangulateCavity(Vertex *u, std::vector<Triangle>& boundaries, TetMeshDataStructure& meshRep, const VolumeVertexInsertionFlags& vifs) {
	Assert(vifs.trueInsertion);
	std::vector<Vertex *> cavityVertices;
	cavityVertices.reserve(boundaries.size());

	TetMeshDataStructure cavityRep;
	cavityRep.Reserve(boundaries.size());
	Vertex *nu = cavityRep.allocVertex(u->point, u->weight);
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
					newVert = cavityRep.allocVertex(f.v[i]->point, f.v[i]->weight);
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
						predicator.inSpherePerturbed(oppo->point, f.v[0]->point, f.v[2]->point,
						f.v[1]->point, newOppo->point) > DelReal(0)) {
						enlarge = true;

						meshRep.deleteTetrahedron(oppo, f.v[0], f.v[2], f.v[1]);
						if (!vifs.trueInsertion) tobeDeletedTets.push_back(Tetrahedron(oppo, f.v[0], f.v[2], f.v[1]));

						if (oppo->getPointedVertex() == NULL) {
							Vertex *newVert = cavityRep.allocVertex(oppo->point, oppo->weight);
							newVert->setVertexPointer(oppo); oppo->setVertexPointer(newVert);
							cavityVertices.push_back(newVert);
							insertVertex(newVert, findPosition(newVert, Tetrahedron(newOppo, face.v[0], face.v[1], face.v[2]),
								cavityRep), cavityRep);

							if (vifs.insertRadiusTest) {
								DelReal minLen2 = (vifs.parent != NULL) ? (vifs.parent->point - u->point).length2()
									: std::numeric_limits<DelReal>::max();
								if (minLen2 > (oppo->point - u->point).length2()) vifs.parent = oppo;
							}
						}

						cavityShells.push_back(Triangle(oppo, back.v[0], back.v[2]));
						cavityShells.push_back(Triangle(oppo, back.v[1], back.v[0]));
						cavityShells.push_back(Triangle(oppo, back.v[2], back.v[1]));

						if (vifs.missingSegTest || vifs.refineSegTest) {
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
					if (vifs.missingFacetTest || vifs.refineFacetTest) {
						if (surfaceRep.Contain(f)) {
							if (vifs.missingFacetTest) mayMissingFacets.push_back(f);
							if (vifs.refineFacetTest) mayRefineFacets.push_back(f);
						}
						else if (surfaceRep.Contain(back)) {
							if (vifs.missingFacetTest) mayMissingFacets.push_back(back);
							if (vifs.refineFacetTest) mayRefineFacets.push_back(back);
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
						Vertex *newVert = cavityRep.allocVertex(oppo->point, oppo->weight);
						newVert->setVertexPointer(oppo); oppo->setVertexPointer(newVert);
						cavityVertices.push_back(newVert);
						Tetrahedron hint;
						bool found = cavityRep.adjacent2Vertex(face.v[0], &hint);
						Assert(found);
						hint = findPosition(newVert, hint, cavityRep);
						insertVertex(newVert, hint, cavityRep);

						if (vifs.insertRadiusTest) {
							DelReal minLen2 = (vifs.parent != NULL) ? (vifs.parent->point - u->point).length2()
								: std::numeric_limits<DelReal>::max();
							if (minLen2 > (oppo->point - u->point).length2()) vifs.parent = oppo;
						}
					}

					cavityShells.push_back(Triangle(oppo, back.v[0], back.v[2]));
					cavityShells.push_back(Triangle(oppo, back.v[1], back.v[0]));
					cavityShells.push_back(Triangle(oppo, back.v[2], back.v[1]));


					if (vifs.missingSegTest || vifs.refineSegTest) {
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
					if (vifs.missingFacetTest || vifs.refineFacetTest) {
						if (surfaceRep.Contain(f)) {
							if (vifs.missingFacetTest) mayMissingFacets.push_back(f);
							if (vifs.refineFacetTest) mayRefineFacets.push_back(f);
						}
						else if (surfaceRep.Contain(back)) {
							if (vifs.missingFacetTest) mayMissingFacets.push_back(back);
							if (vifs.refineFacetTest) mayRefineFacets.push_back(back);
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
	Triangle boundaryFacetOut = Triangle(cavityShells[0].v[0]->getPointedVertex(), cavityShells[0].v[2]->getPointedVertex(),
		cavityShells[0].v[1]->getPointedVertex());
	propagateCleanCavity(boundaryFacetOut, cavityRep, 0);
	//collect tets
	for (auto t : cavityRep) {
		meshRep.addTetrahedron(t.v[0]->getPointedVertex(), t.v[1]->getPointedVertex(),
			t.v[2]->getPointedVertex(), t.v[3]->getPointedVertex());
		if (vifs.skinnyTetTest) 
			maySkinnyTets.push_back(Tetrahedron(t.v[0]->getPointedVertex(), t.v[1]->getPointedVertex(),
				t.v[2]->getPointedVertex(), t.v[3]->getPointedVertex()));
	}

	if (vifs.refineFacetTest) {
		for (auto f : cavityShells) {
			Triangle back(f.v[0], f.v[2], f.v[1]);
			if (surfaceRep.Contain(f)) mayRefineFacets.push_back(f);
			if (surfaceRep.Contain(back)) mayRefineFacets.push_back(back);
		}
	}

	if (vifs.refineSegTest || vifs.missingSegTest) {
		for (auto s : markedSegments) {
			meshRep.unSetMark(s.v[0], s.v[1]);
			if (vifs.refineSegTest) mayRefineSegs.push_back(s);
			if (vifs.missingSegTest) mayMissingSegs.push_back(s);
		}
		markedSegments.clear();
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
					cavityVert = cavityRep.allocVertex(f.v[i]->point, f.v[i]->weight);
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
	Triangle boundaryFacetOut = Triangle(boundary[0].v[0]->getPointedVertex(), boundary[0].v[2]->getPointedVertex(),
		boundary[0].v[1]->getPointedVertex());
	propagateCleanCavity(boundaryFacetOut, cavityRep, 0);

	for (auto t : cavityRep) {
		meshRep.addTetrahedron(t.v[0]->getPointedVertex(), t.v[1]->getPointedVertex(),
			t.v[2]->getPointedVertex(), t.v[3]->getPointedVertex());

		if (vifs.skinnyTetTest) 
			maySkinnyTets.push_back(Tetrahedron(t.v[0]->getPointedVertex(), t.v[1]->getPointedVertex(),
				t.v[2]->getPointedVertex(), t.v[3]->getPointedVertex()));
	}
	
	for (auto v : cavityVertices) (v->getPointedVertex())->setVertexPointer(NULL);
	if (hasGhost) meshRep.getGhostVertex()->setVertexPointer(NULL);
}

void DelMesher::insertSurfaceVertex(Vertex *u, const Triangle& f, const SurfaceVertexInsertionFlags& vifs){
	Vertex *a = f.v[0], *b = f.v[1], *c = f.v[2];
	int faceIndex = surfaceRep.getTriangleIndex(a, b, c);
	if (vifs.trueInsertion) surfaceRep.bindVolumeVertex(u);
	if (vifs.deletionRecord) tobeDeletedFaces.push_back(TriangleWithIndex(a, b, c, faceIndex));
	surfaceRep.deleteTriangle(a, b, c);
	DelVector above = abovePoints[faceIndex];
	digCavity(u, Segment(b, a), above, faceIndex, vifs);
	digCavity(u, Segment(c, b), above, faceIndex, vifs);
	digCavity(u, Segment(a, c), above, faceIndex, vifs);
}

void DelMesher::insertSurfaceSegmentVertex(Vertex *u, const Segment &s, const SurfaceVertexInsertionFlags& vifs) {
	if (vifs.trueInsertion) surfaceRep.bindVolumeVertex(u);
	Vertex *a = s.v[0], *b = s.v[1];
	Vertex *oppo = NULL;
	int index = -1;
	if (surfaceRep.Adjacent(s, &oppo, &index)) {
		surfaceRep.deleteTriangle(a, b, oppo);
		DelVector above = abovePoints[index];
		digCavity(u, Segment(a, oppo), above, index, vifs);
		digCavity(u, Segment(oppo, b), above, index, vifs);
	}
	if (surfaceRep.Adjacent(Segment(b, a), &oppo, &index)) {
		surfaceRep.deleteTriangle(b, a, oppo);
		DelVector above = abovePoints[index];
		digCavity(u, Segment(b, oppo), above, index, vifs);
		digCavity(u, Segment(oppo, a), above, index, vifs);
	}
}

//Also mark all vertices
void DelMesher::detectCoplanarFaces(const Triangle& f, DelReal facetRadianTol,
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

void DelMesher::propagateDetectCoplanarFaces(Vertex *ref, const Segment& s, DelReal facetRadianTol,
	std::vector<Vertex *>& coplanarVertices, std::vector<Segment>& boundaries, TriMeshDataStructure& surfRep, int depth) const {
	Vertex *u = s.v[0], *v = s.v[1], *w = NULL;
	if (!surfRep.Adjacent(s, &w)) {
		if (u->getPointedVertex() && v->getPointedVertex() &&
			findSegment(Segment(u->getPointedVertex(), v->getPointedVertex()))) 
			boundaries.push_back(Segment(v, u));
		return;
	}

	if (Geometer::dihedralAngle(u->point, v->point, ref->point, w->point) > facetRadianTol) {
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
		if (surfaceRep.Contain(f) || !meshRep.Adjacent(f, &a))
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
			if (!surfaceRep.Contain(top) && meshRep.Adjacent(top, &a)) {
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

bool DelMesher::Encroached(const Segment& s, Vertex **encroachedVert) const {
	Vertex *a = s.v[0], *b = s.v[1], *c = NULL;
	auto begin = meshRep.getIncidentFacets(s);
	if (begin == NULL) return true;

	auto cirulator = begin;
	do {
		c = cirulator->v[2];
		if (!c->isGhost() && predicator.inDiametricBall(c->point, a->point, b->point) <= 0) {
			if (encroachedVert) *encroachedVert = c;
			return true;
		}

		++cirulator;
	} while (cirulator != begin);

	return false;
}

bool DelMesher::Encroached(TriangleWithGeometry& f, Vertex **encroachedVert) const {
	f.setGeometricProperties();
	int facetIndex = surfaceRep.getTriangleIndex(f.v[0], f.v[1], f.v[2]);
	Assert(facetIndex >= 0);
	Triangle f2(f.v[0], f.v[2], f.v[1], true);
	Vertex *x = NULL, *y = NULL;
	meshRep.Adjacent(f, &x); meshRep.Adjacent(f2, &y);
	Assert(x != NULL && y!= NULL);

	constexpr DelReal epsilon = DelReal(1e-8);
	DelVector center = f.getCircumcenter();
	DelReal r2 = f.getRadius() * f.getRadius();
	Assert(r2 > 0);

	bool encroached0 = false, encroached1 = false;
	if (Encroached(center, r2, x)) {
		//determine whether the encroached vertex is eligible or not
		if (mayIneligible(facetIndex, x)) {
			DelVector n = Normalize(Geometer::triangleNormal(f.v[0]->point, f.v[1]->point, f.v[2]->point));
			DelReal dis = fabs((x->point - center) * n);
			if (matchVertexFlag(x->getVertexType(), VertexType::Vertex_Segment))
				encroached0 = (dis - DelReal(SQRT_2) * x->relaxedInsetionRadius >= epsilon * dis);
			else
				encroached0 = (dis - x->relaxedInsetionRadius >= epsilon * dis);
		}
		else
			encroached0 = true;
	}
	if (!encroached0 && Encroached(center, r2, y)) {
		//determine whether the encroached vertex is eligible or not
		if (mayIneligible(facetIndex, y)) {
			DelVector n = Normalize(Geometer::triangleNormal(f.v[0]->point, f.v[1]->point, f.v[2]->point));
			DelReal dis = fabs((y->point - center) * n);
			if (matchVertexFlag(y->getVertexType(), VertexType::Vertex_Segment))
				encroached1 = (dis - DelReal(SQRT_2) * y->relaxedInsetionRadius >= epsilon * dis);
			else
				encroached1 = (dis -  y->relaxedInsetionRadius >= epsilon * dis);
		}
		else
			encroached1 = true;
	}

	if (encroached0) *encroachedVert = x;
	else if(encroached1) *encroachedVert = y;

	return encroached0 || encroached1;
}

bool DelMesher::Encroached(TriangleWithGeometry &f, Vertex *v) const {
	f.setGeometricProperties();
	DelVector center = f.getCircumcenter();
	DelReal r2 = f.getRadius() * f.getRadius();
	Assert(r2 > 0);

	constexpr DelReal epsilon = DelReal(1e-8);
	if (Encroached(center, r2, v)) {
		//determine whether the encroached vertex is eligible or not
		int facetIndex = surfaceRep.getTriangleIndex(f.v[0], f.v[1], f.v[2]);
		Assert(facetIndex >= 0);
		if (mayIneligible(facetIndex, v)) {
			DelVector n = Normalize(Geometer::triangleNormal(f.v[0]->point, f.v[1]->point, f.v[2]->point));
			DelReal dis = fabs((v->point - center) * n);
			if (matchVertexFlag(v->getVertexType(), VertexType::Vertex_Segment))
				return dis - DelReal(SQRT_2) * v->relaxedInsetionRadius >= epsilon * dis;
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

bool DelMesher::Adjacent(int facetIndex, Vertex *v) const {
	Vertex **verts = NULL;
	size_t size = getPolygonVertices(facetIndex, &verts);

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
		bool flag = surfaceRep.adjacent2Vertex(v, &f2, &f2Index);

		if (!flag || facetIndex == f2Index) return false;

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

bool DelMesher::mayIneligible(int fIndex, Vertex *v) const {
	Vertex **verts = NULL;
	size_t size = getPolygonVertices(fIndex, &verts);

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
				if (verts[i] == v || findSegment(Segment(v, verts[i]))) {
					adjacent = true;
					break;
				}
			}
		}
	}
	else if (matchVertexFlag(v->getVertexType(), VertexType::Vertex_Facet)) {
		Triangle f2;
		int f2Index = -1;
		surfaceRep.adjacent2Vertex(v, &f2, &f2Index);
		Assert(f2Index >= 0);

		if (fIndex == f2Index) return true;

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

Reference<Mesh> DelMesher::generateMesh(int *vertexLableMap){
	detectAcuteVertices();

	//allocate space
	tobeDeletedFaces.reserve(32);
	newSegsOfFaces.reserve(32);
	tobeDeletedTets.reserve(32);
	newFacesOfTets.reserve(32);

	constrainedTriangulation();

	for (auto t : meshRep) {
		//set vertices relaxed insertion radius
		Vertex *a = t.v[0], *b = t.v[1], *c = t.v[2], *d = t.v[3];
		if (!a->isGhost() && !b->isGhost() && !c->isGhost() && !d->isGhost()) {
			DelVector ab = b->point - a->point, ac = c->point - a->point, ad = d->point - a->point,
				bc = c->point - b->point, bd = d->point - b->point, cd = d->point - c->point;
			DelReal abLen = ab.length(), acLen = ac.length(), adLen = ad.length(),
				bcLen = bc.length(), bdLen = bd.length(), cdLen = cd.length();

			a->relaxedInsetionRadius = std::min({ a->relaxedInsetionRadius, abLen, acLen, adLen });
			b->relaxedInsetionRadius = std::min({ b->relaxedInsetionRadius, abLen, bcLen, bdLen });
			c->relaxedInsetionRadius = std::min({ c->relaxedInsetionRadius, acLen, bcLen, cdLen });
			d->relaxedInsetionRadius = std::min({ d->relaxedInsetionRadius, adLen, bdLen, cdLen });
		}
	}

	constrainedRefinment();

	Triangle out = *(surfaceRep.begin()); std::swap(out.v[0], out.v[1]);
	propagateClean(out, 0);

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
		mesh->setVertex(vertIndices[pair.second], Vector3{ Scalar(vertex->point.x), Scalar(vertex->point.y), Scalar(vertex->point.z) });
	}

	i = 0;
	int surfVertIndices[3];
	for (auto f : polygons){
		std::swap(f.v[0], f.v[1]);
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

	return mesh;
}
}
