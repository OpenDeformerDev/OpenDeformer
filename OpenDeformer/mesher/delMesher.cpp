#include "stdafx.h"
#include "delMesher.h"
#include "geometer.h"
#include <memory>
#include <set>
#include <random>
#include <stack>

namespace ODER{

Predicator<REAL> DelTriangulator::predicator;

void DelTriangulator::generateSubPolygons(Vertex **vertices, int *segments, int vertexCount, int segmentCount, const Face& ref, bool boundaryOnly){
	meshRep.Clear();
	meshRep.Reserve(vertexCount);
	this->segments.clear();

	Vertex *v0, *v1, *v2, *v3;

#ifdef ODER_DEBUG
	if (vertexCount < 3)
		Severe("a subpolygon contains less than 3 vertices");
#endif

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
		this->segments.insert(Segment(vertices[segments[i * 2 + 0]], vertices[segments[i * 2 + 1]]));

	calculateAbovePoint(vertexCount, vertices, ref);
	v0 = vertices[0];
	v1 = vertices[1];
	v2 = vertices[2];

	Face f(v0, v1, v2);
	//start
	meshRep.addTriangle(v0, v1, v2);
	meshRep.addTriangle(v1, v0, ghost);
	meshRep.addTriangle(v2, v1, ghost);
	meshRep.addTriangle(v0, v2, ghost);
	for (int i = 3; i < vertexCount; i++){
		v3 = vertices[i];
		f = findPosition(v3, f);
		v0 = f.v[0]; v1 = f.v[1]; v2 = f.v[2];
		meshRep.deleteTriangle(v0, v1, v2);
		digCavity(v3, Segment(v1, v0), &f, 0);
		digCavity(v3, Segment(v0, v2), NULL, 0);
		digCavity(v3, Segment(v2, v1), NULL, 0);
	}

	//insert segments
	for (auto s : this->segments) {
		Vertex *x = NULL;
		if (!meshRep.Contain(s))
			insertSegment(s);
	}

	//clean out of region triangles
	if (boundaryOnly) {
		for (auto s : this->segments) {
			std::swap(s.v[0], s.v[1]);
			propagateClean(s, 0);
		}
	}
}

void DelTriangulator::generateSubPolygons(Vertex **vertices, Segment *segments, int vertexCount, int segmentCount, const Face& ref, bool boundaryOnly) {
	meshRep.Clear();
	meshRep.Reserve(vertexCount);
	this->segments.clear();

	Vertex *v0, *v1, *v2, *v3;

#ifdef ODER_DEBUG
	if (vertexCount < 3)
		Severe("a subpolygon contains less than 3 vertices");
#endif

	//shuffle input segments
	std::random_device rng;
	std::default_random_engine randomEngine(rng());
	std::shuffle(segments, segments + segmentCount, randomEngine);

	for (int i = 0; i < segmentCount; i++)
		this->segments.insert(segments[i]);

	calculateAbovePoint(vertexCount, vertices, ref);
	v0 = vertices[0];
	v1 = vertices[1];
	v2 = vertices[2];

	Face f(v0, v1, v2);
	//start
	meshRep.addTriangle(v0, v1, v2);
	meshRep.addTriangle(v1, v0, ghost);
	meshRep.addTriangle(v2, v1, ghost);
	meshRep.addTriangle(v0, v2, ghost);
	for (int i = 3; i < vertexCount; i++) {
		v3 = vertices[i];
		f = findPosition(v3, f);
		v0 = f.v[0]; v1 = f.v[1]; v2 = f.v[2];
		meshRep.deleteTriangle(v0, v1, v2);
		digCavity(v3, Segment(v1, v0), &f, 0);
		digCavity(v3, Segment(v0, v2), NULL, 0);
		digCavity(v3, Segment(v2, v1), NULL, 0);
	}

	//insert segments
	for (auto s : this->segments) {
		Vertex *x = NULL;
		if (!meshRep.Contain(s))
			insertSegment(s);
	}

	//clean out of region triangles
	if (boundaryOnly) {
		for (auto s : this->segments) {
			std::swap(s.v[0], s.v[1]);
			propagateClean(s, 0);
		}
	}
}

void DelTriangulator::insertSegments(const Face *triangles, const Segment *segments, int triangleCount, int segmentCount) {
	meshRep.Clear();
	meshRep.Reserve(triangleCount);
	this->segments.clear();

	REAL maxArea2 = 0;
	Face largetsFace;
	for (int i = 0; i < triangleCount; i++) {
		Face f = triangles[i];
		meshRep.addTriangle(f.v[0], f.v[1], f.v[2]);
		REAL area2 = (Geometer::triangleNormal(f.v[0]->vert, f.v[1]->vert, f.v[2]->vert)).length2();
		if (area2 > maxArea2) {
			maxArea2 = area2;
			largetsFace = f;
		}
	}
	//calculate above point
	ghost->vert = Geometer::calculateAbovePoint(largetsFace.v[0]->vert, largetsFace.v[1]->vert, largetsFace.v[2]->vert);

	for (int i = 0; i < segmentCount; i++) {
		Segment s = segments[i];
		this->segments.insert(s);
		if (!meshRep.Contain(s))
			insertSegment(s);
	}
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

void DelTriangulator::outPut(TriMeshDataStructure& meshRep, std::deque<Face>& meshDeque) {
	for (auto f : this->meshRep) {
		meshRep.addTriangle(f.v[0], f.v[1], f.v[2]);
		meshDeque.push_back(f);
	}
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
	
	Vertex *start = vertices[size - 1], *end = vertices[0];
	prev[0] = size - 1; prev[size - 1] = size - 2;
	next[0] = 1; next[size - 1] = 0;
	distance[0] = distance[size - 1] = 0;

	DelVector ab = end->vert - start->vert;
	for (int i = 1; i < size - 1; i++) {
		prev[i] = i - 1;
		next[i] = i + 1;
		distance[i] = ((vertices[i]->vert - start->vert) % ab).length2();
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
	cavityRep.addTriangle(start, end, vertices[indices[0]]);


	std::deque<Vertex *> fanVertices;
	std::vector<Vertex *> convexPoly;
	for (int i = 1; i < indexCount; i++) {
		int trueIndex = indices[i];
		bool convexPolyOri = true;
		insertVertexToCavity(vertices[trueIndex], vertices[next[trueIndex]], vertices[prev[trueIndex]],
			false, true, fanVertices, 0);

		Assert(fanVertices.size() != 1);
		if (fanVertices.size() > 2) {
			convexPoly.reserve(fanVertices.size());
			for (auto vert : fanVertices)
				convexPoly.push_back(vert);
			triangulateConvexPoly(vertices[trueIndex], convexPoly);

			fanVertices.clear();
			convexPoly.clear();
		}
	}

	for (auto f : cavityRep)
		meshRep.addTriangle(f.v[0], f.v[1], f.v[2]);

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
	for (int i = 0; i < convexPoly.size(); i++)
		convexPoly[i]->setMark();

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
	for (int i = 0; i < convexPoly.size(); i++)
		convexPoly[i]->unSetMark();

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

void DelTriangulator::calculateAbovePoint(int vertexCount, Vertex** vertices, const Face& ref) {
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
	ghost->vert = a->vert + n;
	if (predicator.orient3d(ghost->vert, ref.v[0]->vert, ref.v[1]->vert, ref.v[2]->vert) < 0) {
		ghost->vert = a->vert - n;
		std::swap(vertices[1], vertices[2]);
	}
}

Predicator<REAL> DelMesher::predicator;

DelMesher::DelMesher(Vector *surfvs, int *triangles, int numv, int numtri, REAL maxR, REAL maxRa, REAL facetAngleTol) {
	maxRatio = maxRa;
	maxRadius = maxR;
	Tetrahedron::maxREration = maxRa;

	//init ghost vertex
	ghost = allocVertex();
	ghost->setGhost();
	//insert every vertex
	oriVertices.reserve(numv);
	for (int i = 0; i < numv; i++) {
		boundBox.Insert(surfvs[i]);
		Vertex *vert = allocVertex(DelVector{ surfvs[i].x, surfvs[i].y, surfvs[i].z}, REAL(0),
			VertexType(VertexType::Vertex_FixedVolume | VertexType::Vertex_Facet));
		oriVertices.push_back(vert);
	}

	TriMeshDataStructure oriSurface;
	std::vector<Face> oriTriangles;
	oriTriangles.reserve(numtri);
	for (int i = 0; i < numtri; i++) {
		oriSurface.addTriangle(oriVertices[triangles[3 * i]], oriVertices[triangles[3 * i + 1]], oriVertices[triangles[3 * i + 2]]);
		oriTriangles.push_back(Face(oriVertices[triangles[3 * i]], oriVertices[triangles[3 * i + 1]], oriVertices[triangles[3 * i + 2]]));
	}

	REAL facetRadianTol = facetAngleTol / REAL(180) * REAL(M_PI);
	std::vector<Vertex *> coplanarVertices;
	std::vector<Segment> boundaries;
	DelTriangulator triangulator;

	for (auto f : oriTriangles) {
		if (oriSurface.Contain(f)) {
			detectCoplanarFaces(f, facetRadianTol, coplanarVertices, boundaries, oriSurface);
			for (auto s : boundaries) {
				s.v[0]->type = VertexType(s.v[0]->type | VertexType::Vertex_Segment);
				s.v[1]->type = VertexType(s.v[1]->type | VertexType::Vertex_Segment);
				segments.insert(Segment(s.v[0], s.v[1], true));
			}


			//re-triangulation
			if (coplanarVertices.size() > 3) {
				triangulator.generateSubPolygons(&coplanarVertices[0], &boundaries[0], coplanarVertices.size(), boundaries.size(), f, true);
				triangulator.outPut(surfaceRep, mayEncroachedFaces);
			}
			else {
				surfaceRep.addTriangle(f.v[0], f.v[1], f.v[2]);
				mayEncroachedFaces.push_back(f);
			}

			coplanarVertices.clear();
			boundaries.clear();
		}
	}	
	for (auto s : this->segments) mayEncroachedSegs.push_back(s);

	//allocate space
	tobeDeletedFaces.reserve(16);
	newSegsOfFaces.reserve(32);
	tobeDeletedTets.reserve(16);
	newTets.reserve(32);

	constrainedTriangulation();
}

void DelMesher::triangulation3D(std::vector<Vertex *>& vertices, TetMeshDataStructure& meshRep, bool insertToSkinny){
	Vertex *a = vertices[0], *b = vertices[1];
	Vertex *c = NULL, *d = NULL;
	for (int i = 2; i < vertices.size(); i++) {
		c = vertices[i];
		if (!predicator.fastCoLine(a->vert, b->vert, c->vert)) {
			std::swap(vertices[2], vertices[i]);
			break;
		}
	}
	for (int i = 3; i < vertices.size(); i++){
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

	meshRep.addTetrahedron(a, b, c, d);
	meshRep.addTetrahedron(ghost, b, d, c);
	meshRep.addTetrahedron(ghost, a, c, d);
	meshRep.addTetrahedron(ghost, a, b, c);
	meshRep.addTetrahedron(ghost, a, d, b);

	VertexInsertionFlags vifs;
	vifs.skinnyTetTest = insertToSkinny;
	Tetrahedron rt = t;
	for (int i = 4; i < vertices.size(); i++){
		Vertex *v = vertices[i];
		t = findPosition(v, rt, meshRep);
		insertVertex(v, t, meshRep, vifs, &rt);
	}
}

//all segments and faces to be recoverd must be in DelMesher::mayEncroachedSegs and DelMesher::mayEncroachedFaces
void DelMesher::constrainedTriangulation() {
	meshRep.Reserve(oriVertices.size());
	triangulation3D(oriVertices, meshRep, false);

	//working vectors
	std::vector<Vertex *> regionVertices, positiveVertices, negativeVertices;
	std::vector<Segment> regionBoundaries;
	std::vector<Face> regionFaces, positiveFaces, negativeFaces;
	std::vector<Tetrahedron> deleted, inserted;

start:
	while (!mayEncroachedSegs.empty()) {
		Segment s = mayEncroachedSegs.front();
		mayEncroachedSegs.pop_front();
		if (findSegment(s)) {
			Tetrahedron t;
			if (findIntersectedTetrahedron(s.v[0], s.v[1]->vert, &t)) {
				if (s.v[1] != t.v[1] && s.v[1] != t.v[2] && s.v[1] != t.v[3])
					splitSubSegment(s, findSegmentEncroachedReference(s.v[1], t));
			}
			else
				Assert(false);
		}
	}

	while (!mayEncroachedFaces.empty()) {
		Face f = mayEncroachedFaces.front();
		mayEncroachedFaces.pop_front();
		if (surfaceRep.Contain(f) && !meshRep.Contain(f)) {
			if (!faceRecovery(f, regionVertices, regionBoundaries, regionFaces,
				positiveVertices, positiveFaces, negativeVertices, negativeFaces,
				deleted, inserted)) {
				refineRegion(f);
				goto start;
			}
		}
	}

	propagateClean(*(surfaceRep.begin()), 0);
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

	while (maxIter--){
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
					    return Tetrahedron(ghost, d, c, b, true);
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
					return Tetrahedron(ghost, a, b, c, true);
				//Tetrahedron(x, a, b, c)
				d = c; c = b; b = a; a = x;
				break;
			case 2:
				meshRep.Adjacent(Face(a, c, d), &x);
				if (x->isGhost())
					return Tetrahedron(ghost, a, c, d, true);
				//Tetrahedron(x, a, c, d)
				b = a; a = x;
				break;
			case 3:
				if (Randomnation(2) == 0){
					meshRep.Adjacent(Face(a, b, c), &x);
					if (x->isGhost())
						return Tetrahedron(ghost, a, b, c, true);
					//Tetrahedron(x, a, b, c)
					d = c; c = b; b = a; a = x;
				}
				else{
					meshRep.Adjacent(Face(a, c, d), &x);
					if (x->isGhost())
						return Tetrahedron(ghost, a, c, d, true);
					//Tetrahedron(x, a, c, d)
					b = a; a = x;
				}
				break;
			case 4:
				meshRep.Adjacent(Face(a, d, b), &x);
				if (x->isGhost())
					return Tetrahedron(ghost, a, d, b, true);
				//Tetrahedron(x, a, d, b)
				c = d; d = b; b = a; a = x;
				break;
			case 5:
				if (Randomnation(2) == 0){
					meshRep.Adjacent(Face(a, b, c), &x);
					if (x->isGhost())
						return Tetrahedron(ghost, a, b, c, true);
					//Tetrahedron(x, a, b, c)
					d = c; c = b; b = a; a = x;
				}
				else{
					meshRep.Adjacent(Face(a, d, b), &x);
					if (x->isGhost())
						return Tetrahedron(ghost, a, d, b, true);
					//Tetrahedron(x, a, d, b)
					c = d; d = b; b = a; a = x;
				}
				break;
			case 6:
				if (Randomnation(2) == 0){
					meshRep.Adjacent(Face(a, c, d), &x);
					if (x->isGhost())
						return Tetrahedron(ghost, a, c, d, true);
					//Tetrahedron(x, a, c, d)
					b = a; a = x;
				}
				else{
					meshRep.Adjacent(Face(a, d, b), &x);
					if (x->isGhost())
						return Tetrahedron(ghost, a, d, b, true);
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
						return Tetrahedron(ghost, a, b, c, true);
					//Tetrahedron(x, a, b, c)
					d = c; c = b; b = a; a = x;
				}
				else if (term == 1){
					meshRep.Adjacent(Face(a, c, d), &x);
					if (x->isGhost())
						return Tetrahedron(ghost, a, c, d, true);
					b = a; a = x;
				}
				else{
					meshRep.Adjacent(Face(a, d, b), &x);
					if (x->isGhost())
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
	const VertexInsertionFlags& vifs, Tetrahedron *rt){
	Vertex *a;
	if (!meshRep.Adjacent(f, &a)) {
		if (vifs.encroachFaceTest && surfaceRep.Contain(f))
			mayEncroachedFaces.push_back(f);
		return true;
	}

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

	if (deleted){
		if (!vifs.trueInsertion)
			tobeDeletedTets.push_back(Tetrahedron(a, b, c, d));
		meshRep.deleteTetrahedron(a, b, c, d);

		if (vifs.encroachFaceTest && surfaceRep.Contain(f))
			mayEncroachedFaces.push_back(f);
		if (vifs.encroachSegTest) {
			Segment ab(a, b, true), ac(a, c, true), ad(a, d, true);
			if (findSegment(ab)) mayEncroachedSegsSet.insert(ab);
			if (findSegment(ac)) mayEncroachedSegsSet.insert(ac);
			if (findSegment(ad)) mayEncroachedSegsSet.insert(ad);
		}

		success &= digCavity(u, Face(a, c, d, true), meshRep, vifs, rt);
		success &= digCavity(u, Face(a, b, c, true), meshRep, vifs, rt);
		success &= digCavity(u, Face(a, d, b, true), meshRep, vifs, rt);
	}
	else{
		if (constrained) success = predicator.orient3d(u->vert, b->vert, d->vert, c->vert) >= REAL(0);
		if (vifs.trueInsertion) 
			if (success) meshRep.addTetrahedron(u, b, d, c);
		else 
			newTets.push_back(Tetrahedron(u, b, d, c));
		if (rt)
			*rt = Tetrahedron(u, b, d, c, true);
		if (vifs.encroachFaceTest & 2) {
			Face faceIn(b, d, c), faceOut(b, c, d);
			if (surfaceRep.Contain(faceIn)) mayEncroachedFaces.push_back(faceIn);
			if (surfaceRep.Contain(faceOut)) mayEncroachedFaces.push_back(faceOut);
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

void DelMesher::digCavity(Vertex *u, const DelVector& above,const Segment &s, bool insertToQueue, bool trulyDeleteOrAdd){
	Vertex *x;
	if (!surfaceRep.Adjacent(s, &x))
		return;

	Vertex *w = s.v[0];
	Vertex *v = s.v[1];

	if (!findSegment(Segment(w, v, true)) &&
		predicator.inOrthoCirclePerturbed(u->vert, u->weight, v->vert, v->weight, w->vert, w->weight, x->vert, x->weight, above) > 0){
		if (!trulyDeleteOrAdd)
			tobeDeletedFaces.push_back(Face(w, v, x, true));
		surfaceRep.deleteTriangle(w, v, x);
		digCavity(u, above, Segment(x, v), insertToQueue, trulyDeleteOrAdd);
		digCavity(u, above, Segment(w, x), insertToQueue, trulyDeleteOrAdd);
	}
	else{
		if (trulyDeleteOrAdd)
			surfaceRep.addTriangle(u, v, w);
		else
			newSegsOfFaces.push_back(Segment(v, w));
		if (insertToQueue)
			mayEncroachedFaces.push_back(Face(u, v, w, true));
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
					angle = maxAngle;
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
						angle = maxAngle;
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

void DelMesher::splitSubSegment(const Segment& s, Vertex* ref){
	Vertex *a = s.v[0], *b = s.v[1];
	DelVector aa = a->vert, bb = b->vert;
	DelVector ab = bb - aa;
	REAL abLen = ab.length();
	REAL invabLen = REAL(1) / abLen;

	//calulate steiner point
	DelVector steinerPoint;
	Segment oriSegment;
	if (matchVertexFlag(ref->type, VertexType::Vertex_FreeSegment)) {
		auto refOriSegFound = vertSegHash.find(ref);
		Assert(refOriSegFound != vertSegHash.end());
		Segment refOriSeg = refOriSegFound->second;

		Vertex *oriSegStart = a, *oriSegEnd = b;
		Vertex *refOriSegStart = refOriSeg.v[0], *refOriSegEnd = refOriSeg.v[1];

		if (oriSegStart == refOriSegStart || oriSegStart == refOriSegEnd) {
			DelVector oriVec = oriSegEnd->vert - oriSegStart->vert;
			REAL len0 = oriVec.length();
			REAL len1 = (ref->vert - oriSegStart->vert).length();
			steinerPoint = oriSegStart->vert + (len1 / len0) * oriVec;
		}
		else if (oriSegEnd == refOriSegStart || oriSegEnd == refOriSegEnd) {
			DelVector oriVec = oriSegStart->vert - oriSegEnd->vert;
			REAL len0 = oriVec.length();
			REAL len1 = (ref->vert - oriSegEnd->vert).length();
			steinerPoint = oriSegEnd->vert + (len1 / len0) * oriVec;
		}
		else {
			DelVector ar = ref->vert - aa;
			DelVector br = ref->vert - bb;
			REAL arLen = ar.length();
			REAL brLen = br.length();
			if (arLen < REAL(0.5) * abLen)
				steinerPoint = aa + (arLen * invabLen) * ab;
			else if(brLen < REAL(0.5) * abLen)
				steinerPoint = bb - (brLen * invabLen) * ab;
			else
				steinerPoint = REAL(0.5) * (aa + bb);
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
	if (t < REAL(0.2) || t > REAL(0.8))
		steinerPoint = REAL(0.5) * (aa + bb);

	Vertex *steinerVert = allocVertex(steinerPoint, REAL(0), VertexType::Vertex_FreeSegment);

	//insert two sub segments
	Segment s0(a, steinerVert, true);
	Segment s1(steinerVert, b, true);
	segments.erase(s);
	segments.insert(s0);
	segments.insert(s1);
	mayEncroachedSegs.push_back(s0);
	mayEncroachedSegs.push_back(s1);
	vertSegHash.insert(std::make_pair(steinerVert, Segment(a, b)));

	Tetrahedron toBeDeleted;
	meshRep.adjacent2Vertex(a, &toBeDeleted);
	toBeDeleted = findPosition(steinerVert, toBeDeleted, meshRep);

	VertexInsertionFlags vifs;
	vifs.encroachSegTest = true;
	insertVertex(steinerVert, toBeDeleted, meshRep, vifs);

	//inserte it to surface
	insertSurfaceSegmentVertex(steinerVert, s);
}

void DelMesher::splitSubSegment(const Segment &s){
	Vertex *a = s.v[0], *b = s.v[1];
	Vertex *mid = allocVertex((a->vert + b->vert) * REAL(0.5), REAL(0), VertexType::Vertex_FreeSegment);

	Vertex *c = NULL, *d = NULL;
	//should be orient inward
	Face f;
	insertSurfaceSegmentVertex(mid, s, &f);

	Tetrahedron toBeDeleted;
	Vertex *e = NULL;
	if (matchVertexFlag(mid->type, VertexType::Vertex_Facet) && meshRep.Adjacent(f, &e))
		toBeDeleted = Tetrahedron(e, f.v[0], f.v[1], f.v[2], true);
	else{
		meshRep.adjacent2Vertex(s.v[0], &toBeDeleted);
		toBeDeleted = findPosition(mid, toBeDeleted, meshRep);
	}
	VertexInsertionFlags vifs;
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

	digCavity(&center, above, Segment(b, a), false, false);
	digCavity(&center, above, Segment(a, c), false, false);
	digCavity(&center, above, Segment(c, b), false, false);

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
		VertexInsertionFlags vifs;
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

	VertexInsertionFlags vifs;
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
		if (Encroached(nf, center)){
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
	if (!surfaceRep.Adjacent(edge, &w)) return;

	Vertex *u = edge.v[0], *v = edge.v[1];
	Face f(u, v, w);
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
	Vertex *a = boundary.v[0], *b = boundary.v[1], *c = NULL, *d = NULL, *e = NULL;
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
			Severe("Unexpected case in DelMesher::propagateFindCavity");
			break;
		}
	}

	for (auto s : regionBoundaries) surfaceRep.unSetMark(s.v[0], s.v[1]);

	if (newSegments.size() > 0) {
		for (auto f : regionFaces) surfaceRep.deleteTriangle(f.v[0], f.v[1], f.v[2]);
		if (toBeErased.size() != regionFaces.size()) {
			DelTriangulator triangulator;
			triangulator.insertSegments(&regionFaces[0], &newSegments[0], regionFaces.size(), newSegments.size());
			for(auto f : toBeErased) surfaceRep.addTriangle(f.v[0], f.v[1], f.v[2]);
			toBeErased.clear();
			triangulator.outPut(toBeErased);
			regionFaces.clear();
			for (auto f : toBeErased) {
				if (!surfaceRep.Contain(f)) {
					surfaceRep.addTriangle(f.v[0], f.v[1], f.v[2]);
					mayEncroachedFaces.push_back(f);
					regionFaces.push_back(f);
				}
			}
		}
		else {
			for (auto f : toBeErased) 
				surfaceRep.addTriangle(f.v[0], f.v[1], f.v[2]);

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
		posSuccess = triangulateCavity(regionFaces, true, positiveFaces, positiveVertices, deleted, inserted, f);
		if (posSuccess) {
			for (auto& face : regionFaces) std::swap(face.v[1], face.v[2]);
			negativeVertices.reserve(negativeVertices.size() + regionVertices.size());
			for (auto v : regionVertices) negativeVertices.push_back(v);
			negSuccess = triangulateCavity(regionFaces, true, negativeFaces, negativeVertices, deleted, inserted, f);
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
bool DelMesher::triangulateCavity(const std::vector<Face>& regionFaces, int encroachFaceTest, std::vector<Face>& boundaryFaces,
	std::vector<Vertex *>& cavityVertices, std::vector<Tetrahedron>& deleted, std::vector<Tetrahedron>& inserted, Face& encroached) {

	std::vector<Face> shellFaces;
	shellFaces.reserve(boundaryFaces.size());

	//construct triangulation
	cavityRep.Clear();
	cavityRep.Reserve(cavityVertices.size());
	std::random_device rd;
	std::default_random_engine randomEngine(rd());
	std::shuffle(cavityVertices.begin(), cavityVertices.end(), randomEngine);
	triangulation3D(cavityVertices, cavityRep, false);

	int swapCount = 0;
	do {
		shellFaces.clear();
		bool enlarged = false;
		for (auto face : boundaryFaces) {
			if (!cavityRep.Contain(face)) {
				Vertex *oppo = NULL;
				if (meshRep.Adjacent(face, &oppo)) {
					enlarged = true;
					meshRep.deleteTetrahedron(oppo, face.v[0], face.v[1], face.v[2]);
					deleted.push_back(Tetrahedron(oppo, face.v[0], face.v[1], face.v[2]));
					if (!oppo->isMarked()) {
						Tetrahedron hint;
						bool found = cavityRep.adjacent2Vertex(face.v[0], &hint);
						Assert(found);
						hint = findPosition(oppo, hint, cavityRep);
						insertVertex(oppo, hint, cavityRep);
						oppo->setMark();

						cavityVertices.push_back(oppo);
					}
					shellFaces.push_back(Face(oppo, face.v[0], face.v[2]));
					shellFaces.push_back(Face(oppo, face.v[1], face.v[0]));
					shellFaces.push_back(Face(oppo, face.v[2], face.v[1]));
				}

				if (encroachFaceTest) {
					if (surfaceRep.Contain(face)) mayEncroachedFaces.push_back(face);
					Face back(face.v[0], face.v[2], face.v[1]);
					if (surfaceRep.Contain(back)) mayEncroachedFaces.push_back(back);
				}
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
		//check region faces in cavity and mark it
		if (cavityRep.Contain(f))
			cavityRep.setMark(f.v[0], f.v[1], f.v[2]);
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
			cavityRep.setMark(f.v[0], f.v[1], f.v[2]);
		//clean tets out of cavity
		for (auto f : shellFaces)
			propagateCleanCavity(Face(f.v[0], f.v[2], f.v[1]), 0);
		for(auto f : regionFaces)
			propagateCleanCavity(Face(f.v[0], f.v[2], f.v[1]), 0);

		cavityRep.getTetrahedrons(false, inserted);

		if (encroachFaceTest & 2) {
			for (auto f : shellFaces) {
				if (surfaceRep.Contain(f)) mayEncroachedFaces.push_back(f);
				Face back(f.v[0], f.v[2], f.v[1]);
				if (surfaceRep.Contain(back)) mayEncroachedFaces.push_back(back);
			}
		}
	}

	if (swapCount & 1) std::swap(shellFaces, boundaryFaces);
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

void DelMesher::refineRegion(const Face& regionFace) {
	Vertex *steinerVert = allocVertex((regionFace.v[1]->vert + regionFace.v[2]->vert) * REAL(0.5), 
		REAL(0), VertexType::Vertex_FreeFacet);

	Tetrahedron toBeDeleted;
	meshRep.adjacent2Vertex(regionFace.v[1], &toBeDeleted);
	toBeDeleted = findPosition(regionFace.v[1], toBeDeleted, meshRep);

	Assert(mayEncroachedSegs.size() == 0);
	VertexInsertionFlags vifs;
	vifs.cdt = true; vifs.encroachSegTest = true; vifs.encroachFaceTest = 1; vifs.trueInsertion = false;
	insertVertex(steinerVert, toBeDeleted, meshRep, vifs);
	//detect encroachment
	bool encroach = false;
	Segment encrochedSeg;
	while (!mayEncroachedSegs.empty()) {
		encrochedSeg = mayEncroachedSegs.back();
		if (Encroached(encrochedSeg, *steinerVert)) {
			encroach = true;
			break;
		}
		mayEncroachedSegs.pop_back();
	}

	if (!encroach) {
		for (auto t : newTets)
			meshRep.addTetrahedron(t.v[0], t.v[1], t.v[2], t.v[3]);
		insertSurfaceVertex(steinerVert, regionFace, true);
	}
	else {
		//fall back
		for (auto t : tobeDeletedTets)
			meshRep.addTetrahedron(t.v[0], t.v[1], t.v[2], t.v[3]);

		mayEncroachedSegs.clear();

		//split encroch segment in mid point and update related information
		steinerVert->vert = (encrochedSeg.v[0]->vert + encrochedSeg.v[1]->vert) * REAL(0.5);
		steinerVert->weight = REAL(0);
		steinerVert->type = VertexType::Vertex_FreeSegment;

		Segment s0(encrochedSeg.v[0], steinerVert, true), s1(encrochedSeg.v[1], steinerVert, true);
		segments.erase(encrochedSeg);
		segments.insert(s0);
		segments.insert(s1);
		mayEncroachedSegs.push_back(s0);
		mayEncroachedSegs.push_back(s1);

		vertSegHash.insert(std::make_pair(steinerVert, encrochedSeg));

		bool found = findIntersectedTetrahedron(encrochedSeg.v[0], encrochedSeg.v[1]->vert, &toBeDeleted);
		Assert(found);

		VertexInsertionFlags svifs;
		svifs.cdt = true; svifs.encroachSegTest = true; svifs.encroachFaceTest = 1;
		insertVertex(steinerVert, toBeDeleted, meshRep, svifs);

		//insert it to surface
		insertSurfaceSegmentVertex(steinerVert, encrochedSeg);
	}
	newTets.clear();
	tobeDeletedTets.clear();
}

void DelMesher::insertVertex(Vertex *u, const Tetrahedron& tet, TetMeshDataStructure& meshRep,
	const VertexInsertionFlags& vifs, Tetrahedron *rt){
	Vertex *a = tet.v[0], *b = tet.v[1], *c = tet.v[2], *d = tet.v[3];

	if(!vifs.trueInsertion) tobeDeletedTets.push_back(tet);
	meshRep.deleteTetrahedron(a, b, c, d);

	if (vifs.encroachSegTest) {
		mayEncroachedSegsSet.clear();
		Segment ab(a, b, true), ac(a, c, true), ad(a, d, true), bc(b, c, true), bd(b, d, true), cd(c, d, true);
		if (findSegment(ab)) mayEncroachedSegsSet.insert(ab);
		if (findSegment(ac)) mayEncroachedSegsSet.insert(ac);
		if (findSegment(ad)) mayEncroachedSegsSet.insert(ad);
		if (findSegment(bc)) mayEncroachedSegsSet.insert(bc);
		if (findSegment(bd)) mayEncroachedSegsSet.insert(bd);
		if (findSegment(cd)) mayEncroachedSegsSet.insert(cd);
	}

	if (!vifs.cdt) {
		digCavity(u, Face(b, d, c, true), meshRep, vifs, rt);
		digCavity(u, Face(a, c, d, true), meshRep, vifs);
		digCavity(u, Face(a, d, b, true), meshRep, vifs);
		digCavity(u, Face(a, b, c, true), meshRep, vifs);
	}
	else {
		bool success = true;
		VertexInsertionFlags nvifs = vifs;
		nvifs.trueInsertion = false;

		success &= digCavity(u, Face(b, d, c, true), meshRep, nvifs);
		success &= digCavity(u, Face(a, c, d, true), meshRep, nvifs);
		success &= digCavity(u, Face(a, d, b, true), meshRep, nvifs);
		success &= digCavity(u, Face(a, b, c, true), meshRep, nvifs);

		if (success) {
			if (vifs.trueInsertion) {
				for (auto t : newTets) meshRep.addTetrahedron(t.v[0], t.v[1], t.v[2], t.v[3]);
				newTets.clear();
			}
		}
		else {
			std::vector<Vertex *> cavityVertices;
			cavityVertices.reserve(tobeDeletedTets.size());
			std::vector<Face> boundaryFaces;
			boundaryFaces.reserve(tobeDeletedTets.size());

			cavityVertices.push_back(u);
			u->setMark();
			for (auto t : newTets) {
				Face f(t.v[1], t.v[2], t.v[3]);
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
			triangulateCavity(std::vector<Face>(0), vifs.encroachFaceTest, boundaryFaces, cavityVertices, tobeDeletedTets, newTets, Face());
			for (auto v : cavityVertices) v->unSetMark();

			if (vifs.encroachSegTest) {
				int finalCount = tobeDeletedTets.size();
				for (int i = initDeletedCount; i < finalCount; i++) {
					Tetrahedron t = tobeDeletedTets[i];
					Segment ab(t.v[0], t.v[1], true), ac(t.v[0], t.v[2], true), ad(t.v[0], t.v[3], true),
						bc(t.v[1], t.v[2], true), bd(t.v[1], t.v[3], true), cd(t.v[2], t.v[3], true);
					if (findSegment(ab)) mayEncroachedSegsSet.insert(ab);
					if (findSegment(ac)) mayEncroachedSegsSet.insert(ac);
					if (findSegment(ad)) mayEncroachedSegsSet.insert(ad);
					if (findSegment(bc)) mayEncroachedSegsSet.insert(bc);
					if (findSegment(bd)) mayEncroachedSegsSet.insert(bd);
					if (findSegment(cd)) mayEncroachedSegsSet.insert(cd);
				}
			}

			if (vifs.trueInsertion) {
				for (auto t : newTets) meshRep.addTetrahedron(t.v[0], t.v[1], t.v[2], t.v[3]);
				newTets.clear();
			}
		}

		if (rt)
			meshRep.adjacent2Vertex(u, rt);

		tobeDeletedTets.clear();
	}

	if (vifs.encroachSegTest)
		for (auto s : mayEncroachedSegsSet) mayEncroachedSegs.push_back(s);
}

void DelMesher::insertSurfaceVertex(Vertex *u, const Face& f, bool insertToQueue){
	Vertex *a = f.v[0], *b = f.v[1], *c = f.v[2];
	DelVector aboveVert = Geometer::calculateAbovePoint(a->vert, b->vert, c->vert);
	surfaceRep.deleteTriangle(a, b, c);
	digCavity(u, aboveVert, Segment(b, a), insertToQueue);
	digCavity(u, aboveVert, Segment(c, b), insertToQueue);
	digCavity(u, aboveVert, Segment(a, c), insertToQueue);
}

void DelMesher::insertSurfaceSegmentVertex(Vertex *u, const Segment &s, Face *inFace, bool insertToQueue) {
	Vertex *a = s.v[0], *b = s.v[1];
	bool onFacet = false;
	Vertex *oppo = NULL;
	if (surfaceRep.Adjacent(s, &oppo)) {
		onFacet = true;
		surfaceRep.deleteTriangle(a, b, oppo);
		DelVector above = Geometer::calculateAbovePoint(a->vert, b->vert, oppo->vert);
		digCavity(u, above, Segment(a, oppo), insertToQueue);
		digCavity(u, above, Segment(oppo, b), insertToQueue);

		if (inFace) *inFace = Face(b, a, oppo);
	}
	if (surfaceRep.Adjacent(Segment(b, a), &oppo)) {
		onFacet = true;
		DelVector above = Geometer::calculateAbovePoint(b->vert, a->vert, oppo->vert);
		digCavity(u, above, Segment(b, oppo), insertToQueue);
		digCavity(u, above, Segment(oppo, a), insertToQueue);

		if (inFace) *inFace = Face(a, b, oppo);
	}

	if (onFacet) u->type = VertexType(u->type | VertexType::Vertex_Facet);
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
		if (findSegment(Segment(u, v, true))) boundaries.push_back(Segment(v, u));
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

bool DelMesher::Encroached(const Segment& s) const{
	Segment s2(s.v[1], s.v[0]);
	Vertex *x, *y;
	if (!surfaceRep.Adjacent(s, &x) || !surfaceRep.Adjacent(s2, &y))
		return true;

	DelVector center;
	REAL r = 0;
	Geometer::Orthocenter(s.v[0]->vert, s.v[0]->weight, s.v[1]->vert, s.v[1]->weight, &center, &r);
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
	if (!meshRep.Adjacent(f, &x) || !meshRep.Adjacent(f2, &y))
		return true;
	Vertex c;
	REAL r = 0;
	Geometer::Orthocircle(f.v[0]->vert, f.v[0]->weight, f.v[1]->vert, f.v[1]->weight, f.v[2]->vert, f.v[2]->weight, &c.vert, &r);
	REAL r2 = r*r;

	bool encroached0 = false;
	bool encroached1 = false;
	if (!x->isGhost())
		encroached0 = ((x->vert - c.vert).length2() - x->weight) <= r2;
	if (!y->isGhost())
		encroached1 = ((y->vert - c.vert).length2() - y->weight) <= r2;

	if (encroached0 || encroached1)
		return true;
	return false;
}

bool DelMesher::Encroached(const Segment &s, const Vertex &v) const{
	DelVector center;
	REAL r = 0;
	Geometer::Orthocenter(s.v[0]->vert, s.v[0]->weight, s.v[1]->vert, s.v[1]->weight, &center, &r);

	if ((v.vert - center).length2() - v.weight <= r*r)
		return true;

	return false;
}

bool DelMesher::Encroached(const Face &f, const Vertex &v) const{
	Vertex center;
	REAL r;
	Geometer::Orthocircle(f.v[0]->vert, f.v[0]->weight, f.v[1]->vert, f.v[1]->weight, f.v[2]->vert, f.v[2]->weight, &center.vert, &r);
	REAL r2 = r*r;

	return ((v.vert - center.vert).length2() - v.weight) <= r2;
}

Reference<Mesh> DelMesher::generateMesh(int *vertexLableMap){
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
		if (!surfaceRep.Contain(f) || !Encroached(f)){
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
		if (!meshRep.Contain(t)){
			skinnyTets.pop();
			continue;
		}
		splitTetrahedron(t);
		if (!mayEncroachedSegs.empty() || !mayEncroachedFaces.empty())
			goto start;
		skinnyTets.pop();
	}

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
