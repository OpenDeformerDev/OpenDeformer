#include "stdafx.h"
#include "DelMesher.h"
#include <memory>
#include <set>
#include <random>

namespace ODER{

Predicator<REAL> DelTriangulator::predicator;

void DelTriangulator::generateSubPolygons(Vertex **vertices, int *segments, int vertexCount, int segmentCount, bool boundaryOnly){
	meshRep.Clear();
	this->segments.clear();
	invertion = false;

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

	calculateAbovePoint(vertexCount, vertices);
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

	invertion = detectInversion(f);
	//clean out of region triangles
	if(boundaryOnly) cleanRigion(f);
}

bool DelTriangulator::detectInversion(Face& ghostFace) const {
	Face f;
	meshRep.adjacent2Vertex(ghost, &f);
	Vertex *b = f.v[1], *c = f.v[2];
	bool done = false;
	bool flip = false;
	do {
		if (segments.find(Segment(c, b)) != segments.end()) done = true;
		else if (segments.find(Segment(b, c)) != segments.end()) {
			flip = true;
			done = true;
		}
		else {
			meshRep.Adjacent(Segment(ghost, c), &b);
			std::swap(b, c);
		}
	} while (!done && (b != f.v[1]));
	Assert(done);

	ghostFace = Face(ghost, b, c);
	return flip;
}

void DelTriangulator::cleanRigion(const Face& ghostFace) {
	Vertex *a = ghostFace.v[0], *b = ghostFace.v[1], *c = ghostFace.v[2];
	Assert(a == ghost);
	meshRep.deleteTriangle(a, b, c);
	propagateClean(Segment(b, a), 0);
	propagateClean(Segment(a, c), 0);

	for (auto s : segments) {
		if (!invertion) std::swap(s.v[0], s.v[1]);
		propagateClean(s, 0);
	}
}

void DelTriangulator::propagateClean(const Segment& s, int depth) {
	Vertex *w;
	if (!meshRep.Adjacent(s, &w)) return;

	Vertex *u = s.v[0], *v = s.v[1];

	meshRep.deleteTriangle(w, u, v);

	Segment uw(u, w), wv(w, v);
	Segment s0 = uw, s1 = wv;
	if (invertion) {
		s0.v[0] = w; s0.v[1] = u;
		s1.v[0] = v; s1.v[1] = w;
	}
	if (segments.find(s0) == segments.end()) propagateClean(uw, depth + 1);
	if (segments.find(s1) == segments.end()) propagateClean(wv, depth + 1);
}

void DelTriangulator::outPut(DelMesher *mesher){
	auto output = meshRep.getTriangles(false);
	for (auto f : output){
		if (invertion) std::swap(f.v[1], f.v[2]);
		(mesher->surfaceRep).addTriangle(f.v[0], f.v[1], f.v[2]);
		mesher->mayEncroachedFaces.push_back(f);
	}
}

Face DelTriangulator::findPosition(Vertex *u, const Face &f) const{
	size_t maxIter = 16384;
	Vertex *a = f.v[0], *b = f.v[1], *c = f.v[2];
	//ghost triangle
	if (a == ghost){
		if (predicator.inHalfSpace2D(u->vert, b->vert, c->vert, ghost->vert))
			return f;
		Vertex *x = ghost;
		meshRep.Adjacent(Segment(c, b), &x);
		//Face(x, c, b)
		a = x; std::swap(b, c);
	}

	Assert(a != ghost);
	//normal triangle

	while (maxIter--){
		//caution:because NULL->ghost,can't use NULL,
		//so it's might case memory corrupt above
		Vertex *x;
		bool ori0 = predicator.orient2d(b->vert, a->vert, u->vert, ghost->vert) > 0.0;
		bool ori1 = predicator.orient2d(a->vert, c->vert, u->vert, ghost->vert) > 0.0;

		int condition = ori0 + (ori1 << 1);
		switch (condition){
		case 0:
			if (predicator.orient2d(c->vert, b->vert, u->vert, ghost->vert) > 0.0){
				meshRep.Adjacent(Segment(c, b), &x);
				if (x == ghost)
					return Face(ghost, c, b);
				a = x; std::swap(b, c); //Face(x, c, b)
			}
			else
				return Face(a, b, c, true);
			break;
		case 1:
			meshRep.Adjacent(Segment(b, a), &x);
			if (x == ghost)
				return Face(ghost, b, a);
			c = a; a = x; //Face(x, b, a)
			break;
		case 2:
			meshRep.Adjacent(Segment(a, c), &x);
			if (x == ghost)
				return Face(ghost, a, c);
			b = a; a = x; //Face(x, a, c)
			break;
		case 3:
			if (Randomnation(2)){
				meshRep.Adjacent(Segment(b, a), &x);
				if (x == ghost)
					return Face(ghost, b, a);
				c = a; a = x; //Face(x, b, a)
			}
			else{
				meshRep.Adjacent(Segment(a, c), &x);
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

	bool find = meshRep.findIntersectedFace(a, bb, above, &f);
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
	auto triangles = cavityRep.getTriangles(false);
	for (auto f : triangles) {
		if(!meshRep.Contain(f)) 
			meshRep.addTriangle(f.v[0], f.v[1], f.v[2]);
	}

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
	int *prev = new int[size];
	int *next = new int[size];
	int *indices = new int[size];

	for (int i = 0; i < size; i++) {
		prev[i] = i;
		next[i] = i + 2;
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
			                               verts[0]->vert, verts[0]->weight, verts[2]->vert, verts[2]->weight, ghost->vert) > 0.0;
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

void DelTriangulator::calculateAbovePoint(int vertexCount, Vertex** vertices) {
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

	DelVector n = Normalize(vec % (c->vert - a->vert));
	ghost->vert = a->vert + (longest * REAL(0.5)) * n;
}

Predicator<REAL> DelMesher::predicator;

DelMesher::DelMesher(Vector *surfvs, int *segis, int *subpolygons, int numv, int numseg, int numpol, 
	int *numsubpol, REAL maxRa, REAL maxR) : oriVertexIndices(numv){
	maxRatio = maxRa;
	maxRadius = maxR;
	Tetrahedron::maxREration = maxRa;

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
		oriVertexIndices[i] = i;
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
		t.generateSubPolygons(&polygonVerts[0], NULL, numsubpol[i], 0, true);
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
		if (!predicator.fastCoPlane(a->vert, b->vert, c->vert, d->vert)){
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

	meshRep.addTetrahedron(a, b, c, d);
	meshRep.addTetrahedron(ghost, b, d, c);
	meshRep.addTetrahedron(ghost, a, c, d);
	meshRep.addTetrahedron(ghost, a, b, c);
	meshRep.addTetrahedron(ghost, a, d, b);

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

Tetrahedron DelMesher::findPosition(Vertex *u, const Tetrahedron &t, bool boundary) const{
	int maxIter = 16384;
	Vertex *a = t.v[0], *b = t.v[1], *c = t.v[2], *d = t.v[3];

	while (maxIter--){
		Vertex *x = NULL;
		if (a == ghost){
			if (predicator.inOrthoHalfSpace3D(u->vert, u->weight, b->vert, b->weight, c->vert, c->weight, d->vert, d->weight))
				return t;
			meshRep.Adjacent(Face(b, d, c, true), &x);
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
					meshRep.Adjacent(Face(d, c, b, true), &x);
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
				meshRep.Adjacent(Face(a, b, c, true), &x);
				if (x == ghost)
					return Tetrahedron(ghost, a, b, c, true);
				//Tetrahedron(x, a, b, c)
				d = c; c = b; b = a; a = x;
				break;
			case 2:
				meshRep.Adjacent(Face(a, c, d, true), &x);
				if (x == ghost)
					return Tetrahedron(ghost, a, c, d, true);
				//Tetrahedron(x, a, c, d)
				b = a; a = x;
				break;
			case 3:
				if (Randomnation(2) == 0){
					meshRep.Adjacent(Face(a, b, c, true), &x);
					if (x == ghost)
						return Tetrahedron(ghost, a, b, c, true);
					//Tetrahedron(x, a, b, c)
					d = c; c = b; b = a; a = x;
				}
				else{
					meshRep.Adjacent(Face(a, c, d, true), &x);
					if (x == ghost)
						return Tetrahedron(ghost, a, c, d, true);
					//Tetrahedron(x, a, c, d)
					b = a; a = x;
				}
				break;
			case 4:
				meshRep.Adjacent(Face(a, d, b, true), &x);
				if (x == ghost)
					return Tetrahedron(ghost, a, d, b, true);
				//Tetrahedron(x, a, d, b)
				c = d; d = b; b = a; a = x;
				break;
			case 5:
				if (Randomnation(2) == 0){
					meshRep.Adjacent(Face(a, b, c, true), &x);
					if (x == ghost)
						return Tetrahedron(ghost, a, b, c, true);
					//Tetrahedron(x, a, b, c)
					d = c; c = b; b = a; a = x;
				}
				else{
					meshRep.Adjacent(Face(a, d, b, true), &x);
					if (x == ghost)
						return Tetrahedron(ghost, a, d, b, true);
					//Tetrahedron(x, a, d, b)
					c = d; d = b; b = a; a = x;
				}
				break;
			case 6:
				if (Randomnation(2) == 0){
					meshRep.Adjacent(Face(a, c, d, true), &x);
					if (x == ghost)
						return Tetrahedron(ghost, a, c, d, true);
					//Tetrahedron(x, a, c, d)
					b = a; a = x;
				}
				else{
					meshRep.Adjacent(Face(a, d, b, true), &x);
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
					meshRep.Adjacent(Face(a, b, c, true), &x);
					if (x == ghost)
						return Tetrahedron(ghost, a, b, c, true);
					//Tetrahedron(x, a, b, c)
					d = c; c = b; b = a; a = x;
				}
				else if (term == 1){
					meshRep.Adjacent(Face(a, c, d, true), &x);
					if (x == ghost)
						return Tetrahedron(ghost, a, c, d, true);
					b = a; a = x;
				}
				else{
					meshRep.Adjacent(Face(a, d, b, true), &x);
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
		bool ori0 = predicator.orient2d(b->vert, a->vert, u->vert, c->vert + n) > 0;
		bool ori1 = predicator.orient2d(a->vert, c->vert, u->vert, c->vert + n) > 0;

		int condition = ori0 + (ori1 << 1);
		switch (condition){
		case 0:
			if (predicator.orient2d(c->vert, b->vert, u->vert, n) > 0){
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

void DelMesher::digCavity(Vertex *u, const Face& f, Tetrahedron *rt, bool insertToSkinny, bool trulyDeleteOrAdd){
	Vertex *a;
	if (!meshRep.Adjacent(f, &a))
		return;

	Vertex *b = f.v[0];
	Vertex *c = f.v[1];
	Vertex *d = f.v[2];

	bool deleted = false;
	if (a == ghost)
	    deleted = predicator.inOrthoHalfSpace3D(u->vert, u->weight, b->vert, b->weight, c->vert, c->weight, d->vert, d->weight);
	else if (b == ghost){
		deleted = predicator.inOrthoHalfSpace3D(u->vert, u->weight, a->vert, a->weight, d->vert, d->weight, c->vert, c->weight);
	}
	else
		deleted = predicator.inOrthoSpherePerturbed(u->vert, u->weight, b->vert, b->weight, d->vert, d->weight, c->vert, c->weight, a->vert, a->weight) > 0.0;

	if (deleted){
		if (!trulyDeleteOrAdd)
			tobeDeletedTets.push_back(Tetrahedron(a, b, c, d));
		meshRep.deleteTetrahedron(a, b, c, d);
		digCavity(u, Face(a, c, d, true), rt, insertToSkinny, trulyDeleteOrAdd);
		digCavity(u, Face(a, b, c, true), rt, insertToSkinny, trulyDeleteOrAdd);
		digCavity(u, Face(a, d, b, true), rt, insertToSkinny, trulyDeleteOrAdd);
	}
	else{
		if (trulyDeleteOrAdd){
			meshRep.addTetrahedron(u, b, d, c);
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
				if (surfaceRep.Contain(faceOut))
					mayEncroachedFaces.push_back(faceOut);
				if (surfaceRep.Contain(faceIn))
					mayEncroachedFaces.push_back(faceIn);
		}
	}
}

void DelMesher::digCavity(Vertex *u, const Vertex& aboveVert,const Segment &s, bool insertToQueue, bool trulyDeleteOrAdd){
	Vertex *x;
	if (!surfaceRep.Adjacent(s, &x))
		return;

	Vertex *w = s.v[0];
	Vertex *v = s.v[1];
	bool deleted = false;
	if (!findSegment(Segment(w, v, true)))
		deleted = predicator.inOrthoSpherePerturbed(aboveVert.vert, aboveVert.weight, u->vert, u->weight, v->vert, v->weight, w->vert, w->weight, x->vert, x->weight) > 0.0;

	if (deleted){
		if (!trulyDeleteOrAdd)
			tobeDeletedFaces.push_back(Face(w, v, x, true));
		surfaceRep.deleteTriangle(w, v, x);
		digCavity(u, aboveVert, Segment(x, v), insertToQueue, trulyDeleteOrAdd);
		digCavity(u, aboveVert, Segment(w, x), insertToQueue, trulyDeleteOrAdd);
	}
	else{
		if (trulyDeleteOrAdd)
			surfaceRep.addTriangle(u, v, w);
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

	if (!surfaceRep.adjacent2Vertex(center, &startFace)){
		surfaceRep.adjacent2Vertex(center, &startFace);
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
		if (predicator.fastCoPlane(center->vert, b->vert, c->vert, origin.vert)){
			ret = Face(center, b, c);
			break;
		}
		surfaceRep.Adjacent(Segment(center, c), &b);
		std::swap(b, c);
	} while (b != mark);
	return ret;
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
			nextMove = Move::Horizon;
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

Vertex*  DelMesher::findSegmentEncroachedReference(Vertex *end, const Tetrahedron& intersected) const {
	Vertex *a = intersected.v[0], *b = intersected.v[1], *c = intersected.v[2], *d = intersected.v[3];

	DelVector startVert = a->vert, endVert = end->vert;
	Vertex *ref = b;
	REAL maxAngle = interiorAngle(ref->vert, startVert, endVert);

	REAL angle = interiorAngle(c->vert, startVert, endVert);
	if (angle > maxAngle) {
		ref = c;
		maxAngle = angle;
	}
	angle = interiorAngle(d->vert, startVert, endVert);
	if (angle > maxAngle) {
		ref = d;
		maxAngle = angle;
	}

	enum Intersection { AcrossFace, AcrossEdge };
	Intersection intersect = Intersection::AcrossFace;
	if (predicator.orient3d(endVert, a->vert, d->vert, c->vert) == 0) {
		std::swap(d, b);
		std::swap(b, c);
		intersect = Intersection::AcrossEdge;
	}
	if (predicator.orient3d(endVert, a->vert, b->vert, d->vert) == 0) {
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
				angle = interiorAngle(a->vert, startVert, endVert);
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
					std::swap(c, d); std::swap(b, c);
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
				if (!a->isGhost()) {
					angle = interiorAngle(a->vert, startVert, endVert);
					if (angle > maxAngle) {
						angle = maxAngle;
						ref = a;
					}
					std::swap(a, d);
					bool dab = predicator.Intersection(d->vert, a->vert, b->vert, startVert, endVert);
					bool dca = predicator.Intersection(d->vert, c->vert, a->vert, startVert, endVert);
					bool coplane = predicator.orient3d(endVert, b->vert, c->vert, d->vert) == 0;
					int condition = dab + (dca << 1);
					switch (condition) {
					case 0:
						break;
					case 1:
					{
						if (coplane)
							intersect = Intersection::AcrossEdge;
						else
							intersect = Intersection::AcrossFace;
						std::swap(a, b); std::swap(c, a);
						hasIntersect = true;
						break;
					}
					case 2:
						if (coplane)
							intersect = Intersection::AcrossEdge;
						else
							intersect = Intersection::AcrossFace;
						std::swap(a, c); std::swap(b, a);
						hasIntersect = true;
						break;
					case 3:
						if (d == end) terminate = true;
						std::swap(a, b); std::swap(c, d);
						hasIntersect = true;
						break;
					default:
						break;
					}
				}
				else {
					bool found = meshRep.Adjacent(Face(a, b, c), &d);
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

void DelMesher::splitSubSegment(const Segment& s, Vertex*  ref, bool onFace){
	Vertex *a = s.v[0], *b = s.v[1];
	DelVector aa = a->vert, bb = b->vert;
	DelVector ab = bb - aa;
	REAL invabLen = REAL(1) / ab.length();
	//calulate steiner point
	auto oriSegFound = vorisHash.find(a);
	Assert(oriSegFound != vorisHash.end());
	int oriSegIndex = oriSegFound->second;
	DelVector steinerPoint;

	auto refOriSegFound = vorisHash.find(ref);
	if (refOriSegFound != vorisHash.end()) {
		int refOriSegIndex = refOriSegFound->second;

		Segment oriSeg = oriSegments[oriSegIndex];
		Segment refOriSeg = oriSegments[refOriSegIndex];

		Vertex *oriSegStart = oriSeg.v[0], *oriSegEnd = oriSeg.v[1];
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
		else
			steinerPoint = aa + ((ref->vert - aa).length() * invabLen) * ab;
	}
	else 
		steinerPoint = aa + ((ref->vert - aa).length() * invabLen) * ab;

	//check to see steiner point isn't too close to both end point
	REAL t = (steinerPoint - aa).length() * invabLen;
	if (t < REAL(0.2) || t > REAL(0.8))
		steinerPoint = REAL(0.5) * (aa + bb);

	Vertex *steinerVert = allocVertex(steinerPoint);

	//insert two sub segments
	Segment s0(a, steinerVert, true);
	Segment s1(steinerVert, b, true);
	segments.erase(s);
	segments.insert(s0);
	segments.insert(s1);
	mayEncroachedSegs.push_back(s0);
	mayEncroachedSegs.push_back(s1);
	vorisHash.insert(std::pair<Vertex *, int>(steinerVert, oriSegIndex));

	Tetrahedron toBeDeleted;
	meshRep.adjacent2Vertex(a, &toBeDeleted);
	toBeDeleted = findPosition(steinerVert, toBeDeleted);
	insertVertex(steinerVert, toBeDeleted, NULL, true);

	if (onFace) {
		Vertex *c = NULL, *d = NULL;
		if (!surfaceRep.Adjacent(s, &c) || !surfaceRep.Adjacent(Segment(b, a), &d))
			Severe("Unexpected error in DelMesher::splitSubSegment");

		Vertex aboveVert0(aa + triangleNormal(aa, bb, c->vert));
		Vertex aboveVert1(bb + triangleNormal(bb, aa, d->vert));
		surfaceRep.deleteTriangle(a, b, c);
		digCavity(steinerVert, aboveVert0, Segment(c, b));
		digCavity(steinerVert, aboveVert0, Segment(a, c));

		surfaceRep.deleteTriangle(b, a, d);
		digCavity(steinerVert, aboveVert1, Segment(d, a));
		digCavity(steinerVert, aboveVert1, Segment(b, d));
	}
}

void DelMesher::splitSubSegment(const Segment &s){
	Vertex *a = s.v[0], *b = s.v[1];
	Vertex *mid = Cover(s);

	Vertex *c = NULL, *d = NULL;
	//should be orient inward
	Face f0, f1;

	if (!surfaceRep.Adjacent(s, &c) || !surfaceRep.Adjacent(Segment(b, a), &d)){
		Severe("Unexpected error in DelMesher::splitSubSegment");
	}
	else{
		Vertex aboveVert0(a->vert + Normalize((a->vert - c->vert) % (b->vert - c->vert)));
		Vertex aboveVert1(b->vert + Normalize((b->vert - d->vert) % (a->vert - d->vert)));
		surfaceRep.deleteTriangle(a, b, c);
		digCavity(mid, aboveVert0, Segment(c, b));
		digCavity(mid, aboveVert0, Segment(a, c));

		surfaceRep.deleteTriangle(b, a, d);
		digCavity(mid, aboveVert1, Segment(d, a));
		digCavity(mid, aboveVert1, Segment(b, d));

		f0 = Face(a, b, c, true);
		f1 = Face(b, a, d, true);
	}

	Tetrahedron toBeDeleted;
	Vertex *e = NULL;
	if (meshRep.Adjacent(f0, &e))
		toBeDeleted = Tetrahedron(e, f0.v[0], f0.v[1], f0.v[2], true);
	else if (meshRep.Adjacent(f1, &e))
		toBeDeleted = Tetrahedron(e, f1.v[0], f1.v[1], f1.v[2], true);
	else{
		meshRep.adjacent2Vertex(s.v[0], &toBeDeleted);
		toBeDeleted = findPosition(mid, toBeDeleted);
	}
	insertVertex(mid, toBeDeleted, NULL, true);
}

void DelMesher::splitSubPolygon(const Face &f){
	Vertex *a = f.v[0], *b = f.v[1], *c = f.v[2];
	Vertex center; 
	predicator.Orthocircle(a->vert, a->weight, b->vert, b->weight, c->vert, c->weight, &center.vert);
	Vertex aboveVert(a->vert + Normalize((a->vert - c->vert) % (b->vert - c->vert)));

	tobeDeletedFaces.push_back(f);
	surfaceRep.deleteTriangle(a, b, c);

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
			surfaceRep.addTriangle(newPoint, segOfFace.v[0], segOfFace.v[1]);
			mayEncroachedFaces.push_back(Face(newPoint, segOfFace.v[0], segOfFace.v[1], true));
		}
		Tetrahedron toBeDeleted;
		Vertex *d;
		if (meshRep.Adjacent(f, &d))
			toBeDeleted = Tetrahedron(d, a, b, c, true);
		else{
			meshRep.adjacent2Vertex(a, &toBeDeleted);
			toBeDeleted = findPosition(newPoint, toBeDeleted);
		}
		insertVertex(newPoint, toBeDeleted, NULL, true);
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
	predicator.Orthosphere(a->vert, a->weight, b->vert, b->weight, 
		c->vert, c->weight, d->vert, d->weight, &center.vert);

	tobeDeletedTets.push_back(tet);
	meshRep.deleteTetrahedron(a, b, c, d);

	digCavity(&center, Face(b, d, c, true), NULL, false, false);
	digCavity(&center, Face(a, c, d, true), NULL, false, false);
	digCavity(&center, Face(a, d, b, true), NULL, false, false);
	digCavity(&center, Face(a, b, c, true), NULL, false, false);

	bool encrochment = false;
	Face encrocachedFace;
	bool result = false;
	for (auto nf : newFacesOfTets){
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
		Vertex *newVert = allocVertex(center);
		vertices.push_back(newVert);
		for (auto f : newFacesOfTets){
			meshRep.addTetrahedron(newVert, f.v[0], f.v[1], f.v[2]);
			Tetrahedron t(newVert, f.v[0], f.v[1], f.v[2], true);
			t.setRationAndRadius(ghost);
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
	newFacesOfTets.clear();
}

void DelMesher::insertVertex(Vertex *u, const Tetrahedron& tet, Tetrahedron *rt, bool insertToSkinny){
	Vertex *a = tet.v[0], *b = tet.v[1], *c = tet.v[2], *d = tet.v[3];
	meshRep.deleteTetrahedron(a, b, c, d);
	digCavity(u, Face(b, d, c, true), rt, insertToSkinny, true);
	digCavity(u, Face(a, c, d, true), NULL, insertToSkinny, true);
	digCavity(u, Face(a, d, b, true), NULL, insertToSkinny, true);
	digCavity(u, Face(a, b, c, true), NULL, insertToSkinny, true);
}

void DelMesher::insertSurfaceVertex(Vertex *u, const Face& f, bool insertToQueue){
	Vertex *a = f.v[0], *b = f.v[1], *c = f.v[2];
	Vertex aboveVert(a->vert + Normalize((a->vert - c->vert) % (b->vert - c->vert)));
	surfaceRep.deleteTriangle(a, b, c);
	digCavity(u, aboveVert, Segment(b, a), insertToQueue);
	digCavity(u, aboveVert, Segment(c, b), insertToQueue);
	digCavity(u, aboveVert, Segment(a, c), insertToQueue);
}


bool DelMesher::Encroached(const Segment& s) const{
	Segment s2(s.v[1], s.v[0]);
	Vertex *x, *y;
	if (!surfaceRep.Adjacent(s, &x) || !surfaceRep.Adjacent(s2, &y))
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
	if (!meshRep.Adjacent(f, &x) || !meshRep.Adjacent(f2, &y))
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
		if (predicator.fastCoLine(a, b, c)){
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

Reference<Mesh> DelMesher::generateMesh(int *vertexLableMap){
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
	auto tets = meshRep.getTetraherons(false);

	auto polygons = surfaceRep.getTriangles(false);

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

	if (vertexLableMap) {
		for (int i = 0; i < oriVertexIndices.size(); i++)
			vertexLableMap[oriVertexIndices[i]] = vertIndices[i];
	}

	freeAligned(vertIndices);

	return mesh;
}
}

