#include "stdafx.h"
#include "oder.h"
#include "latool.h"
#include "memory.h"
#include "DelMesher.h"
#include "InitShader.h"
#include "sparseMatrix.h"
#include "hooke.h"
#include "forcer.h"
#include "newmark.h"
#include "constrainer.h"
#include "nodeIndexer.h"

using namespace ODER;

Reference<Mesh> meshGenerateTester(float minR){
	Vector vertices[8] =
	{ Vector(0.5f, 0.1f, 0.1f), Vector(-0.5f, 0.1f, 0.1f),
	Vector(-0.5f, 0.1f, -0.1f), Vector(0.5f, 0.1f, -0.1f),
	Vector(0.5f, -0.1f, 0.1f), Vector(-0.5f, -0.1f, 0.1f),
	Vector(-0.5f, -0.1f, -0.1f), Vector(0.5f, -0.1f, -0.1f)
	};

	int segments[24] = { 0, 1, 1, 5, 5, 4, 4, 0, 7, 3, 3, 2, 2, 6, 6, 7, 1, 2, 0, 3, 4, 7, 5, 6 };

	int polygons[24] = {
		0, 3, 2, 1,
		0, 1, 5, 4,
		4, 5, 6, 7,
		2, 3, 7, 6,
		1, 2, 6, 5,
		0, 4, 7, 3
	};

	int polynum[6] = { 4, 4, 4, 4, 4, 4 };

	DelMesher m = DelMesher(vertices, segments, polygons, 8, 12, 6, polynum, 2.0f, minR);
	return m.generateMesh();
}

Reference<Mesh> m;
Mat4 modalView;
GLuint modalViewLoc;
GLuint vao, vertexBuffer, indiceBuffer;
Vector *verts = NULL;
double *displacement = NULL;
NodeIndexer indexer;
Intergrator *intergrator = NULL;

void init(){
	m = meshGenerateTester(0.03f);
	Constrainer constrainer(m);
	constrainer.setConstrains(Vector(-0.5f, 0.1f, 0.1f), Vector(-1.f, 0.f, 0.f));
	indexer = NodeIndexer(constrainer);
	HookeMaterial material(1000.0, 610000.0, 0.490);
	SparseMatrixAssembler Kassembler(indexer.getMatrixOrder(m));
	SparseMatrixAssembler Massembler(indexer.getMatrixOrder(m));
	m->generateMassMatrix(material, indexer, Massembler);
	material.generateStiffnessMatrix(m, indexer, Kassembler);

	double force[3] = { 0.0, -15.0, 0.0 };
	int nodes[4] = { 0, 3, 4, 7 };

	Forcer forcer(m);
	forcer.addSurfaceForceByNode(force, 4, nodes, indexer);

	SparseMatrix M(Massembler);
	SparseMatrix K(Kassembler);

	int total = M.getNumColumns();
	intergrator = new NewmarkIntergrator(total, 0.25, 0.5, M, K, 300, 0.1, 0.01, 0.001);
	intergrator->setExternalVirtualWork(forcer);

	verts = new Vector[m->numNodes];
	displacement = new double[total];

	const Vector eye = Vector(0.2f, 0.2f, 0.5f);
	const Vector at = Vector(0.f, 0.f, 0.f);
	const Vector up = Vector(0.f, 1.f, 0.f);

	modalView = lookAt(eye, at, up);
	modalView = Orient(-1.f, 1.f, -1.f, 1.f, 0.01f, 100.f) * modalView;

	int *indices = new int[12 * m->numElements];
	for (int i = 0; i < m->numElements; i++){
		int a = m->elements[4 * i];
		int b = m->elements[4 * i + 1];
		int c = m->elements[4 * i + 2];
		int d = m->elements[4 * i + 3];

		indices[12 * i] = a;
		indices[12 * i + 1] = c;
		indices[12 * i + 2] = d;

		indices[12 * i + 3] = a;
		indices[12 * i + 4] = b;
		indices[12 * i + 5] = c;

		indices[12 * i + 6] = a;
		indices[12 * i + 7] = d;
		indices[12 * i + 8] = b;

		indices[12 * i + 9] = b;
		indices[12 * i + 10] = d;
		indices[12 * i + 11] = c;
	}

	glEnable(GL_DEPTH_TEST);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

	glGenVertexArrays(1, &vao);
	glBindVertexArray(vao);

	glGenBuffers(1, &vertexBuffer);
	glBindBuffer(GL_ARRAY_BUFFER, vertexBuffer);
	glBufferData(GL_ARRAY_BUFFER, m->numNodes*sizeof(Vector), verts, GL_DYNAMIC_DRAW);

	glGenBuffers(1, &indiceBuffer);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indiceBuffer);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, 12 * m->numElements*sizeof(int), indices, GL_STATIC_DRAW);

	GLuint shader = initShader("shader/test.vs", "shader/test.fs");
	glUseProgram(shader);

	GLuint vPosition = glGetAttribLocation(shader, "vPosition");
	glEnableVertexAttribArray(vPosition);
	glVertexAttribPointer(vPosition, 3, GL_FLOAT, GL_FALSE, 0, GL_BUFFER_OFFSET(0));

	modalViewLoc = glGetUniformLocation(shader, "ModalView");

	glClearColor(1.0, 1.0, 1.0, 1.0);

	delete[] indices;

}


void display(){
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glUniformMatrix4fv(modalViewLoc, 1, GL_TRUE, &(modalView.m[0][0]));
	glBindBuffer(GL_ARRAY_BUFFER, vertexBuffer);
	glBufferData(GL_ARRAY_BUFFER, m->numNodes*sizeof(Vector), verts, GL_DYNAMIC_DRAW);

	glDrawElements(GL_TRIANGLES, 12 * m->numElements, GL_UNSIGNED_INT, 0);
	glutSwapBuffers();
}

void idle(){
	intergrator->runOneTimeStep();
	intergrator->getVertexPositions(m, indexer, m->numNodes, verts, displacement);
	glutPostRedisplay();
}

int main(int argc, char** argv){
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
	glutInitWindowSize(800, 600);
	glutInitContextVersion(4, 3);
	glutInitContextProfile(GLUT_CORE_PROFILE);
	glutCreateWindow("Test");

	glewInit();

	init();

	glutDisplayFunc(display);
	glutIdleFunc(idle);
	glutMainLoop();
	return 0;
}
