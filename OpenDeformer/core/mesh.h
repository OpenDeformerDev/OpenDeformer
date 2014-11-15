#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_CORE_MESH_H
#define ODER_CORE_MESH_H

#include "oder.h"
#include "memory.h"
#include "latool.h"
#include "sparseMatrix.h"
#include "mechMaterial.h"

namespace ODER{
	class Mesh : public ReferenceCounted{
	public:
		Mesh() :numNodes(0), numElements(0), numSurfaces(0), numNodesPerElement(0), numVertPerSur(0)
		        ,vertices(NULL), elements(NULL), surfaces(NULL){}
		Mesh(int numN, int numE, int numS, int numNPE, int numVPS = 3);
		Mesh(const Mesh& m) = delete;
		Mesh& operator=(const Mesh& m) = delete;

		virtual Element* getEmptyElement() const = 0;
		virtual Element* getEmptyMaterialElement(MarterialType type) const = 0;
		virtual Facet* getEmptyFacet() const = 0;
		virtual Element* getElement(int elementIndex, MarterialType type) const = 0;

		virtual int getCloestNode(const Vector &v)  const { return -1; }
		virtual int getCloestElement(const Vector &v) const { return -1; }

		int getNodeCount() const { return numNodes; }
		int getNodePerElementCount() const{ return numNodesPerElement; }
		int getElementCount() const{ return numElements; }
		int getVertPerFacetCount() const { return numVertPerSur; }
		int getFacetCount() const{ return numSurfaces; }
		Vector getVertex(int vertIndex) const{ return vertices[vertIndex]; }
		void setVertex(int vertIndex, const Vector& vert){ vertices[vertIndex] = vert; }
		int getElementNodeIndex(int elementIndex, int nodeIndex){ return elements[elementIndex*numNodesPerElement + nodeIndex]; }
		const int* getElementNodeReference(int elementIndex) const{ return &elements[elementIndex*numNodesPerElement]; }
		const int* getFacetVertReference(int facetIndex) const{ return &surfaces[facetIndex*numVertPerSur]; }
		void setElement(int elementIndex, int *nodeIndices);
		void setFacet(int facetIndex, int *facetIndices);

		virtual ~Mesh();

	private:
		int numNodes;
		int numElements;
		int numSurfaces;
		const int numNodesPerElement;
		const int numVertPerSur;

		Vector *vertices;
		int *elements;
		int *surfaces;
	};
}

#endif
