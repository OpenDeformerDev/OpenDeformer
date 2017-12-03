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
		        ,vertices(NULL), elements(NULL), surfaces(NULL), displacements(NULL), facetElementMap(NULL) {}
		Mesh(int numN, int numE, int numS, int numNPE, int numVPS = 3);
		Mesh(const Mesh& m) = delete;
		Mesh& operator=(const Mesh& m) = delete;
		Mesh(Mesh&& m);
		Mesh& operator=(Mesh&& m) = delete;

		virtual GeometricElement* getGeometricElement() const = 0;
		virtual Element* getMaterialElement(MarterialType type) const = 0;
		virtual Facet* getFacet() const = 0;

		virtual int getCloestNode(const Vector3f &v)  const { return -1; }
		virtual int getCloestElement(const Vector3f &v) const { return -1; }

		int getNodeCount() const { return numNodes; }
		int getNodePerElementCount() const{ return numNodesPerElement; }
		int getElementCount() const{ return numElements; }
		int getVertPerFacetCount() const { return numVertPerSur; }
		int getFacetCount() const{ return numSurfaces; }
		const Vector3& getVertex(int vertIndex) const{ return vertices[vertIndex]; }
		Vector3& getVertexDisplacement(int vertIndex) { return displacements[vertIndex]; }
		const Vector3& getVertexDisplacementConst(int vertIndex) const { return displacements[vertIndex]; }
		void setVertex(int vertIndex, const Vector3& vert){ vertices[vertIndex] = vert; }
		int getElementNodeIndex(int elementIndex, int nodeIndex){ return elements[elementIndex*numNodesPerElement + nodeIndex]; }
		const int* getElementNodeReference(int elementIndex) const{ return &elements[elementIndex*numNodesPerElement]; }
		const int* getFacetVertReference(int facetIndex) const{ return &surfaces[facetIndex*numVertPerSur]; }
		void setElement(int elementIndex, int *nodeIndices);
		void setFacet(int facetIndex, int *facetIndices);

		void allocFacetElementMap() { facetElementMap = new int[numSurfaces]; }
		void setFacetElementIndex(int facetIndex, int elementIndex) { facetElementMap[facetIndex] = elementIndex; }
		int getMappedElementIndex(int facetIndex) const { return facetElementMap[facetIndex]; }

		virtual ~Mesh();

	private:
		int numNodes;
		int numElements;
		int numSurfaces;
		const int numNodesPerElement;
		const int numVertPerSur;

		Vector3 *vertices;
		int *elements;
		int *surfaces;
		int *facetElementMap;
		Vector3 *displacements;
	};

	template <class MeshSubClass> void setMeshFacetElementMap(MeshSubClass& m) {
		static_assert(std::is_base_of<Mesh, MeshSubClass>::value, "MeshSubClass must be sub-class of ODER::Mesh for ODER::setMeshFacetElementMap");
		static_assert(false, "Unimplemented function for MeshSubClass");
	}
}

#endif
