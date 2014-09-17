#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_CORE_MESH_H
#define ODER_CORE_MESH_H

#include "oder.h"
#include "memory.h"
#include "latool.h"
#include "sparseMatrix.h"
#include "mecMaterial.h"

namespace ODER{
	class Mesh : public ReferenceCounted{
	public:
		Mesh() :numNodes(0), numElements(0), numSurfaces(0), numNodesPerElement(0), numVertPerSur(0){}
		Mesh(int numN, int numE, int numS, int numNPE, int numVPS = 3);
		virtual void generateMassMatrix(const MecMaterial &material, const NodeIndexer &indexer, SparseMatrixAssembler &m) const = 0;
		virtual Element* getEmptyElement() const = 0;
		virtual Element* getEmptyMaterialElement(MarterialType type) const = 0;
		virtual Facet* getEmptyFacet() const = 0;
		virtual Element* getElement(int elementIndex, MarterialType type) const = 0;
		void updateElementInfo(int elementIndex, Element *element) const;
		void updateFacetInfo(int elementIndex, Facet *facet) const;

		virtual int getCloestNode(const Vector &v)  const { return -1; }
		virtual int getCloestElement(const Vector &v) const { return -1; }
		~Mesh();

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
