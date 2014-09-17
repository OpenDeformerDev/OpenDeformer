#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_MESH_TETMESH_H
#define ODER_MESH_TETMESH_H

#include "oder.h"
#include "mesh.h"

namespace ODER{
	class TetMesh : public Mesh{
	public:
		TetMesh(int numN, int numE, int numS) : Mesh(numN, numE, numS, 4, 3){}
		void generateMassMatrix(const MecMaterial &material, const NodeIndexer &indexer, SparseMatrixAssembler& m) const;
		Element* getEmptyElement() const;
		Element* getEmptyMaterialElement(MarterialType type) const;
		Facet* getEmptyFacet() const;
		Element* getElement(int elementIndex, MarterialType type) const;
		int getCloestNode(const Vector &v)  const { return -1; }
		int getCloestElement(const Vector &v) const { return -1; }
	};
}

#endif

