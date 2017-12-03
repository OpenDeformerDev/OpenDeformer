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
		TetMesh(const TetMesh&) = delete;
		TetMesh(TetMesh&& mesh) : Mesh(std::move(mesh)){}
		TetMesh& operator=(const TetMesh&) = delete;
		TetMesh& operator=(TetMesh&& mesh) = delete;
		GeometricElement* getGeometricElement() const;
		Element* getMaterialElement(MarterialType type) const;
		Facet* getFacet() const;
		Element* getElement(int elementIndex, MarterialType type) const;
		int getCloestNode(const Vector3f &v)  const { return -1; }
		int getCloestElement(const Vector3f &v) const { return -1; }
		~TetMesh() = default;
	};

	template<> void setMeshFacetElementMap(TetMesh& m);
}

#endif

