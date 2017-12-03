#include "stdafx.h"
#include "tetmesh.h"
#include "tetelement.h"

namespace ODER{
	GeometricElement* TetMesh::getGeometricElement() const{
		return new TetElement((TetMesh *)this);
	}

	Element* TetMesh::getMaterialElement(MarterialType type) const{
		if (matchMaterialFlag(type, MarterialType(Marterial_Isotropic | Marterial_Linear)))
			return new LinearIsotropicTetElement((TetMesh *)this);
		else if (matchMaterialFlag(type, MarterialType(Marterial_Reduced | Marterial_NonLinear | Marterial_Isotropic)))
			return new ReducedIsotropicTetElement((TetMesh *)this);
		else if (matchMaterialFlag(type, MarterialType(Marterial_Anisortropic | Marterial_Linear)))
			return new LinearAnisortropicTetElement((TetMesh *)this);
		else if (matchMaterialFlag(type, MarterialType(Marterial_Invertible | Marterial_Isotropic | Marterial_NonLinear)))
			return new InvertibleHyperelasticTetElement((TetMesh *)this);
		else if (matchMaterialFlag(type, MarterialType(Marterial_Corotational | Marterial_Isotropic | Marterial_NonLinear | Marterial_Plastic)))
			return new CorotationalPlasticTetElement((TetMesh *)this);
		else if (matchMaterialFlag(type, MarterialType(Marterial_Corotational | Marterial_Isotropic | Marterial_NonLinear)))
			return new CorotationalHyperelasticTetElement((TetMesh *)this);
		else{
			Severe("Unimplemented features in TetMesh::getMaterialElement");
			return NULL;
		}
	}

	Element* TetMesh::getElement(int elementIndex, MarterialType type) const {
		Element *ret = getMaterialElement(type);
		ret->setNodeIndices(elementIndex);
		return ret;
	}

	Facet* TetMesh::getFacet() const{
		return new TetFacet((TetMesh *)this);
	}

	template<> void setMeshFacetElementMap(TetMesh& m) {
		int triangleCount = m.getFacetCount();
		m.allocFacetElementMap();

		std::tuple<int, int, int> triangle;
		std::map<std::tuple<int, int, int>, int> triIndexMap;

		for (int i = 0; i < triangleCount; i++) {
			const int *tri = m.getFacetVertReference(i);
			std::get<0>(triangle) = tri[0];
			std::get<1>(triangle) = tri[1];
			std::get<2>(triangle) = tri[2];

			//sort index
			if (std::get<0>(triangle) > std::get<1>(triangle)) std::swap(std::get<0>(triangle), std::get<1>(triangle));
			if (std::get<1>(triangle) > std::get<2>(triangle)) std::swap(std::get<1>(triangle), std::get<2>(triangle));
			if (std::get<0>(triangle) > std::get<1>(triangle)) std::swap(std::get<0>(triangle), std::get<1>(triangle));

			triIndexMap.insert(std::make_pair(triangle, i));
		}


		int tetCount = m.getElementCount();
		for (int i = 0; i < tetCount; i++) {
			const int *tet = m.getElementNodeReference(i);
			for (int j = 0; j < 4; j++) {
				std::get<0>(triangle) = tet[j % 4];
				std::get<1>(triangle) = tet[(j + 1) % 4];
				std::get<2>(triangle) = tet[(j + 2) % 4];

				//sort index
				if (std::get<0>(triangle) > std::get<1>(triangle)) std::swap(std::get<0>(triangle), std::get<1>(triangle));
				if (std::get<1>(triangle) > std::get<2>(triangle)) std::swap(std::get<1>(triangle), std::get<2>(triangle));
				if (std::get<0>(triangle) > std::get<1>(triangle)) std::swap(std::get<0>(triangle), std::get<1>(triangle));

				auto found = triIndexMap.find(triangle);
				if (found != triIndexMap.end()) m.setFacetElementIndex(found->second, i);
			}
		}
	}
}
