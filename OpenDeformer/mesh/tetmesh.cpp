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
		else{
			Severe("Unimplemented features in TetMesh::getMaterialElement");
			return NULL;
		}
	}


	Facet* TetMesh::getFacet() const{
		return new TetFacet((TetMesh *)this);
	}
}
