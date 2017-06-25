#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_CORE_MECHMATERIAL_H
#define ODER_CORE_MECHMATERIAL_H

#include "oder.h"
#include "memory.h"

namespace ODER{
	enum MarterialType{
		Marterial_Null = 0,
		Marterial_Isotropic = 1 << 0,
		Marterial_Anisortropic = 1 << 1,
		Marterial_Linear = 1 << 2,
		Marterial_NonLinear = 1 << 3,
		Marterial_Reduced = 1 << 4,
		Marterial_Invertible = 1 << 5,
		Marterial_Corotational = 1 << 6,
		Marterial_Plastic = 1 << 7
	};

	class SparseSymMatrixIndicesPerElementCache;

	class MechMaterial{
	public:
		MechMaterial(Scalar rho, MarterialType t) :density(rho), type(t){}
		void generateMassMatrix(const Reference<Mesh> &mesh, const Reference<NodeIndexer> &indexer, SparseMatrixAssembler& matrix) const;
		void generateMassMatrix(const Reference<Mesh> &mesh, const Reference<NodeIndexer> &indexer, BlockedSymSpMatrixAssembler& matrix) const;
		void generateMassMatrix(const Reference<Mesh> &mesh, const Reference<NodeIndexer> &indexer, 
			const SparseSymMatrixIndicesPerElementCache& matrixIndices, BlockedSymSpMatrix& matrix) const;

		Scalar getDensity() const{ return density; }
		MarterialType getMaterialType() const{ return type; }
		
		virtual ~MechMaterial() = default;
	protected:
		Scalar density;
		MarterialType type;
	};

	inline bool matchMaterialFlag(MarterialType type, MarterialType flag){
		return (type & flag) == flag;
	}
}

#endif



