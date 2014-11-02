#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_CORE_MECMATERIAL_H
#define ODER_CORE_MECMATERIAL_H

#include "oder.h"
#include "memory.h"

namespace ODER{
	enum MarterialType{
		Marterial_Null = 0,
		Marterial_Isotropic = 1 << 0,
		Marterial_Anisortropic = 1 << 1,
		Marterial_Linear = 1 << 2,
		Marterial_NonLinear = 1 << 3
	};

	class MecMaterial : public ReferenceCounted{
	public:
		MecMaterial(double rho, MarterialType t) :density(rho), type(t){}
		void generateMassMatrix(const Reference<Mesh> &mesh, const Reference<NodeIndexer> &indexer, SparseMatrixAssembler& matrix) const;
		virtual void generateStiffnessMatrix(const Reference<Mesh> &mesh, const Reference<NodeIndexer> &indexer, SparseMatrixAssembler& matrix) const = 0;
		virtual void getNodeForces(const Reference<Mesh> &mesh, const Reference<NodeIndexer> &indexer, int order, const double *ds, double *forces){
			Severe("MecMaterial::getNodeForces is not implemented");
		}
		virtual void preprocessWithReduction(const Reference<Mesh> &mesh, const NodeIndexer &indexer){
			Severe("MecMaterial::preprocessWithReduction is not implemented");
		}
		virtual int getNonlinearAsymptoticOrder(){
			return 0;
		}

		double getDensity() const{ return density; }
		MarterialType getMaterialType() const{ return type; }
		
		virtual ~MecMaterial(){}
	protected:
		double density;
		MarterialType type;
	};

	inline bool matchMaterialFlag(MarterialType type, MarterialType flag){
		return (type & flag) == flag;
	}
}

#endif



