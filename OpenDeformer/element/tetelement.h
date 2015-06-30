#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_ELEMENT_TETELEMENT_H
#define ODER_ELEMENT_TETELEMENT_H

#include "oder.h"
#include "element.h"
#include "tetmesh.h"

namespace ODER{
	struct TetElement : public virtual Element{
		TetElement(TetMesh *m) :Element(m){}

		void generateSubMassMatrix(double *result) const;
		void getBodyVirtualWorks(double bodyForce[3], double *result) const;
		~TetElement() = default;
	protected:
		float getVolume() const;
		void getShapeFunctionDerivatives(double *dndx, double *dndy, double *dndz) const;
	};

	struct TetFacet : public Facet{
		TetFacet(TetMesh *m) :Facet(m){}
		void getSurfVirtualWorks(double surfForce[3], double *result) const;
		~TetFacet() = default;
	private:
		float getArea() const;
	};


	struct LinearIsotropicTetElement : public TetElement, public LinearIsotropicElement{
		LinearIsotropicTetElement(TetMesh *m) :Element(m), TetElement(m), LinearIsotropicElement(m, 12){}
		void setBMatrixs();
		void generateSubStiffnessMatrix(int aNodeIndex, int bNodeIndex, const double *D, double *result) const;
		~LinearIsotropicTetElement() = default;
	};

	struct ReducedIsotropicTetElement : public TetElement, public ReducedIsotropicElement{
		ReducedIsotropicTetElement(TetMesh *m) :Element(m), TetElement(m), ReducedIsotropicElement(m, 12){}
		void setBMatrixs();
		void Intergration(const double *C, double *nlpart, double *nnpart) const;
		void generateSubStiffnessMatrix(int aNodeIndex, int bNodeIndex, const double *D, double *result) const;
		~ReducedIsotropicTetElement() = default;
	};

	struct LinearAnisortropicTetElement : public TetElement, public LinearAnisortropicElement{
		LinearAnisortropicTetElement(TetMesh *m) :Element(m), TetElement(m), LinearAnisortropicElement(m, 12){}
		void setBMatrixs();
		void generateSubStiffnessMatrix(int aNodeIndex, int bNodeIndex, const double *D, double *result) const;
	};


}

#endif