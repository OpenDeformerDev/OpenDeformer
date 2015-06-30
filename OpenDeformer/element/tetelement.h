#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_ELEMENT_TETELEMENT_H
#define ODER_ELEMENT_TETELEMENT_H

#include "oder.h"
#include "element.h"
#include "tetmesh.h"

namespace ODER{
	struct TetElement : public GeometricElement{
		TetElement(TetMesh *m) :GeometricElement(m){}

		void generateSubMassMatrix(double *result) const;
		void getBodyVirtualWorks(double bodyForce[3], double *result) const;
		~TetElement() = default;
	private:
		float getVolume() const;
	};

	struct TetFacet : public Facet{
		TetFacet(TetMesh *m) :Facet(m){}
		void getSurfVirtualWorks(double surfForce[3], double *result) const;
		~TetFacet() = default;
	private:
		float getArea() const;
	};


	struct LinearIsotropicTetElement : public LinearIsotropicElement{
		LinearIsotropicTetElement(TetMesh *m) : LinearIsotropicElement(m){}
		void setBMatrixs();
		void generateSubStiffnessMatrix(int aNodeIndex, int bNodeIndex, const double *D, double *result) const;
		~LinearIsotropicTetElement() = default;
	private:
		float getVolume() const;
		void getShapeFunctionDerivatives(double *dndx, double *dndy, double *dndz) const;

		double BMatrixs[12];
	};

	struct ReducedIsotropicTetElement : public ReducedIsotropicElement{
		ReducedIsotropicTetElement(TetMesh *m) : ReducedIsotropicElement(m){}
		void setBMatrixs();
		void Intergration(const double *C, double *nlpart, double *nnpart) const;
		void generateSubStiffnessMatrix(int aNodeIndex, int bNodeIndex, const double *D, double *result) const;
		~ReducedIsotropicTetElement() = default;
	private:
		float getVolume() const;
		void getShapeFunctionDerivatives(double *dndx, double *dndy, double *dndz) const;

		double BMatrixs[12];
	};

	struct LinearAnisortropicTetElement : public LinearAnisortropicElement{
		LinearAnisortropicTetElement(TetMesh *m) : LinearAnisortropicElement(m){}
		void setBMatrixs();
		void generateSubStiffnessMatrix(int aNodeIndex, int bNodeIndex, const double *D, double *result) const;
	private:
		float getVolume() const;
		void getShapeFunctionDerivatives(double *dndx, double *dndy, double *dndz) const;

		double BMatrixs[12];
	};

}

#endif