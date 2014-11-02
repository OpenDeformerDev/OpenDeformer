#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_ELEMENT_TETELEMENT_H
#define ODER_ELEMENT_TETELEMENT_H

#include "oder.h"
#include "element.h"
#include "tetmesh.h"

namespace ODER{
	struct TetElement : public Element{
		TetElement(TetMesh *m, int n, MarterialType t);
		TetElement(TetMesh *m, MarterialType t);
		TetElement(TetMesh *m) :Element(m){}
		void setBMatrixs();
		void generateSubStiffnessMatrix(int aNodeIndex, int bNodeIndex, const double *D, double *result) const;
		void generateSubMassMatrix(double *result) const;
		void Intergration(const double *C, double *nlpart, double *nnpart) const;
		void getBodyVirtualWorks(double bodyForce[3], double *result) const;
	private:
		float getVolume() const;
		void getShapeFunctionDerivatives(double *dndx, double *dndy, double *dndz) const;
	};

	struct TetFacet : public Facet{
		TetFacet(TetMesh *m) :Facet(m){}
		TetFacet(TetMesh *m, int index) :Facet(m, index){}
		void getSurfVirtualWorks(double surfForce[3], double *result) const;
	private:
		float getArea() const;
	};
}

#endif