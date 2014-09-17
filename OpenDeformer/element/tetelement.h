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
		void setVolume();
		void setBMatrixs();
		void generateSubStiffnessMatrix(int aNodeIndex, int bNodeIndex, const double *D, double *result) const;
		void getBodyVirtualWorks(double bodyForce[3], double *result) const;
	private:
		void getShapeFunctionDerivatives(double *dndx, double *dndy, double *dndz) const;
	};

	struct TetFacet : public Facet{
		TetFacet(TetMesh *m) :Facet(m){}
		TetFacet(TetMesh *m, int index) :Facet(m, index){}
		float getArea() const;
		void getSurfVirtualWorks(double surfForce[3], double *result) const;
	};
}

#endif