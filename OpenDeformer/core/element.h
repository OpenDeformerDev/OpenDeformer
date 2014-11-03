#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_CORE_ELEMENT_H
#define ODER_CORE_ELEMENT_H

#include "oder.h"
#include "mecMaterial.h"

namespace ODER{
	struct Element{
		Element(Mesh *m, int index, MarterialType t);
		Element(Mesh *m, MarterialType t);
		Element(Mesh *m) :mesh(m), nodeIndexs(NULL), type(Marterial_Null), BMatrixs(NULL){}

		virtual void setBMatrixs() = 0;
		virtual void generateSubStiffnessMatrix(int aNodeIndex, int bNodeIndex, const double *D, double *result) const = 0;
		virtual void generateSubMassMatrix(double *result) const = 0;
		virtual void getBodyVirtualWorks(double bodyForce[3], double *result) const = 0;
		virtual void Intergration(const double *C, double *nlpart, double *nnpart) const = 0;
		void setNodeIndexs(int elementIndex);
		int getNodeIndex(int localIndex) const;
		virtual ~Element();

	protected:
		const int *nodeIndexs;
		Reference<Mesh> mesh;
		MarterialType type;
		double *BMatrixs;
	};

	struct Facet{
		Facet(Mesh *m) :mesh(m), vertIndexs(NULL){}
		Facet(Mesh *m, int index);
		virtual void setVertIndexs(int vertIndex);
		virtual void getSurfVirtualWorks(double surfForce[3], double *result) const = 0;
		int getVertIndex(int localIndex) const;

	protected:
		const int *vertIndexs;
		Reference<Mesh> mesh;
	};
}

#endif