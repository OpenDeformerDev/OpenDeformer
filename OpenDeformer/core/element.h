#pragma once
#include "oder.h"
#include "mecMaterial.h"

namespace ODER{
	struct Element{
		Element(Mesh *m, int index, MarterialType t);
		Element(Mesh *m, MarterialType t);
		Element(Mesh *m) :mesh(m), nodeIndexs(NULL), volume(-1.f), type(Marterial_Null), BMatrixs(NULL){}

		virtual void setVolume() = 0;
		virtual void setBMatrixs() = 0;
		virtual void generateSubStiffnessMatrix(int aNodeIndex, int bNodeIndex, const double *D, double *result) const = 0;
		virtual void getBodyVirtualWorks(double bodyForce[3], double *result) const = 0;
		virtual void setNodeIndexs(int elementIndex);
		~Element();

		int *nodeIndexs;
		float volume;

	protected:
		Reference<Mesh> mesh;
		MarterialType type;
		double *BMatrixs;
	};

	struct Facet{
		Facet(Mesh *m) :mesh(m), vertIndexs(NULL){}
		Facet(Mesh *m, int index);
		virtual float getArea() const { return -1.0; };
		virtual void setVertIndexs(int vertIndex);
		virtual void getSurfVirtualWorks(double surfForce[3], double *result) const = 0;

		int *vertIndexs;

	protected:
		Reference<Mesh> mesh;
	};
}