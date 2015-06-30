#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_CORE_ELEMENT_H
#define ODER_CORE_ELEMENT_H

#include "oder.h"
#include "mechMaterial.h"
#include "mesh.h"

namespace ODER{
	struct Element{
		Element(Mesh *m) :mesh(m), nodeIndexs(NULL){}

		virtual void generateSubMassMatrix(double *result) const = 0;
		virtual void getBodyVirtualWorks(double bodyForce[3], double *result) const = 0;
		void setNodeIndexs(int elementIndex){ nodeIndexs = mesh->getElementNodeReference(elementIndex); }
		int getNodeIndex(int localIndex) const{ return nodeIndexs[localIndex]; }
		virtual ~Element() = default;

	protected:
		const int *nodeIndexs;
		Reference<Mesh> mesh;
	};

	struct Facet{
		Facet(Mesh *m) :mesh(m), vertIndexs(NULL){}
		virtual void getSurfVirtualWorks(double surfForce[3], double *result) const = 0;
		void setVertIndexs(int vertIndex){ vertIndexs = mesh->getFacetVertReference(vertIndex); }
		int getVertIndex(int localIndex) const { return vertIndexs[localIndex]; }

		virtual ~Facet() = default;
	protected:
		const int *vertIndexs;
		Reference<Mesh> mesh;
	};

	struct LinearIsotropicElement : public virtual Element{
		LinearIsotropicElement(Mesh *m, int entries);
		virtual void setBMatrixs() = 0;
		virtual void generateSubStiffnessMatrix(int aNodeIndex, int bNodeIndex, const double *D, double *result) const = 0;
		virtual ~LinearIsotropicElement();
	protected:
		double *BMatrixs;
	};

	struct ReducedIsotropicElement : public virtual Element{
		ReducedIsotropicElement(Mesh *m, int entries);
		virtual void setBMatrixs() = 0;
		virtual void Intergration(const double *C, double *nlpart, double *nnpart) const = 0;
		virtual void generateSubStiffnessMatrix(int aNodeIndex, int bNodeIndex, const double *D, double *result) const = 0;
		virtual ~ReducedIsotropicElement();
	protected:
		double *BMatrixs;
	};

	struct LinearAnisortropicElement : public virtual Element{
		LinearAnisortropicElement(Mesh *m, int entries);
		virtual void setBMatrixs() = 0;
		virtual void generateSubStiffnessMatrix(int aNodeIndex, int bNodeIndex, const double *D, double *result) const = 0;
		virtual ~LinearAnisortropicElement();
	protected:
		double *BMatrixs;
	};

}

#endif