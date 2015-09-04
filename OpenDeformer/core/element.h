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

	struct GeometricElement : public Element{
		GeometricElement(Mesh *m) : Element(m){}
		virtual void generateSubMassMatrix(double *result) const = 0;
		virtual void getBodyVirtualWorks(double bodyForce[3], double *result) const = 0;
		
		virtual ~GeometricElement() = default;
	};

	struct LinearIsotropicElement : public Element{
		LinearIsotropicElement(Mesh *m): Element(m){}
		virtual void setBMatrixs() = 0;
		virtual void generateSubStiffnessMatrix(int aNodeIndex, int bNodeIndex, const double *D, double *result) const = 0;
		virtual ~LinearIsotropicElement() = default;
	};

	struct ReducedIsotropicElement : public Element{
		ReducedIsotropicElement(Mesh *m): Element(m){}
		virtual void setBMatrixs() = 0;
		virtual void Intergration(const double *C, double *nlpart, double *nnpart) const = 0;
		virtual void generateSubStiffnessMatrix(int aNodeIndex, int bNodeIndex, const double *D, double *result) const = 0;
		virtual ~ReducedIsotropicElement() = default;
	};

	struct LinearAnisortropicElement : public Element{
		LinearAnisortropicElement(Mesh *m): Element(m){}
		virtual void setBMatrixs() = 0;
		virtual void generateSubStiffnessMatrix(int aNodeIndex, int bNodeIndex, const double *D, double *result) const = 0;
		virtual ~LinearAnisortropicElement() = default;
	};

	struct InvertibleHyperelasticElement : public Element {
		InvertibleHyperelasticElement(Mesh *m) : Element(m) {}
		virtual void getPrecomputes(double *drivates, double *deforamtionGradients, double *virtualWorks) const = 0;
		virtual int getDirvateEntryCount() const = 0;
		virtual int getNodalVirtualWorksPreEntryCount() const = 0;
		virtual int getDeformGradientsPreEntryCount() const = 0;

		virtual void generateDeformationGradient(const double *precompute, double *gradients) const = 0;
		virtual void generateSubStiffnessMatrix(const double *drivates, const double *diags, const double *hassians, double *result) const = 0;
		virtual void generateNodalVirtualWorks(const double *precompute, double *result) const = 0;
		virtual ~InvertibleHyperelasticElement() = default;
	};

}

#endif