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
		Element(Mesh *m) :mesh(m), nodeIndices(NULL){}

		void setNodeIndices(int elementIndex){ nodeIndices = mesh->getElementNodeReference(elementIndex); }
		int getNodeIndex(int localIndex) const{ return nodeIndices[localIndex]; }
		int getLocalMatrixIndex(int aNodeIndex, int aNodeAxis, int bNodeIndex, int bNodeAxis) const;
		int getLocalMatrixIndex(int aNodeDofIndex, int bNodeDofIndex) const;
		virtual ~Element() = default;

	protected:
		const int *nodeIndices;
		Reference<Mesh> mesh;
	};

	struct Facet{
		Facet(Mesh *m) :mesh(m), vertIndexs(NULL){}
		virtual void getSurfVirtualWorks(Scalar surfForce[3], Scalar *result) const = 0;
		void setVertIndexs(int vertIndex){ vertIndexs = mesh->getFacetVertReference(vertIndex); }
		int getVertIndex(int localIndex) const { return vertIndexs[localIndex]; }

		virtual ~Facet() = default;
	protected:
		const int *vertIndexs;
		Reference<Mesh> mesh;
	};

	struct GeometricElement : public Element{
		GeometricElement(Mesh *m) : Element(m){}
		virtual void generateSubMassMatrix(Scalar *result) const = 0;
		virtual void getBodyVirtualWorks(Scalar bodyForce[3], Scalar *result) const = 0;
		
		virtual ~GeometricElement() = default;
	};

	struct LinearIsotropicElement : public Element{
		LinearIsotropicElement(Mesh *m): Element(m){}
		virtual void setBMatrixs() = 0;
		virtual void generateSubStiffnessMatrix(int aNodeIndex, int bNodeIndex, const Scalar *D, Scalar *result) const = 0;
		virtual ~LinearIsotropicElement() = default;
	};

	struct ReducedIsotropicElement : public Element{
		ReducedIsotropicElement(Mesh *m): Element(m){}
		virtual void setBMatrixs() = 0;
		virtual void Intergration(const Scalar *C, Scalar *nlpart, Scalar *nnpart) const = 0;
		virtual void generateSubStiffnessMatrix(int aNodeIndex, int bNodeIndex, const Scalar *D, Scalar *result) const = 0;
		virtual ~ReducedIsotropicElement() = default;
	};

	struct LinearAnisortropicElement : public Element{
		LinearAnisortropicElement(Mesh *m): Element(m){}
		virtual void setBMatrixs() = 0;
		virtual void generateSubStiffnessMatrix(int aNodeIndex, int bNodeIndex, const Scalar *D, Scalar *result) const = 0;
		virtual ~LinearAnisortropicElement() = default;
	};

	struct InvertibleHyperelasticElement : public Element {
		InvertibleHyperelasticElement(Mesh *m) : Element(m) {}
		virtual void getPrecomputes(Scalar *drivates, Scalar *deforamtionGradients) const = 0;
		virtual int getDirvateEntryCount() const = 0;
		virtual int getDeformGradientsPreEntryCount() const = 0;
		virtual int getQuadraturePointCount() const = 0;

		virtual void generateDeformationGradient(const Scalar *precompute, Scalar *gradients) const = 0;
		virtual void generateSubStiffnessMatrix(const Scalar *drivates, const Scalar *diags, const Scalar *leftOrthoMats,
			const Scalar *rightOrthoMats, const Scalar *energyGradients, const Scalar *energyHassians, Scalar *result) const = 0;
		virtual void generateNodalVirtualWorks(const Scalar *precompute, const Scalar *stresses, Scalar *result) const = 0;
		virtual ~InvertibleHyperelasticElement() = default;
	};

	struct CorotationalHyperelasticElement : public Element {
		CorotationalHyperelasticElement(Mesh *m) : Element(m) {}
		virtual void getPrecomputes(const Scalar *D, Scalar *initSubStiffMat, Scalar *deforamtionGradients) const = 0;
		virtual int getInitSubStiffMatEntryCount() const = 0;
		virtual int getDeformGradientsPreEntryCount() const = 0;
		virtual int getQuadraturePointCount() const = 0;

		virtual void generateProperOrthoMats(const Scalar *deformationGradientPrecomputed, Scalar threshold, Scalar *properOrthoMats) const = 0;
		virtual void generateSubStiffnessMatrixNodalVirtualWorks(const Scalar *orthoMats, const Scalar *initStiffMat,
			Scalar *subStiffMat, Scalar *nodalVirtualWorks) const = 0;
	};

	struct CorotationalPlasticElement : public Element {
		CorotationalPlasticElement(Mesh *m) : Element(m) {};
		virtual void getPrecomputes(const Scalar *D, Scalar *initSubStiffMat, Scalar *deforamtionGradients, Scalar *drivates) const = 0;
		virtual int getInitSubStiffMatEntryCount() const = 0;
		virtual int getDeformGradientsPreEntryCount() const = 0;
		virtual int getDirvateEntryCount() const = 0;
		virtual int getQuadraturePointCount() const = 0;

		virtual void generateDecomposedDeformationGradient(const Scalar *deformationGradientPrecomputed, Scalar threshold, 
			Scalar *properOrthoparts, Scalar *factoredParts) const = 0;
		virtual void generateSubStiffnessMatrix(const Scalar *orthoMats, const Scalar *initStiffMat, Scalar *subStiffMat) const = 0;
		virtual void generateNodalVirtualWorks(const Scalar *precompute, const Scalar *stresses, Scalar *result) const = 0;
	};

	inline int Element::getLocalMatrixIndex(int aNodeDofIndex, int bNodeDofIndex) const {
		int numNodeElement = mesh->getNodePerElementCount();
		if (aNodeDofIndex > bNodeDofIndex) std::swap(aNodeDofIndex, bNodeDofIndex);
		int diagIndex = aNodeDofIndex * numNodeElement * 3 - (((aNodeDofIndex - 1) * aNodeDofIndex) / 2);
		return diagIndex + bNodeDofIndex - aNodeDofIndex;
	}

	inline int Element::getLocalMatrixIndex(int aNodeIndex, int aNodeAxis, int bNodeIndex, int bNodeAxis) const {
		return getLocalMatrixIndex(aNodeIndex * 3 + aNodeAxis, bNodeIndex * 3 + bNodeAxis);
	}

}

#endif