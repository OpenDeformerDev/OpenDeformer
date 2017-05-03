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

		void generateSubMassMatrix(Scalar *result) const;
		void getBodyVirtualWorks(Scalar bodyForce[3], Scalar *result) const;
		~TetElement() = default; 
	};

	struct TetFacet : public Facet{
		TetFacet(TetMesh *m) :Facet(m){}
		void getSurfVirtualWorks(Scalar surfForce[3], Scalar *result) const;
		~TetFacet() = default;
	};


	struct LinearIsotropicTetElement : public LinearIsotropicElement{
		LinearIsotropicTetElement(TetMesh *m) : LinearIsotropicElement(m){}
		void setBMatrixs();
		void generateSubStiffnessMatrix(int aNodeIndex, int bNodeIndex, const Scalar *D, Scalar *result) const;
		~LinearIsotropicTetElement() = default;
	private:
		Scalar BMatrixs[12];
	};

	struct ReducedIsotropicTetElement : public ReducedIsotropicElement{
		ReducedIsotropicTetElement(TetMesh *m) : ReducedIsotropicElement(m){}
		void setBMatrixs();
		void Intergration(const Scalar *C, Scalar *nlpart, Scalar *nnpart) const;
		void generateSubStiffnessMatrix(int aNodeIndex, int bNodeIndex, const Scalar *D, Scalar *result) const;
		~ReducedIsotropicTetElement() = default;
	private:
		Scalar BMatrixs[12];
	};

	struct LinearAnisortropicTetElement : public LinearAnisortropicElement{
		LinearAnisortropicTetElement(TetMesh *m) : LinearAnisortropicElement(m){}
		void setBMatrixs();
		void generateSubStiffnessMatrix(int aNodeIndex, int bNodeIndex, const Scalar *D, Scalar *result) const;
	private:
		Scalar BMatrixs[12];
	};

	struct InvertibleHyperelasticTetElement : public InvertibleHyperelasticElement {
		InvertibleHyperelasticTetElement(TetMesh *m): InvertibleHyperelasticElement(m){}
		void getPrecomputes(Scalar *drivates, Scalar *deforamtionGradients) const;
		int getDirvateEntryCount() const { return 12;}
		int getDeformGradientsPreEntryCount() const { return 9; }
		int getQuadraturePointCount() const { return 1; }

		void generateDeformationGradient(const Scalar *precompute, Scalar *gradients) const;
		void generateSubStiffnessMatrix(const Scalar *drivates, const Scalar *diags, const Scalar *leftOrthoMats, 
			const Scalar *rightOrthoMats, const Scalar *energyGradients, const Scalar *energyHassians, Scalar *result) const;
		void generateNodalVirtualWorks(const Scalar *precompute, const Scalar *stress, Scalar *result) const;
	private:
		void getdPdF(const Scalar *diag, const Scalar *energyGradient, const Scalar *energyHassian, Scalar dPdF[81]) const;
		void forceSemidefinite3x3(Scalar mat[6]) const;
		void forceSemidefinite2x2(Scalar& diag, Scalar& offdiag) const;
		int tensorIndex(int i, int j, int k, int l) const { return i * 27 + j * 9 + k * 3 + l; }
	};

	void getTetShapeFunctionDerivatives(const Vector3& a, const Vector3& b, const Vector3& c, const Vector3& d,
		Scalar *dn0, Scalar *dn1, Scalar *dn2, Scalar *dn3);

	inline Scalar getTetVolume(const Vector3& a, const Vector3& b, const Vector3& c, const Vector3& d) {
		Vector3 ab = b - a;
		Vector3 ac = c - a;
		Vector3 ad = d - a;

		return fabs(ab*(ac%ad)) / 6.f;
	}

	inline Scalar getTriArea(const Vector3& a, const Vector3& b, const Vector3& c) {
		return ((b - a) % (c - a)).length() * 0.5f;
	}
}

#endif