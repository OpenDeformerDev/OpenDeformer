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
	};

	struct TetFacet : public Facet{
		TetFacet(TetMesh *m) :Facet(m){}
		void getSurfVirtualWorks(double surfForce[3], double *result) const;
		~TetFacet() = default;
	};


	struct LinearIsotropicTetElement : public LinearIsotropicElement{
		LinearIsotropicTetElement(TetMesh *m) : LinearIsotropicElement(m){}
		void setBMatrixs();
		void generateSubStiffnessMatrix(int aNodeIndex, int bNodeIndex, const double *D, double *result) const;
		~LinearIsotropicTetElement() = default;
	private:
		double BMatrixs[12];
	};

	struct ReducedIsotropicTetElement : public ReducedIsotropicElement{
		ReducedIsotropicTetElement(TetMesh *m) : ReducedIsotropicElement(m){}
		void setBMatrixs();
		void Intergration(const double *C, double *nlpart, double *nnpart) const;
		void generateSubStiffnessMatrix(int aNodeIndex, int bNodeIndex, const double *D, double *result) const;
		~ReducedIsotropicTetElement() = default;
	private:
		double BMatrixs[12];
	};

	struct LinearAnisortropicTetElement : public LinearAnisortropicElement{
		LinearAnisortropicTetElement(TetMesh *m) : LinearAnisortropicElement(m){}
		void setBMatrixs();
		void generateSubStiffnessMatrix(int aNodeIndex, int bNodeIndex, const double *D, double *result) const;
	private:
		double BMatrixs[12];
	};

	struct InvertibleHyperelasticTetElement : public InvertibleHyperelasticElement {
		InvertibleHyperelasticTetElement(TetMesh *m): InvertibleHyperelasticElement(m){}
		void getPrecomputes(double *drivates, double *deforamtionGradients) const;
		int getDirvateEntryCount() const { return 12;}
		int getDeformGradientsPreEntryCount() const { return 12; }
		int getQuadraturePointCount() const { return 1; }

		void generateDeformationGradient(const double *precompute, const double *u, double *gradients) const;
		void generateSubStiffnessMatrix(const double *drivates, const double *diags, const double *leftOrthoMats, 
			const double *rightOrthoMats, const double *energyGradients, const double *energyHassians, double *result) const;
		void generateNodalVirtualWorks(const double *precompute, const double *stress, double *result) const;
	private:
		void getdPdF(const double *diag, const double *leftOrthoMat,
			const double *rightOrthoMat, const double *energyGradient, const double *energyHassian, double dPdF[81]) const;
	};

	void getTetShapeFunctionDerivatives(const Vector& a, const Vector& b, const Vector& c, const Vector& d,
		double *dn0, double *dn1, double *dn2, double *dn3);

	inline float getTetVolume(const Vector& a, const Vector& b, const Vector& c, const Vector& d) {
		Vector ab = b - a;
		Vector ac = c - a;
		Vector ad = d - a;

		return fabsf(ab*(ac%ad)) / 6.f;
	}

	inline float getTriArea(const Vector& a, const Vector& b, const Vector& c) {
		return ((b - a) % (c - a)).length() * 0.5f;
	}
}

#endif