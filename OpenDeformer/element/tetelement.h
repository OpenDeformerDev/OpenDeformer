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

	void getTetShapeFunctionDerivatives(const Vector& a, const Vector& b, const Vector& c, const Vector& d,
		double *dndx, double *dndy, double *dndz);

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