#include "stdafx.h"
#include "element.h"

namespace ODER{
	LinearIsotropicElement::LinearIsotropicElement(Mesh *m, int entries): Element(m){
		BMatrixs = new double[entries];
	}

	LinearIsotropicElement::~LinearIsotropicElement(){
		delete[] BMatrixs;
	}

	ReducedIsotropicElement::ReducedIsotropicElement(Mesh *m, int entries) : Element(m){
		BMatrixs = new double[entries];
	}

	ReducedIsotropicElement::~ReducedIsotropicElement(){
		delete[] BMatrixs;
	}

	LinearAnisortropicElement::LinearAnisortropicElement(Mesh *m, int entries) : Element(m){
		BMatrixs = new double[entries];
	}

	LinearAnisortropicElement::~LinearAnisortropicElement(){
		delete[] BMatrixs;
	}

}