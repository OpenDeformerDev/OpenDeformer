#include "forcer.h"
#include "intergrator.h"

using namespace ODER;

Intergrator::Intergrator(int DOFS, double massDamp, double stiffDamp, double ts)
	:dofs(DOFS), massDamping(massDamp), stiffnessDamping(stiffDamp), timeStep(ts){
	d = allocAligned<double>(dofs);
	v = allocAligned<double>(dofs);
	a = allocAligned<double>(dofs);
	pre_d = allocAligned<double>(dofs);
	pre_v = allocAligned<double>(dofs);
	pre_a = allocAligned<double>(dofs);
	externalVirtualWork = allocAligned<double>(dofs);

	memset(d, 0, dofs*sizeof(double));
	memset(v, 0, dofs*sizeof(double));
	memset(a, 0, dofs*sizeof(double));
	memset(pre_d, 0, dofs*sizeof(double));
	memset(pre_v, 0, dofs*sizeof(double));
	memset(pre_a, 0, dofs*sizeof(double));
	memset(externalVirtualWork, 0, dofs*sizeof(double));
}

Intergrator::~Intergrator(){
	freeAligned(a);
	freeAligned(d);
	freeAligned(v);
	freeAligned(pre_a);
	freeAligned(pre_d);
	freeAligned(pre_v);
	freeAligned(externalVirtualWork);
}