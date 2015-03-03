#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_CORE_PRECONDITIONER
#define ODER_CORE_PRECONDITIONER

#include "oder.h"

namespace ODER{
	class Preconditioner{
	public:
		virtual void solvePreconditionerSystem(const DenseVector& rhs, DenseVector& result);
		virtual ~Preconditioner() = default;
	};
}

#endif