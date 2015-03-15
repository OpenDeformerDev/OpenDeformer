#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_CORE_PRECONDITIONER_H
#define ODER_CORE_PRECONDITIONER_H

#include "oder.h"

namespace ODER{
	class Preconditioner{
	public:
		virtual void solvePreconditionerSystem(const DenseVector& rhs, DenseVector& result) const = 0;
		virtual ~Preconditioner() = default;
	};
}

#endif