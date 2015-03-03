#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_CORE_LINEARSLOVER_H
#define ODER_CORE_LINEARSLOVER_H

#include "oder.h"

namespace ODER{
    class LinearSolver{
	public:
		virtual void solveLinearSystem(const DenseVector& rhs, DenseVector& result) const = 0;
		virtual ~LinearSolver() = default;
	};
}

#endif