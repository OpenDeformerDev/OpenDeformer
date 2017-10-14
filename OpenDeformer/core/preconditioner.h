#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_CORE_PRECONDITIONER_H
#define ODER_CORE_PRECONDITIONER_H

#include "oder.h"
#include "memory.h"

namespace ODER{
	class Preconditioner : public ReferenceCounted {
	public:
		virtual void solvePreconditionerSystem(int width, const Scalar *rhs, Scalar *result) const {
			memcpy(result, rhs, width * sizeof(Scalar));
		}
		virtual void resetPreconditionerSystem(const BlockedSymSpMatrix& mat) {}
		virtual void resetPreconditionerSystem(const SparseMatrix& mat) {}
		virtual void Preprocess(const BlockedSymSpMatrix& mat) {}
		virtual void Preprocess(const SparseMatrix& mat) {}
		virtual ~Preconditioner() = default;
	};
}

#endif