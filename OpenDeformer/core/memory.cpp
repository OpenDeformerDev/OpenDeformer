#include "stdafx.h"
#include "memory.h"

namespace ODER{
	void freeAligned(void *p){
#if defined(ODER_IS_WINDOWS)
		_aligned_free(p);
#elif defined(ODER_IS_APPLE) || defined(ODER_IS_OPENBSD)
		ptrdiff_t adjustment = (reinterpret_cast<uint8_t *>(p))[-1];
		delete[] (reinterpret_cast<uint8_t *>((reinterpret_cast<uintptr_t>(p) - adjustment)));
#else
		free(p);
#endif
	}
}