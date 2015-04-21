#include "stdafx.h"
#include "memory.h"
#if !defined(ODER_IS_APPLE) && !defined(ODER_IS_OPENBSD) && !defined(ODER_IS_WINDOWS)
#include <malloc.h>
#endif

namespace ODER{
	void *allocAligned(size_t size){
#if defined(ODER_IS_WINDOWS)
		return _aligned_malloc(size, ODER_L1_CACHE_LINE_SIZE);
#elif defined(ODER_IS_APPLE) || defined(ODER_IS_OPENBSD)
#if ODER_L1_CACHE_LINE_SIZE < 256
		using adj_type = uint8_t;
#else
		using adj_type = uint16_t;
#endif
		constexpr size_t offset = ODER_L1_CACHE_LINE_SIZE + sizeof(adj_type) - 1;
		uintptr_t mem = reinterpret_cast<uintptr_t>(new int8_t[size + offset]);
		ptrdiff_t adjustment = ODER_L1_CACHE_LINE_SIZE - (mem & (ODER_L1_CACHE_LINE_SIZE - 1));
		uintptr_t alignedMem = mem + adjustment;
		((adj_type *)alignedMem)[-1] = static_cast<adj_type>(adjustment);
		return (void *)alignedMem;
#else
		return memalign(ODER_L1_CACHE_LINE_SIZE, size);
#endif
	}

	void freeAligned(void *p){
#if defined(ODER_IS_WINDOWS)
		_aligned_free(p);
#elif defined(ODER_IS_APPLE) || defined(ODER_IS_OPENBSD)
#if ODER_L1_CACHE_LINE_SIZE < 256
		using adj_type = uint8_t;
#else
		using adj_type = uint16_t;
#endif
		ptrdiff_t adjustment = ((adj_type *)p)[-1];
		free((void *)(uintptr_t(p) - adjustment));
#else
		free(p);
#endif
	}
}