#include "stdafx.h"
#include "memory.h"

namespace ODER{
	void *allocAligned(size_t size){
		return _aligned_malloc(size, L1_CACHE_LINE_SIZE);
	}
	void *reallocAligned(void* memory, size_t size){
		return _aligned_realloc(memory, size, L1_CACHE_LINE_SIZE);
	}
	void freeAligned(void *p){
		_aligned_free(p);
	}
}