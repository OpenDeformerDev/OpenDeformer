#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_CORE_MEMORY_H
#define ODER_CORE_MEMORY_H

#include "oder.h"
#include <cstdint>
#include <vector>
#include <atomic>

#if !defined(ODER_IS_APPLE) && !defined(ODER_IS_OPENBSD) && !defined(ODER_IS_WINDOWS)
#include <malloc.h>
#endif

namespace ODER{
	class ReferenceCounted{
	public:
		ReferenceCounted() : nReference(0){}
		std::atomic<int> nReference;
	};

	template <class T> class Reference{
	public:
		Reference(T *ptr = NULL) noexcept{
			p = ptr;
			if (p) p->nReference.fetch_add(1, std::memory_order_relaxed);
		}
		Reference(const Reference &r) noexcept{
			p = r.p;
			if (p) p->nReference.fetch_add(1, std::memory_order_relaxed);
		}
		Reference(Reference &&r) noexcept{
			p = r.p;
			r.p = NULL;
		}
		Reference &operator=(T *ptr) noexcept{
			if (ptr) ptr->nReference.fetch_add(1, std::memory_order_relaxed);
			if (p && p->nReference.fetch_sub(1, std::memory_order_release) == 1) {
				p->nReference.load(std::memory_order_acquire);
				delete p;
			}
			p = ptr;
			return *this;
		}
		Reference &operator=(const Reference &r) noexcept{
			if (r.p) r.p->nReference.fetch_add(1, std::memory_order_relaxed);
			if (p && p->nReference.fetch_sub(1, std::memory_order_release) == 1) {
				p->nReference.load(std::memory_order_acquire);
				delete p;
			}
			p = r.p;
			return *this;          
		}
		Reference &operator=(Reference &&r) noexcept{
			std::swap(p, r.p);
			return *this;
		}
		~Reference(){
			if (p && p->nReference.fetch_sub(1, std::memory_order_release) == 1) {
				p->nReference.load(std::memory_order_acquire);
				delete p;
			}
		}
		T *operator->() { return p; }
		const T *operator->() const{ return p; }
		operator bool() const{ return p != NULL; }

	private:
		T* p;
	};

	template<class T, size_t Align = std::alignment_of<T>::value> class MemoryArena{
	public:
		static_assert(std::alignment_of<T>::value <= Align, "ODER::MemoryArena requires larger Align parameter");
		static_assert((Align & (Align - 1)) == 0, "ODER::MemoryArena Align parameter should be power of 2");
		MemoryArena(unsigned int bs = ODER_DEFAULT_ARENA_OBJ_SIZE){
			blockSize = bs * sizeof(T);
			curBlockPos = 0;
			curBlock = allocAligned<int8_t>(blockSize);
		}
		MemoryArena(const MemoryArena& arena) = delete;
		MemoryArena(MemoryArena&& arena) = delete;
		MemoryArena& operator=(const MemoryArena& arena) = delete;
		MemoryArena& operator=(MemoryArena&& arena) = delete;
		template<class... Args> T* Alloc(unsigned int size = 1, Args&&... inits){
			constexpr unsigned int complement = Align - 1;
			unsigned int bytes = ((size*sizeof(T) + complement) & (~complement));
			if (curBlockPos + bytes > blockSize){
				usedBlocks.push_back(curBlock);
				if (avaBlocks.size() && size <= blockSize){
					curBlock = avaBlocks.back();
					avaBlocks.pop_back();
				}
				else
					curBlock = allocAligned<int8_t>(std::max(size, blockSize));
				curBlockPos = 0;
			}
			T *ret = (T *)(curBlock + curBlockPos);
			curBlockPos += bytes;

			Initiation(ret, size, std::forward<Args>(inits)...);
			return ret;
		}

		void freeAll(){
			curBlockPos = 0;
			unsigned int us = usedBlocks.size();
			for (unsigned int i = 0; i < us; i++){
				avaBlocks.push_back(usedBlocks.back());
				usedBlocks.pop_back();
			}
		}

		~MemoryArena(){
			freeAligned(curBlock);
			unsigned int us = usedBlocks.size(), as = avaBlocks.size();
			for (unsigned int i = 0; i < us; i++){
				freeAligned(usedBlocks[i]);
			}
			for (unsigned int i = 0; i < as; i++){
				freeAligned(avaBlocks[i]);
			}
		}

	private:
		unsigned int curBlockPos, blockSize;
		int8_t *curBlock;
		std::vector<int8_t *> usedBlocks, avaBlocks;
	};

	template<class T, size_t Align = 8> class MemoryPool{
	public:
		static_assert(std::alignment_of<T>::value <= Align, "ODER::MemoryPool requires larger Align parameter");
		static_assert(sizeof(void *) <= Align, "ODER::MemoryPool requires Align parameter larger than the size of a pointer");
		static_assert((Align & (Align - 1)) == 0, "ODER::MemoryPool Align parameter should be power of 2");

		MemoryPool(unsigned int count = ODER_DEFAULT_POOL_OBJ_COUNT) : curBlockPos(0){
			constexpr unsigned int complement = Align - 1;
			constexpr unsigned int objectByteCounts = ((sizeof(T) + complement) & (~complement));
			blockSize = count*objectByteCounts;
			curBlock = allocAligned<int8_t>(blockSize);
			deadStack = NULL;
		}
		MemoryPool(const MemoryPool&) = delete;
		MemoryPool(MemoryPool&&) = delete;
		MemoryPool& operator=(const MemoryPool&) = delete;
		MemoryPool& operator=(MemoryPool&&) = delete;
		template<class... Args> T* Alloc(Args&&... inits){
			constexpr unsigned int complement = Align - 1;
			constexpr unsigned int objectByteCounts = ((sizeof(T) + complement) & (~complement));

			T *ret = NULL;
			if (deadStack != NULL){
				ret = (T *)deadStack;
				deadStack = *(int8_t **)deadStack;
			}
			else{
				if (curBlockPos + objectByteCounts > blockSize){
					usedBlocks.push_back(curBlock);
					if (avaBlocks.size() > 0){
						curBlock = avaBlocks.back();
						avaBlocks.pop_back();
					}
					else
						curBlock = allocAligned<int8_t>(blockSize);
					curBlockPos = 0;
				}
				ret = (T *)(curBlock + curBlockPos);
				curBlockPos += objectByteCounts;
			}
			Construct(ret, std::forward<Args>(inits)...);
			return ret;
		}
		void Dealloc(T *item){
			Destory(item);
			*(int8_t **)item = deadStack;
			deadStack = (int8_t *)item;
		}
		void freeAll(){
			curBlockPos = 0;
			unsigned int us = usedBlocks.size();
			unsigned int deleteSize = (us + 1) >> 1;
			unsigned int remainSize = us - deleteSize;
			for (unsigned int i = 0; i < deleteSize; i++){
				freeAligned(usedBlocks.back());
				usedBlocks.pop_back();
			}
			for (unsigned int i = 0; i < remainSize; i++){
				avaBlocks.push_back(usedBlocks.back());
				usedBlocks.pop_back();
			}
			deadStack = NULL;
		}
		~MemoryPool(){
			freeAligned(curBlock);
			unsigned int us = usedBlocks.size(), as = avaBlocks.size();
			for (unsigned int i = 0; i < us; i++)
				freeAligned(usedBlocks[i]);
			for (unsigned int i = 0; i < as; i++)
				freeAligned(avaBlocks[i]);
		}
	private:
		unsigned int curBlockPos, blockSize;
		int8_t *curBlock;
		int8_t *deadStack;
		std::vector<int8_t *> usedBlocks, avaBlocks;
	};



	template<size_t Align = ODER_L1_CACHE_LINE_SIZE> void *allocAligned(size_t size){
		static_assert((Align & (Align - 1)) == 0, "ODER::allocAligned Align parameter should be power of 2");
		static_assert(Align >= 1 && Align <= 128, "ODER::allocAligned Align parameter should be set between 1 and 128");

#if defined(ODER_IS_WINDOWS)
		return _aligned_malloc(size, Align);
#elif defined(ODER_IS_APPLE) || defined(ODER_IS_OPENBSD)
		uintptr_t mem = reinterpret_cast<uintptr_t>(new uint8_t[size + Align]);
		ptrdiff_t adjustment = Align - (mem & (Align - 1));
		uintptr_t alignedMem = mem + adjustment;
		(reinterpret_cast<uint8_t *>(alignedMem))[-1] = static_cast<uint8_t>(adjustment);
		return reinterpret_cast<void *>(alignedMem);
#else
		return memalign(Align, size);
#endif
	}

	template <class T, size_t Align = ODER_L1_CACHE_LINE_SIZE> T *allocAligned(size_t num = 1){
		static_assert((Align & (Align - 1)) == 0, "ODER::allocAligned Align parameter should be power of 2");
		static_assert(Align >= 1 && Align <= 128, "ODER::allocAligned Align parameter should be set between 1 and 128");

		return (T *)allocAligned<Align>(num * sizeof(T));
	}

	void freeAligned(void *);

	template<class T, class... Args> inline void Initiation(T *vals, unsigned int size, Args&&... inits){
		for (unsigned int i = 0; i < size; i++){
			new (&vals[i]) T(std::forward<Args>(inits)...);
		}
	}
	template<> inline void Initiation(double *vals, unsigned int size){
		memset(vals, 0, sizeof(double) * size);
	}
	template<> inline void Initiation(float *vals, unsigned int size){
		memset(vals, 0, sizeof(float) * size);
	}
	template<> inline void Initiation(int *vals, unsigned int size){
		memset(vals, 0, sizeof(int) * size);
	}
	template<> inline void Initiation(unsigned int *vals, unsigned int size){
		memset(vals, 0, sizeof(unsigned int) * size);
	}

	template<class T> inline void Destory(T* p){
		p->~T();
	}
	template<class T> inline void Destory(T* p, unsigned int size){
		for (unsigned int i = 0; i < size; i++){
			p[i].~T();
		}
	}

	template<class T, class... Args> inline void Construct(T* p, Args&&... inits){
		new (p) T(std::forward<Args>(inits)...);
	}
}

#endif