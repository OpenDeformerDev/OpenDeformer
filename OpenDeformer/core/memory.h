#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_CORE_MEMORY_H
#define ODER_CORE_MEMORY_H

#include "oder.h"

namespace ODER{
	class ReferenceCounted{
	public:
		ReferenceCounted(){ nReference = 0; }
		atomic<int> nReference;
	};

	template <class T> class Reference{
	public:
		Reference(T *ptr = NULL){
			p = ptr;
			if (p) p->nReference++;
		}
		Reference(const Reference<T> &r){
			p = r.p;
			if (p) p->nReference++;
		}
		Reference(Reference<T> &&r){
			p = r.p;
			r.p = NULL;
		}
		Reference &operator=(T *ptr){
			if (ptr) ptr->nReference++;
			if (p && --p->nReference == 0) delete p;
			p = ptr;
			return *this;
		}
		Reference &operator=(const Reference<T> &r){
			if (r.p) r.p->nReference++;
			if (p && --p->nReference) delete p;
			p = r.p;
			return *this;
		}
		Reference &operator=(Reference<T> &&r){
			std::swap(p, r.p);
			return *this;
		}
		~Reference(){
			if (p && --p->nReference == 0) delete p;
		}
		T *operator->() { return p; }
		const T *operator->() const{ return p; }
		operator bool() const{ return p != NULL; }

	private:
		T* p;
	};

	template<class T> class MemoryArena{
	public:
		MemoryArena(unsigned int bs = DEFAULT_ARENA_OBJ_SIZE * sizeof(T)) :blockSize(bs), curBlockPos(0){
			curBlock = allocAligned<char>(blockSize);
		}
		T* Alloc(unsigned int size = 1){
			unsigned int bytes = ((size*sizeof(T) + 7) & (~7));
			if (curBlockPos + bytes > blockSize){
				usedBlocks.push_back(curBlock);
				if (avaBlocks.size() && size <= blockSize){
					curBlock = avaBlocks.back();
					avaBlocks.pop_back();
				}
				else
					curBlock = allocAligned<char>(max(size, blockSize));
				curBlockPos = 0;
			}
			T *ret = (T *)(curBlock + curBlockPos);
			curBlockPos += bytes;

			for (unsigned int i = 0; i < size; i++){
				new (&ret[i]) T();
			}
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
		char *curBlock;
		std::vector<char *> usedBlocks, avaBlocks;
	};

	template<class T> class MemoryPool{
	public:
		MemoryPool(unsigned int count = DEFAULT_POOL_OBJ_COUNT) : curBlockPos(0){
			objectByteCounts = ((sizeof(T) + 15) & (~15));
			blockSize = count*objectByteCounts;
			curBlock = allocAligned<char>(blockSize);
			deadStack = (char *)NULL;
		}
		T* Alloc(){
			T *ret = NULL;
			if (deadStack != (char *)NULL){
				ret = (T *)deadStack;
				deadStack = *(char **)deadStack;
			}
			else{
				if (curBlockPos + objectByteCounts > blockSize){
					usedBlocks.push_back(curBlock);
					if (avaBlocks.size() > 0){
						curBlock = avaBlocks.back();
						avaBlocks.pop_back();
					}
					else
						curBlock = allocAligned<char>(blockSize);
					curBlockPos = 0;
				}
				ret = (T *)(curBlock + curBlockPos);
				curBlockPos += objectByteCounts;
			}
			new (ret)T();
			return ret;
		}
		void Dealloc(T *item){
			*(char **)item = deadStack;
			deadStack = (char *)item;
		}
		void freeAll(){
			curBlockPos = 0;
			unsigned int us = usedBlocks.size();
			for (unsigned int i = 0; i < us; i++){
				avaBlocks.push_back(usedBlocks.back());
				usedBlocks.pop_back();
			}
			deadStack = (char *)NULL;
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
		unsigned int curBlockPos, objectByteCounts, blockSize;
		char *curBlock;
		char *deadStack;
		std::vector<char *> usedBlocks, avaBlocks;
	};

	void *allocAligned(size_t size);
	template <class T> T *allocAligned(size_t num = 1){
		return (T *)allocAligned(num*sizeof(T));
	}
	void *reallocAligned(void* memory, size_t size);
	template <class T> T *reallocAligned(T* memory, size_t num){
		return (T *)reallocAligned(memory, num*sizeof(T));
	}
	void freeAligned(void *);
}

#endif