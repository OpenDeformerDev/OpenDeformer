#include "stdafx.h"
#include "allocator.h"

namespace ODER{
	ThreadSafeFreelist::ThreadSafeFreelist(unsigned int bs){
		blockSize = bs + std::alignment_of<std::max_align_t>::value;
		perInitBytes = std::min(512U, bs);
		MemBlock memBlock;
		memBlock.usedBlockBytes = 0;

		uintptr_t block = (uintptr_t)malloc(blockSize);
		*(uintptr_t *)block = 0;
		usedBlocks = block;
		memBlock.blockHead = block + std::alignment_of<std::max_align_t>::value;
		curBlock.store(memBlock, std::memory_order_relaxed);
		for (int i = 0; i < ODER_FREELIST_COUNT; i++)
			freeLists[i].store(0, std::memory_order_relaxed);
	}

	void *ThreadSafeFreelist::Alloc(size_t size){
		if (size <= ODER_FREELIST_MAXBYTES){
			unsigned int bytes = roundUp(size);
			size_t index = getListIndex(bytes);

			uintptr_t ret = 0;
			uintptr_t head = freeLists[index].load(std::memory_order_acquire);
			do{
				if (head != 0)
					ret = head;
				else
					return (void *)initList(index);
			} while (!freeLists[index].compare_exchange_weak(head, *(uintptr_t *)ret, 
				std::memory_order_acquire, std::memory_order_relaxed));

			return (void *)ret;
		}
		else
			return malloc(size);
	}

	uintptr_t ThreadSafeFreelist::initList(size_t index){
		unsigned int size = (index + 1) * ODER_FREELIST_ALIGN;
		unsigned int objectCount = perInitBytes / size;

		MemBlock preBlock = curBlock.load(std::memory_order_relaxed), newBlock;
		//thread-safe(which depends on malloc) but might waste some memory space
		do{
		nextTry:
			newBlock.usedBlockBytes = preBlock.usedBlockBytes + perInitBytes;
			newBlock.blockHead = preBlock.blockHead + perInitBytes;
			if (newBlock.usedBlockBytes > blockSize){
				uintptr_t raw = (uintptr_t)malloc(blockSize);
				newBlock.usedBlockBytes = perInitBytes + std::alignment_of<std::max_align_t>::value;
				newBlock.blockHead = raw + perInitBytes + std::alignment_of<std::max_align_t>::value;

				uintptr_t usedBlockHead = usedBlocks.load(std::memory_order_relaxed);
				*(uintptr_t *)raw = usedBlockHead;
				if (!usedBlocks.compare_exchange_strong(usedBlockHead, raw, std::memory_order_relaxed)){
					free((void *)raw);
					preBlock = curBlock.load(std::memory_order_relaxed);
					goto nextTry;
				}
			}
		} while (!curBlock.compare_exchange_strong(preBlock, newBlock, std::memory_order_relaxed));


		uintptr_t ret = newBlock.blockHead - perInitBytes;
		uintptr_t iter = ret + size;
		for (unsigned int i = 0; i < objectCount - 2; i++){
			*(uintptr_t *)iter = iter + size;
			iter += size;
		}

		*(uintptr_t *)iter = freeLists[index].load(std::memory_order_relaxed);
		while (!freeLists[index].compare_exchange_weak(*(uintptr_t *)iter, ret + size,
			std::memory_order_release, std::memory_order_relaxed));

		return ret;
	}

	void ThreadSafeFreelist::Dealloc(void *p, size_t size){
		if (size <= ODER_FREELIST_MAXBYTES){
			unsigned int bytes = roundUp(size);
			size_t index = getListIndex(bytes);
			*(uintptr_t *)p = freeLists[index].load(std::memory_order_relaxed);
			while (!freeLists[index].compare_exchange_weak(*(uintptr_t *)p, (uintptr_t)p, 
				std::memory_order_release, std::memory_order_relaxed));
		}
		else
			free(p);
	}

	ThreadSafeFreelist::~ThreadSafeFreelist(){
		uintptr_t p = usedBlocks.load(std::memory_order_relaxed);
		while (p != 0){
			uintptr_t freed = p;
			p = *(uintptr_t *)p;
			free((void *)freed);
		}
	}

	ThreadUnsafeFreelist::ThreadUnsafeFreelist(unsigned int bs){
		blockSize = bs;
		curBlockPos = 0;
		perInitBytes = std::min(512U, blockSize);
		curBlock = (uint8_t *)malloc(blockSize);
		memset(freeLists, 0, sizeof(uint8_t *) * ODER_FREELIST_COUNT);
	}

	void *ThreadUnsafeFreelist::Alloc(size_t size){
		if (size <= ODER_FREELIST_MAXBYTES){
			unsigned int bytes = roundUp(size);
			size_t index = getListIndex(bytes);
			uint8_t *list = freeLists[index];
			if (!list)
				initList(list, index);
			freeLists[index] = *((uint8_t **)list);
			return list;
		}
		else
			return malloc(size);
	}

	void ThreadUnsafeFreelist::initList(uint8_t *&list, size_t index){
		unsigned int size = (index + 1) * ODER_FREELIST_ALIGN;
		unsigned int objectCount = perInitBytes / size;

		if (curBlockPos + perInitBytes > blockSize){
			usedBlocks.push_back(curBlock);
			curBlock = (uint8_t *)malloc(blockSize);
			curBlockPos = 0;
		}
		uint8_t *ret = curBlock + curBlockPos;
		curBlockPos += perInitBytes;

		uint8_t *iter = ret;
		for (unsigned int i = 0; i < objectCount - 1; i++){
			*(uint8_t **)iter = iter + size;
			iter += size;
		}
		*(uint8_t **)iter = NULL;
		list = ret;
	}

	void ThreadUnsafeFreelist::Dealloc(uint8_t *p, size_t size){
		if (size <= ODER_FREELIST_MAXBYTES){
			unsigned int bytes = roundUp(size);
			size_t index = getListIndex(bytes);
			*(uint8_t **)p = freeLists[index];
			freeLists[index] = (uint8_t *)p;
		}
		else
			free(p);
	}

	ThreadUnsafeFreelist::~ThreadUnsafeFreelist(){
		free(curBlock);
		unsigned int us = usedBlocks.size();
		for (unsigned int i = 0; i < us; i++)
			free(usedBlocks[i]);
	}
}