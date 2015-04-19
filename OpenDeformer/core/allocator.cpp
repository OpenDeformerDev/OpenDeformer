#include "stdafx.h"
#include "allocator.h"

namespace ODER{
	Freelist::Freelist(unsigned int bs){
		blockSize = bs;
		curBlockPos = 0;
		perInitBytes = std::min(512U, blockSize);
		curBlock = (uint8_t *)malloc(blockSize);
		memset(freeLists, 0, sizeof(uint8_t *) * ODER_FREELIST_COUNT);
	}

	void *Freelist::Alloc(size_t size){
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

	void Freelist::initList(uint8_t *&list, size_t index){
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

	void Freelist::Dealloc(uint8_t *p, size_t size){
		if (size <= ODER_FREELIST_MAXBYTES){
			unsigned int bytes = roundUp(size);
			size_t index = getListIndex(bytes);
			*(uint8_t **)p = freeLists[index];
			freeLists[index] = (uint8_t *)p;
		}
		else
			free(p);
	}

	Freelist::~Freelist(){
		free(curBlock);
		unsigned int us = usedBlocks.size();
		for (unsigned int i = 0; i < us; i++)
			free(usedBlocks[i]);
	}
}