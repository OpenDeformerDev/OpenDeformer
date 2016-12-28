#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_CORE_ALLOCATOR_H
#define ODER_CORE_ALLOCATOR_H

#include <memory>
#include <vector>
#include <limits>
#include "memory.h"

#define ODER_DEFAULT_FREELIST_BLOCK_SIZE 4096
#define ODER_FREELIST_ALIGN 8
#define ODER_FREELIST_MAXBYTES 128

namespace ODER{
	template<size_t Align, size_t MaxBytes>
	class ThreadSafeFreelist{// might be
		static_assert((Align & (Align - 1)) == 0, "ODER::ThreadSafeFreelist Align parameter should be power of 2");
		static_assert(MaxBytes % Align == 0, "ODER::ThreadSafeFreelist MaxBytes parameter should be power of Align");
	public:
		ThreadSafeFreelist(unsigned int bs = ODER_DEFAULT_FREELIST_BLOCK_SIZE);
		ThreadSafeFreelist(const ThreadSafeFreelist&) = delete;
		ThreadSafeFreelist(ThreadSafeFreelist&& list) = delete;
		ThreadSafeFreelist& operator=(const ThreadSafeFreelist&) = delete;
		ThreadSafeFreelist& operator=(ThreadSafeFreelist&& list) = delete;
		template<class T> T *Alloc(size_t size){
			T *ret = (T *)Alloc(sizeof(T) * size);
			return ret;
		}
		template<class T> void Dealloc(T *p, size_t size){
			Dealloc((void *)p, sizeof(T) * size);
		}
		void freeAll();
		~ThreadSafeFreelist();
	private:
		void *Alloc(size_t size);
		void Dealloc(void *p, size_t size);
		size_t roundUp(size_t bytes) const{ 
			constexpr size_t complement = Align - 1;
			return (bytes + complement) & (~complement);
		}
		size_t getListIndex(size_t bytes) const{
			return (bytes + Align - 1) / Align - 1;
		}
		uintptr_t initList(size_t index);

		std::atomic<uintptr_t> freeLists[MaxBytes / Align];

		struct MemBlock{
			MemBlock() noexcept :usedBlockBytes(std::numeric_limits<unsigned int>::max()), blockHead(0){}
			unsigned int usedBlockBytes;
			uintptr_t blockHead;
		};

		unsigned int blockSize;
		unsigned int perInitBytes;
		std::atomic<MemBlock> curBlock;
		std::atomic<uintptr_t> usedBlocks;
	};

	template<size_t Align, size_t MaxBytes>
	class ThreadUnsafeFreelist{
		static_assert((Align & (Align - 1)) == 0, "ODER::ThreadUnsafeFreelist Align parameter should be power of 2");
		static_assert(MaxBytes % Align == 0, "ODER::ThreadUnsafeFreelist MaxBytes parameter should be power of Align");
	public:
		ThreadUnsafeFreelist(unsigned int bs = ODER_DEFAULT_FREELIST_BLOCK_SIZE);
		ThreadUnsafeFreelist(const ThreadUnsafeFreelist&) = delete;
		ThreadUnsafeFreelist(ThreadUnsafeFreelist&& list) = delete;
		ThreadUnsafeFreelist& operator=(const ThreadUnsafeFreelist&) = delete;
		ThreadUnsafeFreelist& operator=(ThreadUnsafeFreelist&& list) = delete;
		template<class T> T *Alloc(size_t size){
			T *ret = (T *)Alloc(sizeof(T) * size);
			return ret;
		}
		template<class T> void Dealloc(T *p, size_t size){
			Dealloc((uint8_t *)p, sizeof(T) * size);
		}
		void freeAll();
		~ThreadUnsafeFreelist();
	private:
		void *Alloc(size_t size);
		void Dealloc(uint8_t *p, size_t size);
		size_t roundUp(size_t bytes){
			constexpr size_t complement = Align - 1;
			return (bytes + complement) & (~complement);
		}
		size_t getListIndex(size_t bytes){
			return (bytes + Align - 1) / Align - 1;
		}
		void initList(uint8_t *&list, size_t index);

		uint8_t *freeLists[MaxBytes / Align];

		unsigned int blockSize;
		unsigned int curBlockPos;
		unsigned int perInitBytes;
		uint8_t *curBlock;
		std::vector<uint8_t *> usedBlocks;
	};

	template<class T, class Freelist> class FreelistAllocator{
	public:
		using size_type = size_t;
		using difference_type = ptrdiff_t;
		using value_type = T;
		using pointer = T*;
		using const_pointer = const T*;
		using reference = T&;
		using const_reference = const T&;

		using propagate_on_container_copy_assignment = std::true_type;
		using propagate_on_container_move_assignment = std::true_type;
		using propagate_on_container_swap = std::true_type;

		template<class U> struct rebind { using other = FreelistAllocator<U, Freelist>; };

		FreelistAllocator() noexcept : alloc(std::make_shared<Freelist>()){}
		FreelistAllocator(const FreelistAllocator&) = default;
		FreelistAllocator(FreelistAllocator&& allocator) = default;
		FreelistAllocator& operator=(const FreelistAllocator&) = default;
		FreelistAllocator& operator=(FreelistAllocator&& allocator) = default;

		template<class U> FreelistAllocator(const FreelistAllocator<U, Freelist>& allocator) noexcept : alloc(allocator.alloc){}
		template<class U> FreelistAllocator(FreelistAllocator<U, Freelist>&& allocator) noexcept : alloc(std::move(allocator.alloc)){}
		pointer address(reference x) const{ return &x;}
		const_pointer address(const_reference x) const{ return &x; }

		pointer allocate(size_type n, const void *hint = 0){
			pointer ret = alloc->Alloc<value_type>(n);
			return ret;
		}

		void deallocate(pointer p, size_type n){
			alloc->Dealloc<value_type>(p, n);
		}

		size_type max_size() const{ return std::numeric_limits<size_type>::max(); }

		template<class U, class... Args> void construct(U* p, Args&&... args){
			Construct(p, std::forward<Args>(args)...);
		}
		template<class U> void destroy(U* p){
			Destory(p);
		}
		~FreelistAllocator() = default;

		std::shared_ptr<Freelist> alloc;
	};

	template<size_t Align, size_t MaxBytes> ThreadSafeFreelist<Align, MaxBytes>::ThreadSafeFreelist(unsigned int bs) {
		blockSize = bs + std::alignment_of<std::max_align_t>::value;
		perInitBytes = std::min(512U, bs);
		MemBlock memBlock;
		memBlock.usedBlockBytes = 0;

		uintptr_t block = (uintptr_t)malloc(blockSize);
		*(uintptr_t *)block = 0;
		usedBlocks.store(block, std::memory_order_relaxed);
		memBlock.blockHead = block + std::alignment_of<std::max_align_t>::value;
		curBlock.store(memBlock, std::memory_order_relaxed);
		for (int i = 0; i < MaxBytes / Align; i++)
			freeLists[i].store(0, std::memory_order_relaxed);
	}

	template<size_t Align, size_t MaxBytes> void ThreadSafeFreelist<Align, MaxBytes>::freeAll() {
		uintptr_t block = *(uintptr_t *)(usedBlocks.load(std::memory_order_relaxed));
		while (block != 0) {
			uintptr_t freed = block;
			block = *(uintptr_t *)block;
			free((void *)freed);
		}

		MemBlock memBlock;
		memBlock.usedBlockBytes = 0;
		memBlock.blockHead = block + std::alignment_of<std::max_align_t>::value;
		usedBlocks.store(block, std::memory_order_relaxed);
		curBlock.store(memBlock, std::memory_order_relaxed);
		for (int i = 0; i < MaxBytes / Align; i++)
			freeLists[i].store(0, std::memory_order_relaxed);
	}

	template<size_t Align, size_t MaxBytes> void *ThreadSafeFreelist<Align, MaxBytes>::Alloc(size_t size) {
		if (size <= MaxBytes) {
			unsigned int bytes = roundUp(size);
			size_t index = getListIndex(bytes);

			uintptr_t ret = 0;
			uintptr_t head = freeLists[index].load(std::memory_order_acquire);
			do {
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

	template<size_t Align, size_t MaxBytes> uintptr_t ThreadSafeFreelist<Align, MaxBytes>::initList(size_t index) {
		unsigned int size = (index + 1) * Align;
		unsigned int objectCount = perInitBytes / size;

		MemBlock preBlock = curBlock.load(std::memory_order_relaxed), newBlock;
		//thread-safe(which depends on malloc) but might waste some memory space
		do {
		nextTry:
			newBlock.usedBlockBytes = preBlock.usedBlockBytes + perInitBytes;
			newBlock.blockHead = preBlock.blockHead + perInitBytes;
			if (newBlock.usedBlockBytes > blockSize) {
				uintptr_t raw = (uintptr_t)malloc(blockSize);
				newBlock.usedBlockBytes = perInitBytes + std::alignment_of<std::max_align_t>::value;
				newBlock.blockHead = raw + perInitBytes + std::alignment_of<std::max_align_t>::value;

				uintptr_t usedBlockHead = usedBlocks.load(std::memory_order_relaxed);
				*(uintptr_t *)raw = usedBlockHead;
				if (!usedBlocks.compare_exchange_strong(usedBlockHead, raw, std::memory_order_relaxed)) {
					free((void *)raw);
					preBlock = curBlock.load(std::memory_order_relaxed);
					goto nextTry;
				}
			}
		} while (!curBlock.compare_exchange_strong(preBlock, newBlock, std::memory_order_relaxed));


		uintptr_t ret = newBlock.blockHead - perInitBytes;
		uintptr_t iter = ret + size;
		for (unsigned int i = 0; i < objectCount - 2; i++) {
			*(uintptr_t *)iter = iter + size;
			iter += size;
		}

		*(uintptr_t *)iter = freeLists[index].load(std::memory_order_relaxed);
		while (!freeLists[index].compare_exchange_weak(*(uintptr_t *)iter, ret + size,
			std::memory_order_release, std::memory_order_relaxed));

		return ret;
	}

	template<size_t Align, size_t MaxBytes> void ThreadSafeFreelist<Align, MaxBytes>::Dealloc(void *p, size_t size) {
		if (size <= MaxBytes) {
			unsigned int bytes = roundUp(size);
			size_t index = getListIndex(bytes);
			*(uintptr_t *)p = freeLists[index].load(std::memory_order_relaxed);
			while (!freeLists[index].compare_exchange_weak(*(uintptr_t *)p, (uintptr_t)p,
				std::memory_order_release, std::memory_order_relaxed));
		}
		else
			free(p);
	}

	template<size_t Align, size_t MaxBytes> ThreadSafeFreelist<Align, MaxBytes>::~ThreadSafeFreelist() {
		uintptr_t p = usedBlocks.load(std::memory_order_relaxed);
		while (p != 0) {
			uintptr_t freed = p;
			p = *(uintptr_t *)p;
			free((void *)freed);
		}
	}

	template<size_t Align, size_t MaxBytes> ThreadUnsafeFreelist<Align, MaxBytes>::ThreadUnsafeFreelist(unsigned int bs) {
		blockSize = bs;
		curBlockPos = 0;
		perInitBytes = std::min(512U, blockSize);
		curBlock = (uint8_t *)malloc(blockSize);
		memset(freeLists, 0, sizeof(uint8_t *) * (MaxBytes / Align));
	}

	template<size_t Align, size_t MaxBytes> void ThreadUnsafeFreelist<Align, MaxBytes>::freeAll() {
		curBlockPos = 0;
		size_t us = usedBlocks.size();
		memset(freeLists, 0, sizeof(uint8_t *) * (MaxBytes / Align));
		for (auto block : usedBlocks) free(block);
	}

	template<size_t Align, size_t MaxBytes> void *ThreadUnsafeFreelist<Align, MaxBytes>::Alloc(size_t size) {
		if (size <= MaxBytes) {
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

	template<size_t Align, size_t MaxBytes> void ThreadUnsafeFreelist<Align, MaxBytes>::initList(uint8_t *&list, size_t index) {
		unsigned int size = (index + 1) * Align;
		unsigned int objectCount = perInitBytes / size;

		if (curBlockPos + perInitBytes > blockSize) {
			usedBlocks.push_back(curBlock);
			curBlock = (uint8_t *)malloc(blockSize);
			curBlockPos = 0;
		}
		uint8_t *ret = curBlock + curBlockPos;
		curBlockPos += perInitBytes;

		uint8_t *iter = ret;
		for (unsigned int i = 0; i < objectCount - 1; i++) {
			*(uint8_t **)iter = iter + size;
			iter += size;
		}
		*(uint8_t **)iter = NULL;
		list = ret;
	}

	template<size_t Align, size_t MaxBytes> void ThreadUnsafeFreelist<Align, MaxBytes>::Dealloc(uint8_t *p, size_t size) {
		if (size <= MaxBytes) {
			unsigned int bytes = roundUp(size);
			size_t index = getListIndex(bytes);
			*(uint8_t **)p = freeLists[index];
			freeLists[index] = (uint8_t *)p;
		}
		else
			free(p);
	}

	template<size_t Align, size_t MaxBytes> ThreadUnsafeFreelist<Align, MaxBytes>::~ThreadUnsafeFreelist() {
		free(curBlock);
		for (auto block : usedBlocks) free(block);
	}

	template<class T> using MultiThreadFreelistAllocator = FreelistAllocator<T, ThreadSafeFreelist<ODER_FREELIST_ALIGN, ODER_FREELIST_MAXBYTES>>;
	template<class T> using SingleThreadFreelistAllocator = FreelistAllocator<T, ThreadUnsafeFreelist<ODER_FREELIST_ALIGN, ODER_FREELIST_MAXBYTES>>;

	template<class T1, class T2, class Freelist1, class Freelist2>
	bool operator==(const FreelistAllocator<T1, Freelist1>& lhs, const FreelistAllocator<T2, Freelist2>& rhs){
		return std::is_same<Freelist1, Freelist2>::value && lhs.alloc == rhs.alloc;
	}
	template<class T1, class T2, class Freelist1, class Freelist2>
	bool operator!=(const FreelistAllocator<T1, Freelist1>& lhs, const FreelistAllocator<T2, Freelist2>& rhs){
		return !std::is_same<Freelist1, Freelist2>::value || lhs.alloc != rhs.alloc;
	}
}

#endif