#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_CORE_ALLOCATOR_H
#define ODER_CORE_ALLOCATOR_H

#include <memory>
#include <vector>
#include <cstdint>
#include <limits>
#include "memory.h"

#define ODER_DEFAULT_FREELIST_BLOCK_SIZE 4096
#define ODER_FREELIST_ALIGN 8
#define ODER_FREELIST_MAXBYTES 128
#define ODER_FREELIST_COUNT (ODER_FREELIST_MAXBYTES / ODER_FREELIST_ALIGN)

namespace ODER{
	class ThreadSafeFreelist{// might be
	public:
		ThreadSafeFreelist(unsigned int bs = ODER_DEFAULT_FREELIST_BLOCK_SIZE);
		ThreadSafeFreelist(const ThreadSafeFreelist&) = delete;
		ThreadSafeFreelist(ThreadSafeFreelist&& list) = delete;
		ThreadSafeFreelist& operator=(const ThreadSafeFreelist&) = delete;
		ThreadSafeFreelist& operator=(ThreadSafeFreelist&& list) = delete;
		void *Alloc(size_t size);
		template<class T> T *Alloc(size_t size){
			T *ret = (T *)Alloc(sizeof(T) * size);
			return ret;
		}
		void Dealloc(void *p, size_t size);
		template<class T> void Dealloc(T *p, size_t size){
			Dealloc((void *)p, sizeof(T) * size);
		}
		~ThreadSafeFreelist();
	private:
		size_t roundUp(size_t bytes) const{ 
			constexpr unsigned int complement = ODER_FREELIST_ALIGN - 1;
			return (bytes + complement) & (~complement);
		}
		size_t getListIndex(size_t bytes) const{
			return (bytes + ODER_FREELIST_ALIGN - 1) / ODER_FREELIST_ALIGN - 1;
		}
		uintptr_t initList(size_t index);

		std::atomic<uintptr_t> freeLists[ODER_FREELIST_COUNT];

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

	class ThreadUnsafeFreelist{
	public:
		ThreadUnsafeFreelist(unsigned int bs = ODER_DEFAULT_FREELIST_BLOCK_SIZE);
		ThreadUnsafeFreelist(const ThreadUnsafeFreelist&) = delete;
		ThreadUnsafeFreelist(ThreadUnsafeFreelist&& list) = delete;
		ThreadUnsafeFreelist& operator=(const ThreadUnsafeFreelist&) = delete;
		ThreadUnsafeFreelist& operator=(ThreadUnsafeFreelist&& list) = delete;
		void *Alloc(size_t size);
		template<class T> T *Alloc(size_t size){
			T *ret = (T *)Alloc(sizeof(T) * size);
			return ret;
		}
		void Dealloc(uint8_t *p, size_t size);
		template<class T> void Dealloc(T *p, size_t size){
			Dealloc((uint8_t *)p, sizeof(T) * size);
		}
		~ThreadUnsafeFreelist();
	private:
		size_t roundUp(size_t bytes){
			constexpr unsigned int complement = ODER_FREELIST_ALIGN - 1;
			return (bytes + complement) & (~complement);
		}
		size_t getListIndex(size_t bytes){
			return (bytes + ODER_FREELIST_ALIGN - 1) / ODER_FREELIST_ALIGN - 1;
		}
		void initList(uint8_t *&list, size_t index);

		uint8_t *freeLists[ODER_FREELIST_COUNT];

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

	template<class T> using MultiThreadFreelistAllocator = FreelistAllocator<T, ThreadSafeFreelist>;
	template<class T> using SingleThreadFreelistAllocator = FreelistAllocator<T, ThreadUnsafeFreelist>;

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