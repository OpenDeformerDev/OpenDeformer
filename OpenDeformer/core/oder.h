#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_CORE_ODER_H
#define ODER_CORE_ODER_H

#include <vector>
#include <map>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <deque>
#include <functional>
#include <assert.h>
#include <atomic>
#include <set>
#include <cmath>

#define ODER_L1_CACHE_LINE_SIZE 64
#define ODER_DEFAULT_ARENA_OBJ_SIZE 1024
#define ODER_DEFAULT_POOL_OBJ_COUNT 128

#ifndef NDEBUG
#define ODER_DEBUG
#endif

#ifdef ODER_DEBUG
#define Assert(expr) assert(expr)
#else
#define Assert(expr) ((void)0)
#endif

#if defined(_WIN32) || defined(_WIN64)
#define ODER_IS_WINDOWS
#elif defined(__linux__)
#define ODER_IS_LINUX
#elif defined(__APPLE__)
#define ODER_IS_APPLE
#elif defined(__OpenBSD__)
#define ODER_IS_OPENBSD
#endif

namespace ODER{
	using std::vector;
	using std::map;
	using std::unordered_set;
	using std::unordered_map;
	using std::priority_queue;
	using std::queue;
	using std::deque;
	using std::set;
	using std::atomic;
	using std::min;
	using std::max;

	static unsigned int randomSeed = 23;

	class SparseMatrixAssembler;
	class SparseMatrix;
	template<int blockLength, int blockWidth> class BlockedSymSparseMatrix;
	template<int blockLength, int blockWidth> class BlockedSymSparseMatrixAssembler;
	class Mesh;
	struct Element;
	struct Facet;
	template<class FT> struct VectorBase;
	struct Vector;
	class DenseVector;
	class SparseVector;
	template<class FT> struct Tensor2;
	struct Quaternion;
	class MechMaterial;
	class Mesher;
	class MeshRelabeler;
	class Intergrator;
	class Forcer;
	class EigenSolver;
	class LinearSolver;
	class Preconditioner;
	class Constrainer;
	class NodeIndexer;
	class Simulator;

	template<class T> class MemoryArena;
	template<class T> class MemoryPool;

	using BlockedSymSpMatrix = BlockedSymSparseMatrix< 2, 4 >;
	using BlockedSymSpMatrixAssembler = BlockedSymSparseMatrixAssembler< 2, 4 > ;

	inline unsigned int Randomnation(unsigned int choices){
		randomSeed = (randomSeed * 1366l + 150889l) % 714025l;
		return randomSeed % choices;
	}

	template<class T> inline void hashCombine(size_t& seed, T val){
		seed ^= std::hash<T>()(val)+0x9e3779b9 + (seed << 6) + (seed >> 2);
	}

	inline void Severe(const char *error){
		perror(error);
	}

	template<class FT> inline FT Clamp(FT val, FT low, FT high) {
		if (val<low) return low;
		else if (val>high) return high;
		else return val;
	}
}

#endif