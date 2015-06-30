#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_CORE_ODER_H
#define ODER_CORE_ODER_H

#include <algorithm>
#include <functional>
#include <assert.h>

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
	static unsigned int randomSeed = 23;

	class SparseMatrixAssembler;
	class SparseMatrix;
	template<int blockLength, int blockWidth> class BlockedSymSparseMatrix;
	template<int blockLength, int blockWidth> class BlockedSymSparseMatrixAssembler;
	class Mesh;
	struct Element;
	struct GeometricElement;
	struct Facet;
	template<class FT> struct VectorBase;
	using Vector = VectorBase<float>;
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
	template<class SpMatrix> class LinearSolver;
	class Preconditioner;
	class Constrainer;
	class NodeIndexer;
	class Simulator;
	template<class FT> class ExactArthmeticer;

	template<class T, unsigned int Align> class MemoryArena;
	template<class T, unsigned int Align> class MemoryPool;

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