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

#define GLEW_STATIC
#include <GL/glew.h>
#include <GL/freeglut.h>

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

#define L1_CACHE_LINE_SIZE 64
#define DEFAULT_ARENA_OBJ_SIZE 1024
#define DEFAULT_POOL_OBJ_COUNT 128
#define GL_BUFFER_OFFSET( offset ) ((GLvoid*) (offset))

#ifndef NDEBUG
#define ODER_DEBUG
#endif

#ifdef ODER_DEBUG
#define Assert(expr) assert(expr)
#else
#define Assert(expr) ((void)0)
#endif

	static unsigned int randomSeed = 23;

	class SparseMatrixAssembler;
	class SparseMatrix;
	template<int, int> class BlockedSymSparseMatrix;
	template<int, int> class BlockedSymSparseMatrixAssembler;
	class Mesh;
	struct Element;
	struct Facet;
	template<class FT> struct VectorBase;
	struct Vector;
	class DenseVector;
	template<class FT> struct Tensor2;
	struct Quaternion;
	class MechMaterial;
	class Mesher;
	class MeshRelabeler;
	class Intergrator;
	class Forcer;
	class EigenSlover;
	class Constrainer;
	class NodeIndexer;
	class Simulator;

	template<class T> class MemoryArena;
	template<class T> class MemoryPool;

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