#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_CORE_ODER_H
#define ODER_CORE_ODER_H

#include <math.h>
#include <algorithm>
#include <functional>
#include <assert.h>
#include <stddef.h>
#include <stdint.h>

#ifndef NDEBUG
#define ODER_DEBUG
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

#if defined(_MSC_VER)
#define _ENABLE_ATOMIC_ALIGNMENT_FIX
#endif

#ifdef ODER_DEBUG
#define Assert(expr) assert(expr)
#else
#define Assert(expr) ((void)0)
#endif

#define ODER_L1_CACHE_LINE_SIZE 64
#define ODER_DEFAULT_ARENA_OBJ_SIZE 1024
#define ODER_DEFAULT_POOL_OBJ_COUNT 128

#ifdef M_PI
#undef M_PI
#endif
#define M_PI 3.141592653589793238462643383279502884197169399375105820974944592308

namespace ODER{
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

	using BlockedSymSpMatrix = BlockedSymSparseMatrix< 3, 3 >;
	using BlockedSymSpMatrixAssembler = BlockedSymSparseMatrixAssembler< 3, 3 > ;

	namespace RandomnationInternal { static unsigned int randomSeed = 23u; }
	template<unsigned int choices> inline unsigned int Randomnation() {
		using RandomnationInternal::randomSeed;
		if (choices < 714025u) {
			randomSeed = (randomSeed * 1366u + 150889u) % 714025u;
			return randomSeed % choices;
		}
		else {
			unsigned int newRandom = (randomSeed * 1366u + 150889u) % 714025u;
			randomSeed = (newRandom * 1366u + 150889u) % 714025u;
			newRandom = newRandom * (choices / 714025u) + randomSeed;

			if (newRandom >= choices) return newRandom - choices;
			else return newRandom;
		}
	}

	template<class T> inline void hashCombine(size_t& seed, T val) {
		seed ^= std::hash<T>()(val) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
	}

	template<class T, class... Rest> inline void hashCombine(size_t& seed, T val, Rest... rest) {
		seed ^= std::hash<T>()(val) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		hashCombine(seed, rest...);
	}

	inline void Severe(const char *error){
		perror(error);
	}

	template<class FT> inline FT Clamp(FT val, FT low, FT high) {
		if (val < low) return low;
		else if (val > high) return high;
		else return val;
	}

	inline float nearestPowerOfTwo(float val) {
		int exp = 0;
		float significand = frexp(val, &exp);
		return fabs(significand) < 0.75f ? ldexp(1.f, exp - 1) : ldexp(1.f, exp);
	}

	inline double nearestPowerOfTwo(double val) {
		int exp = 0;
		double significand = frexp(val, &exp);
		return fabs(significand) < 0.75 ? ldexp(1.0, exp - 1) : ldexp(1.0, exp);
	}

	constexpr float const_pow(float base, int exp) noexcept {
		return exp < 0 ? 1.f / const_pow(base, -exp) :
			(exp == 0 ? 1.f :
				exp % 2 == 0 ? const_pow(base * base, exp / 2) :
				const_pow(base * base, (exp - 1) / 2) * base);
	}

	constexpr double const_pow(double base, int exp) noexcept {
		return exp < 0 ? 1.0 / const_pow(base, -exp) :
			(exp == 0 ? 1.0 :
				exp % 2 == 0 ? const_pow(base * base, exp / 2) :
				const_pow(base * base, (exp - 1) / 2) * base);
	}
}

#endif