#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_CORE_LATOOL_H
#define ODER_CORE_LATOOL_H

#include "oder.h"
#include "datastructure.h"

namespace ODER{
	template<class FT> struct VectorBase{
		VectorBase(){ x = y = z = FT(0); }
		VectorBase(FT xx, FT yy, FT zz) :x(xx), y(yy), z(zz){}

		VectorBase operator+(const VectorBase &v) const{
			return VectorBase(x + v.x, y + v.y, z + v.z);
		}

		VectorBase operator-(const VectorBase &v) const{
			return VectorBase(x - v.x, y - v.y, z - v.z);
		}

		VectorBase &operator+=(const VectorBase &v) {
			x += v.x; y += v.y; z += v.z;
			return *this;
		}

		VectorBase &operator-=(const VectorBase &v){
			x -= v.x; y -= v.y; z -= v.z;
			return *this;
		}

		VectorBase operator*(FT f) const{
			return VectorBase(f*x, f*y, f*z);
		}

		friend inline VectorBase operator*(FT f, const VectorBase &v){
			return v*f;
		}

		FT operator*(const VectorBase &v) const{
			return x*v.x + y*v.y + z*v.z;
		}

		VectorBase &operator*=(FT f){
			x *= f; y *= f; z *= f;
			return *this;
		}

		VectorBase operator/(FT f) const{
			Assert(f != 0);
			FT inf = FT(1.0) / f;
			return VectorBase(inf*x, inf*y, inf*z);
		}


		VectorBase &operator/=(FT f){
			Assert(f != 0);
			FT inf = FT(1.0) / f;
			x *= inf; y *= inf; z *= inf;
			return *this;
		}

		VectorBase operator-() const{
			return VectorBase(-x, -y, -z);
		}

		VectorBase operator%(const VectorBase& v) const{ //this is cross
			return VectorBase(y*v.z - z*v.y, z*v.x - x*v.z, x*v.y - y*v.x);
		}


		Tensor2<FT> operator^(const VectorBase& v) const;

		FT operator[](int i) const{
			Assert(i > -1 && i<3);
			return (&x)[i];
		}
		FT &operator[](int i){
			Assert(i>-1 && i<3);
			return (&x)[i];
		}

		FT length() const{
			return sqrt(x*x + y*y + z*z);
		}
		FT length2() const{
			return x*x + y*y + z*z;
		}

		bool operator==(const VectorBase& v) const{
			return x == v.x && y == v.y && z == v.z;
		}
		bool operator!=(const VectorBase& v) const{
			return x != v.x || y != v.y || z != v.z;
		}

		bool operator>(const VectorBase &v) const{
			if (x == v.x){
				if (y == v.y)
					return z > v.z;
				return y > v.y;
			}
			return x > v.x;
		}

		bool hasNaNs() const{
			return isnan(x) || isnan(y) || isnan(z);
		}

		bool hasInfs() const {
			return isinf(x) || isinf(y) || isinf(z);
		}

		FT x, y, z;
	};

	template<class FT> struct Tensor2{
		Tensor2(){
			m[0][0] = m[0][1] = m[0][2] =
			m[1][0] = m[1][1] = m[1][2] = 
			m[2][0] = m[2][1] = m[2][2] = FT(0);
		}
		Tensor2(const FT *mat) {
			memcpy(m, mat, 9 * sizeof(FT));
		}
		Tensor2(FT m00, FT m01, FT m02,
			FT m10, FT m11, FT m12,
			FT m20, FT m21, FT m22){
			m[0][0] = m00; m[0][1] = m01; m[0][2] = m02;
			m[1][0] = m10; m[1][1] = m11; m[1][2] = m12;
			m[2][0] = m20; m[2][1] = m21; m[2][2] = m22;
		}
		bool operator==(const Tensor2& mat) const{
			for (int i = 0; i < 3; i++){
				for (int j = 0; j < 3; j++){
					if (m[i][j] != mat.m[i][j])
						return false;
				}
			}
			return true;
		}
		bool operator!=(const Tensor2& mat) const{
			for (int i = 0; i < 3; i++){
				for (int j = 0; j < 3; j++){
					if (m[i][j] != mat.m[i][j])
						return true;
				}
			}
			return false;
		}
		Tensor2 operator*(const Tensor2& mat) const{
			Tensor2 r;
			for (int i = 0; i < 3; i++){
				for (int j = 0; j < 3; j++){
					r.m[i][j] = m[i][0] * mat.m[0][j] + m[i][1] * mat.m[1][j] +
						m[i][2] * mat.m[2][j];
				}
			}
			return r;
		}
		VectorBase<FT> operator*(const VectorBase<FT>& v) const{
			FT x = v.x, y = v.y, z = v.z;
			return VectorBase<FT>(m[0][0] * x + m[0][1] * y + m[0][2] * z, m[1][0] * x + m[1][1] * y + m[1][2] * z, m[2][0] * x + m[2][1] * y + m[2][2] * z);
		}
		FT operator&(const Tensor2& t) const{
			FT ret = 0.0;
			for (int i = 0; i < 3; i++){
				for (int j = 0; j < 3; j++){
					ret += m[i][j] * t.m[i][j];
				}
			}
			return ret;
		}
		FT operator()(int row, int column) const{
			Assert(row > -1 && row < 3 && column > -1 && column < 3);
			return m[row][column];
		}
		FT& operator()(int row, int column){
			Assert(row > -1 && row < 3 && column > -1 && column < 3);
			return m[row][column];
		}
		FT Determinant() const{
			return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
				- m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
				+ m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
		}
		friend Tensor2 Inverse(const Tensor2& tensor){
			Tensor2 r;

			FT det = tensor.Determinant();

			Assert(fabs(det) > 1e-7);

			FT invDet = 1.f / det;
			r.m[0][0] = (tensor.m[1][1] * tensor.m[2][2] - tensor.m[1][2] * tensor.m[2][1])*invDet;
			r.m[0][1] = (tensor.m[0][2] * tensor.m[2][1] - tensor.m[0][1] * tensor.m[2][2])*invDet;
			r.m[0][2] = (tensor.m[0][1] * tensor.m[1][2] - tensor.m[0][2] * tensor.m[1][1])*invDet;
			r.m[1][0] = (tensor.m[1][2] * tensor.m[2][0] - tensor.m[1][0] * tensor.m[2][2])*invDet;
			r.m[1][1] = (tensor.m[0][0] * tensor.m[2][2] - tensor.m[0][2] * tensor.m[2][0])*invDet;
			r.m[1][2] = (tensor.m[0][2] * tensor.m[1][0] - tensor.m[0][0] * tensor.m[1][2])*invDet;
			r.m[2][0] = (tensor.m[1][0] * tensor.m[2][1] - tensor.m[1][1] * tensor.m[2][0])*invDet;
			r.m[2][1] = (tensor.m[0][1] * tensor.m[2][0] - tensor.m[0][0] * tensor.m[2][1])*invDet;
			r.m[2][2] = (tensor.m[0][0] * tensor.m[1][1] - tensor.m[0][1] * tensor.m[1][0])*invDet;

			return r;
		}
		friend Tensor2 Transpose(const Tensor2& m){
			Tensor2 ret;
			for (int i = 0; i < 3; i++){
				for (int j = 0; j < 3; j++){
					ret.m[i][j] = m.m[j][i];
				}
			}
			return ret;
		}
	private:
		FT m[3][3];
	};

	struct Mat4{
		Mat4(){
			m[0][0] = m[1][1] = m[2][2] = m[3][3] = 1.0f;
			m[0][1] = m[0][2] = m[0][3] =
				m[1][0] = m[1][2] = m[1][3] =
				m[2][0] = m[2][1] = m[2][3] =
				m[3][0] = m[3][1] = m[3][2] = 0.0f;
		}
		Mat4(float mat[4][4]);
		Mat4(float m00, float m01, float m02, float m03,
			float m10, float m11, float m12, float m13,
			float m20, float m21, float m22, float m23,
			float m30, float m31, float m32, float m33);
		Mat4(const Quaternion &q);
		bool operator==(const Mat4& mat) const{
			for (int i = 0; i < 4; i++){
				for (int j = 0; j < 4; j++){
					if (m[i][j] != mat.m[i][j])
						return false;
				}
			}
			return true;
		}
		bool operator!=(const Mat4& mat) const{
			for (int i = 0; i < 4; i++){
				for (int j = 0; j < 4; j++){
					if (m[i][j] != mat.m[i][j])
						return true;
				}
			}
			return false;
		}
		Mat4 operator*(const Mat4& mat) const{
			Mat4 r;
			for (int i = 0; i < 4; i++){
				for (int j = 0; j < 4; j++){
					r.m[i][j] = m[i][0] * mat.m[0][j] + m[i][1] * mat.m[1][j] +
						m[i][2] * mat.m[2][j] + m[i][3] * mat.m[3][j];
				}
			}
			return r;
		}
		Mat4 operator*(float f) const{
			Mat4 r;
			for (int i = 0; i < 4; i++){
				for (int j = 0; j < 4; j++){
					r.m[i][j] = m[i][j] * f;
				}
			}
			return r;
		}
		friend inline Mat4 operator*(float f, const Mat4 &mat){
			return mat*f;
		}
		Vector3f operator*(const Vector3f& v) const{
			float x = v.x, y = v.y, z = v.z;
			return Vector3f(m[0][0] * x + m[0][1] * y + m[0][2] * z, m[1][0] * x + m[1][1] * y + m[1][2] * z, m[2][0] * x + m[2][1] * y + m[2][2] * z);
		}
		friend Mat4 Transpose(const Mat4&);
		friend Mat4 Inverse(const Mat4&);
		float m[4][4];
	};

	struct Quaternion{
		Quaternion(){ v = Vector3f(0.0f, 0.0f, 0.0f); w = 1.0f; }
		Quaternion(const Vector3f& vv, float ww) :v(vv), w(ww){};
		Quaternion(const Mat4 &m);
		Quaternion &operator~(){ //conjugate
			v = -v;
			return *this;
		}
		Quaternion &operator+=(const Quaternion& q){
			v += q.v;
			w += q.w;
			return *this;
		}
		friend Quaternion operator+(const Quaternion &q1, const Quaternion &q2){
			Quaternion rq = q1;
			return rq += q2;
		}
		Quaternion &operator-=(const Quaternion& q){
			v -= q.v;
			w -= q.w;
			return *this;
		}
		friend Quaternion operator-(const Quaternion &q1, const Quaternion &q2){
			Quaternion rq = q1;
			return rq -= q2;
		}
		Quaternion operator*(const Quaternion& q) const{
			Quaternion rq;
			rq.v = v%q.v + w*q.v + q.w*v;
			rq.w = w*q.w - v*q.v;
			return rq;
		}
		Quaternion operator*(float f) const{
			Quaternion rq = *this;
			rq.v *= f;
			rq.w *= f;
			return rq;
		}
		friend inline Quaternion operator*(float f, const Quaternion& q){
			return q*f;
		}

		Quaternion &operator*=(float f){
			v *= f;
			w *= f;
			return *this;
		}
		Quaternion operator/(float f) const{
			Quaternion rq = *this;
			rq.v /= f;
			rq.w /= f;
			return rq;
		}
		Quaternion &operator/=(float f){
			v /= f;
			w /= f;
			return *this;
		}

		Vector3f v;
		float w;
	};

	class SparseVector {
	public:
		SparseVector(): indices(NULL), values(NULL), size(0), capability(0) {}
		SparseVector(const SparseVector& vec) = delete;
		SparseVector& operator=(const SparseVector& vec) = delete;
		SparseVector(SparseVector&& vec);
		SparseVector& operator=(SparseVector&& vec);
		double operator*(const FastSparseVector &vec) const;
		void Set(int index, double val){
			for (size_t i = 0; i < size; i++) {
				if (index == indices[i]) {
					values[i] = val;
					return;
				}
			}
			emplaceBack(index, val);
		}

		void Add(int index, double val){
			for (size_t i = 0; i < size; i++) {
				if (index == indices[i]) {
					values[i] += val;
					return;
				}
			}
			emplaceBack(index, val);
		}

		void emplaceBack(int index, double val){
			if (size >= capability) enlargeVector();
			indices[size] = index; values[size] = val;
			size += 1;
		}

		size_t Size() const { return size; }

		template<class indexPtr, class valPtr> struct IndexValueIterator {
			using indexPointer = indexPtr;
			using valuePointer = valPtr;

			IndexValueIterator(indexPointer idexptr, valuePointer valptr): indexIterator(idexptr), valueIterator(valptr) {}
			IndexValueIterator(const IndexValueIterator&) = default;
			IndexValueIterator& operator=(IndexValueIterator&) = default;
			IndexValueIterator(IndexValueIterator&&) = default;
			IndexValueIterator& operator=(IndexValueIterator&&) = default;

			IndexValueIterator& operator++() {
				++indexIterator; ++valueIterator;
				return *this;
			}
			IndexValueIterator operator++(int) {
				IndexValueIterator temp = *this;
				++(*this);
				return temp;
			}
			IndexValueIterator& operator--() {
				--indexIterator; --valueIterator;
				return *this;
			}
			IndexValueIterator operator--(int) {
				IndexValueIterator temp = *this;
				--(*this);
				return temp;
			}

			bool operator==(const IndexValueIterator& x) const { return indexIterator == x.indexIterator && valueIterator == x.valueIterator; }
			bool operator!=(const IndexValueIterator& x) const { return indexIterator != x.indexIterator || valueIterator != x.valueIterator;}

			indexPointer indexIterator;
			valuePointer valueIterator;
		};

		using IndexValConstIter = IndexValueIterator<const int*, const double*>;
		using IndexValIter = IndexValueIterator<int*, double*>;
		IndexValConstIter cbegin() const { return IndexValConstIter(indices, values); }
		IndexValConstIter cend() const { return IndexValConstIter(indices + size, values + size); }
		IndexValIter begin() { return IndexValIter(indices, values); }
		IndexValIter end() { return IndexValIter(indices + size, values + size); }
		void Clear(){ size = 0; }
		void Reserve(size_t count) { 
			if (count > capability) {
				capability = count;
				reallocSpace(count);
			}
		}
		~SparseVector() {
			free(indices);
			free(values);
		}
	private:
		void reallocSpace(size_t capability);
		void enlargeVector() {
			capability = std::max(size_t(1), capability + capability / 2 + 1);
			reallocSpace(capability);
		}


		int *indices;
		double *values;

		size_t size;
		size_t capability;
		friend class FastSparseVector;
	};

	class FastSparseVector {
	public:
		FastSparseVector(): indexIndicators(NULL), column(0) {}
		FastSparseVector(int col);
		FastSparseVector(const FastSparseVector&) = delete;
		FastSparseVector& operator=(const FastSparseVector&) = delete;
		FastSparseVector(FastSparseVector&& vec);
		FastSparseVector& operator=(FastSparseVector&& vec);
		double operator*(const SparseVector &vec) const;
		void Set(int index, double val) {
			int indicator = indexIndicators[index];
			if (indicator != 0)
				vector.values[indicator] = val;
			else
				emplaceBack(index, val);
		}

		void Add(int index, double val) {
			int indicator = indexIndicators[index];
			if (indicator != 0)
				vector.values[indicator] += val;
			else
				emplaceBack(index, val);
		}

		void Clear() {
			for (size_t i = 1; i < vector.size; i++)
				indexIndicators[vector.indices[i]] = 0;
			vector.size = 1;
		}

		double operator[](int index) const { 
			return vector.values[indexIndicators[index]];
		}

		void emplaceBack(int index, double val) {
			indexIndicators[index] = vector.size;
			vector.indices[vector.size] = index;
			vector.values[vector.size] = val;
			vector.size += 1;
		}

		size_t Size() const { return vector.Size() - 1; }

		using IndexValConstIter = SparseVector::IndexValConstIter;
		using IndexValIter = SparseVector::IndexValIter;
		IndexValConstIter cbegin() const { return ++vector.cbegin(); }
		IndexValConstIter cend() const { return vector.cend(); }
		IndexValIter begin() { return ++vector.begin(); }
		IndexValIter end() { return vector.end(); }

		~FastSparseVector() {
			delete[] indexIndicators;
		}
	private:
		SparseVector vector;

		int *indexIndicators;
		int column;
	};

	/*template<size_t N, size_t...index> struct IndexSequenceGenerator : public IndexSequenceGenerator <N - 1, N - 1, index... > {};
	template<size_t... index> struct IndexSequenceGenerator<0, index...> {};
	template<size_t...index> using IndexSequence = IndexSequenceGenerator<0, index...>;
	template<class T, class FUN, size_t... index> std::array<T, sizeof...(index)> generateSequence(FUN f, IndexSequence<index...> seqs){
		return {  f(index)...  };
	}

	//can be further optimized with constexpr
	template<class T, int N, class FUN> std::array<T, N> generateSequence(FUN f){
		return generateSequence<T>(f, IndexSequenceGenerator<N>{});
	}*/
	

	//inline funtion for vec
	template<class FT> inline VectorBase<FT> Normalize(const VectorBase<FT> &v){
		return v / v.length();
	}


	//inline function for quaternion
	inline float Dot(const Quaternion& q1, const Quaternion& q2){
		return q1.v*q2.v + q1.w*q2.w;
	}

	inline Quaternion Normalize(const Quaternion& q){
		return q / sqrtf(Dot(q, q));
	}

	template<class FT> inline void coordinateSystem(const VectorBase<FT> &v0, VectorBase<FT> &v1, VectorBase<FT> &v2){
		if (fabs(v0.x) > fabs(v0.y)){
			FT invlen = FT(1.0) / sqrt(v0.x * v0.x + v0.z * v0.z);
			v1 = VectorBase<FT>(-v0.z * invlen, FT(0.0), v0.x * invlen);
		}
		else{
			FT invlen = 1.f / sqrt(v0.y * v0.y + v0.z * v0.z);
			v1 = VectorBase<FT>(FT(0.0), -v0.z * invlen, v0.y * invlen);
		}
		v2 = v0 % v1;
	}

	template<class FT> inline Tensor2<FT> VectorBase<FT>::operator^(const VectorBase<FT>& v) const{
		Tensor2<FT> ret;
		for (int i = 0; i < 3; i++){
			for (int j = 0; j < 3; j++){
				ret(i, j) = (*this)[i] * v[j];
			}
		}
		return ret;
	}

	inline Mat4 lookAt(const Vector3f& eye, const Vector3f& at, const Vector3f& up){
		Vector3f n = Normalize(eye - at);
		Vector3f u = Normalize(up%n);
		Vector3f v = Normalize(n%u);

		Vector3f trans = -eye;

		float xTrans = u*trans;
		float yTrans = v*trans;
		float zTrans = n*trans;

		Mat4 m = Mat4(u.x, u.y, u.z, xTrans,
			v.x, v.y, v.z, yTrans,
			n.x, n.y, n.z, zTrans,
			0.f, 0.f, 0.f, 1.f);
		return m;
	}

	inline Mat4 Orient(float left, float right, float bottom, float top, float znear, float zfar){
		float invrdl = 1.f / (right - left);
		float invtdb = 1.f / (top - bottom);
		float invfdn = 1.f / (zfar - znear);
		return  Mat4(2.f * invrdl, 0.f, 0.f, -(left + right) * invrdl,
			0.f, 2.f * invtdb, 0.f, -(top + bottom) * invtdb,
			0.f, 0.f, -2.f * invfdn, -(zfar + znear) * invfdn,
			0.f, 0.f, 0.f, 1.f);
	}



	inline double SparseVector::operator*(const FastSparseVector &vec) const {
		double dot = 0;
		for (size_t i = 0; i < size; i++)
			dot += values[i] * vec[indices[i]];

		return dot;
	}

	inline double FastSparseVector::operator*(const SparseVector &vec) const {
		return vec * (*this);
	}

}

#endif