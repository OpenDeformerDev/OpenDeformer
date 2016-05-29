#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_CORE_LATOOL_H
#define ODER_CORE_LATOOL_H

#include "oder.h"
#include "datastructure.h"
#include <array>

namespace ODER{
	template<class FT> struct VectorBase{
		VectorBase(){ x = y = z = FT(0.0); }
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
			return _isnan(x) || _isnan(y) || _isnan(z);
		}
		FT x, y, z;
	};

	template<class FT> struct Tensor2{
		Tensor2(){
			m[0][0] = m[0][1] = m[0][2] =
			m[1][0] = m[1][1] = m[1][2] = 
			m[2][0] = m[2][1] = m[2][2] = 0.0;
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
		Vector operator*(const Vector& v) const{
			float x = v.x, y = v.y, z = v.z;
			return Vector(m[0][0] * x + m[0][1] * y + m[0][2] * z, m[1][0] * x + m[1][1] * y + m[1][2] * z, m[2][0] * x + m[2][1] * y + m[2][2] * z);
		}
		friend Mat4 Transpose(const Mat4&);
		friend Mat4 Inverse(const Mat4&);
		float m[4][4];
	};

	struct Quaternion{
		Quaternion(){ v = Vector(0.0f, 0.0f, 0.0f); w = 1.0f; }
		Quaternion(const Vector& vv, float ww) :v(vv), w(ww){};
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

		Vector v;
		float w;
	};

	//not thread safe
	class SparseVector{
	public:
		using IndexValPair = std::pair<int, double>;
		using IndexValConstIter = RecycledList<IndexValPair>::const_iterator;
		using IndexValIter = RecycledList<IndexValPair>::iterator;
		SparseVector(){}
		SparseVector(const SparseVector&) = delete;
		SparseVector& operator=(const SparseVector&) = delete;
		SparseVector(SparseVector&&) = default;
		SparseVector& operator=(SparseVector&&) = default;
		double operator*(const SparseVector& vec) const;
		double operator*(const double *vec) const;
		void Set(int index, double val){
			auto iter = std::lower_bound(indices.begin(), indices.end(), index, 
				[](const IndexValPair& lhs, int rhs){ return lhs.first < rhs; });
			if (iter == indices.end() || iter->first != index)
				indices.insert(iter, IndexValPair(index, val));
			else
				iter->second = val;
		}
		void Add(int index, double val){
			auto iter = std::lower_bound(indices.begin(), indices.end(), index, 
				[](const IndexValPair& lhs, int rhs){ return lhs.first < rhs; });
			if (iter == indices.end() || iter->first != index)
				indices.insert(iter, IndexValPair(index, val));
			else
				iter->second += val;
		}
		IndexValIter Set(const IndexValIter& pos, int index, double val){
			return indices.insert(pos, IndexValPair(index, val));
		}
		IndexValConstIter Set(const IndexValConstIter& pos, int index, double val){
			return indices.insert(pos, IndexValPair(index, val));
		}
		void emplaceBack(int index, double val){
			indices.emplace_back(index, val);
		}
		IndexValConstIter Delete(const IndexValConstIter& iter){ return indices.erase(iter); }
		IndexValIter Delete(const IndexValIter& iter){ return indices.erase(iter); }
		IndexValConstIter cbegin() const { return indices.cbegin(); }
		IndexValConstIter cend() const { return indices.cend(); }
		IndexValIter begin() { return indices.begin(); }
		IndexValIter end() { return indices.end(); }
		void Clear(){ indices.clear(); }
	private:
		RecycledList<IndexValPair> indices;
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

	inline Mat4 lookAt(const Vector& eye, const Vector& at, const Vector& up){
		Vector n = Normalize(eye - at);
		Vector u = Normalize(up%n);
		Vector v = Normalize(n%u);

		Vector trans = -eye;

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

	inline double SparseVector::operator*(const SparseVector& vec) const{
		auto rightIter = vec.cbegin();
		auto rightEnd = vec.cend();
		auto leftIter = cbegin();
		auto leftEnd = cend();

		double dot = 0.0;
		while (rightIter != rightEnd && leftIter != leftEnd){
			int gap = rightIter->first - leftIter->first;
			if (gap == 0){
				dot += rightIter->second * leftIter->second;
				++rightIter;
				++leftIter;
			}
			else if (gap > 0)
				++leftIter;
			else
				++rightIter;
		}

		return dot;
	}

	inline double SparseVector::operator*(const double* vec) const {
		auto iter = cbegin();
		auto end = cend();

		double dot = 0.0;
		while (iter != end) {
			dot += vec[iter->first] * iter->second;
			++iter;
		}

		return dot;
	}

	inline double operator*(const double *left, const SparseVector& right) {
		return right * left;
	}

	template<class FT> inline VectorBase<FT> triangleNormal(const VectorBase<FT>& a,
		const VectorBase<FT>& b, const VectorBase<FT>& c) {
		VectorBase<FT> ab = b - a;
		VectorBase<FT> ac = c - a;
		VectorBase<FT> bc = c - b;
		FT abLen2 = ab.length2();
		FT acLen2 = ac.length2();
		FT bcLen2 = bc.length2();

		VectorBase<FT> *u = &ab;
		VectorBase<FT> *v = &ac;

		if (bcLen2 < abLen2) {
			if (abLen2 < acLen2) v = &bc;
			else {
				u = &ac; v = &bc;
			}
		}
		else if (bcLen2 < acLen2) v = &bc;

		return *u % *v;
	}


	constexpr float const_pow(float base, int exp) noexcept{
		return exp < 0 ? 1.f / const_pow(base, -exp) : 
			(exp == 0 ? 1.f :
				exp % 2 == 0 ? const_pow(base * base, exp / 2) :
				const_pow(base * base, (exp - 1) / 2) * base);
	}

	constexpr double const_pow(double base, int exp) noexcept{
		return exp < 0 ? 1.0 / const_pow(base, -exp) : 
			(exp == 0 ? 1.0 :
				exp % 2 == 0 ? const_pow(base * base, exp / 2) :
				const_pow(base * base, (exp - 1) / 2) * base);
	}


}

#endif