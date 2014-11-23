#if defined(_MSC_VER)
#pragma once
#endif

#ifndef ODER_CORE_LATOOL_H
#define ODER_CORE_LATOOL_H

#include "oder.h"

namespace ODER{
	template<class FT> struct VectorBase{
		VectorBase(){ x = y = z = FT(0.0); }
		VectorBase(FT xx, FT yy, FT zz) :x(xx), y(yy), z(zz){}
		VectorBase(const VectorBase &v) :x(v.x), y(v.y), z(v.z){}

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


		inline Tensor2<FT> operator^(const VectorBase& v) const;

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

	struct Vector : public VectorBase <float> {
		Vector() : VectorBase(){}
		Vector(float xx, float yy, float zz) :VectorBase(xx, yy, zz){}
		Vector(const Vector &v) : VectorBase(v){}
		Vector(const VectorBase<float> &v) : VectorBase(v){}
		Vector(const VectorBase<double> &v){
			x = float(v.x); y = float(v.y); z = float(v.z);
		}
	};

	template<class FT> struct Tensor2{
		Tensor2(){
			m[0][0] = m[0][1] = m[0][2] =
			m[1][0] = m[1][1] = m[1][2] = 
			m[2][0] = m[2][1] = m[2][2] = 0.0;
		}
		Tensor2(FT mat[3][3]){
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
		friend Tensor2 Inverse(const Tensor2&){
			Tensor2 r;

			FT det = tensor.Determinant();

			Assert(fabsf(det) > 1e-7);

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
		friend Tensor2 Transpose(const Tensor2&){
			Tensor2 ret;
			for (int i = 0; i < 3; i++){
				for (int j = 0; j < 3; j++){
					ret.m[i][j] = m.m[j][i];
				}
			}
			return ret;
		}
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

	inline void coordinateSystem(const Vector &v0, Vector *v1, Vector *v2){
		if (fabs(v0.x) > fabs(v0.y)){
			float invlen = 1.f / sqrtf(v0.x*v0.x + v0.z*v0.z);
			*v1 = Vector(-v0.z*invlen, 0.f, v0.x*invlen);
		}
		else{
			float invlen = 1.f / sqrtf(v0.y*v0.y + v0.z*v0.z);
			*v1 = Vector(0.f, -v0.z*invlen, v0.y*invlen);
		}
		*v2 = v0 % *v1;
	}

	template<class FT> inline Tensor2<FT> VectorBase<FT>::operator^(const VectorBase<FT>& v) const{
		Tensor2<FT> ret;
		for (int i = 0; i < 3; i++){
			for (int j = 0; j < 3; j++){
				ret.m[i][j] = (*this)[i] * v[j];
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

	template<class FT> void gaussianElimination3x3(FT *A, FT *rhs, FT *result){
		//a.k.a LU decomposition
		int rowPermute[3] = { 0, 1, 2 };
		int columnPermute[3] = { 0, 1, 2 };

		constexpr int order = 3;
		constexpr int count = order - 1;
		//forward substitution
		for (int i = 0; i < count; i++){
			FT largest = 0.0;
			int rowIndex = i;
			int columnIndex = i;
			for (int j = i; j < order; j++){
				int subRow = rowPermute[j];
				for (int k = i; k < order; k++){
					int subColumn = columnPermute[k];
					FT entry = fabs(A[order*subRow + subColumn]);
					if (entry > largest){
						largest = entry;
						rowIndex = j;
						columnIndex = k;
					}
				}
			}
			Assert(largest > 0.0);
			if (rowPermute[i] != rowPermute[rowIndex])
				std::swap(rowPermute[i], rowPermute[rowIndex]);
			if (columnPermute[i] != columnPermute[columnIndex])
				std::swap(columnPermute[i], columnPermute[columnIndex]);

			int row = rowPermute[i];
			int column = columnPermute[i];
			FT scale = FT(1.0) / A[order*row + column];
			for (int j = i + 1; j < order; j++){
				int subRow = rowPermute[j];
				FT lower = scale*A[order*subRow + column];
				for (int k = i + 1; k < order; k++){
					int subColumn = columnPermute[k];
					A[order*subRow + subColumn] -= lower * A[order*row + subColumn];
				}
				rhs[subRow] -= lower * rhs[row];
			}
		}

		//backward substitution
		for (int i = count; i >= 0; i--){
			int row = rowPermute[i];
			int column = columnPermute[i];
			FT dot = 0.0;
			for (int j = i + 1; j < order; j++){
				int subColumn = columnPermute[j];
				dot += A[order*row + subColumn] * result[subColumn];
			}
			result[column] = (rhs[row] - dot) / A[order*row + column];
		}
	}

	template<class FUN, class FT> FT findRoot(FUN f, FT start, FT end, FT tolerance = 0.0, FT epsilon = FT(2e-8)){
		FT a = start, b = end, c = start, d = 0.0, e = 0.0;
		FT fa = f(a), fb = f(b), fc = fa;

		constexpr int maxIter = 1000;
		for (int i = 0; i < maxIter; i++){
			if (fb > 0.0 && fc > 0.0 || fb < 0.0 && fc < 0.0){
				c = a; fc = fa;
				d = e = b - a;
			}
			if (fabs(fc) < fabs(fb)){
				a = b; b = c; c = a;
				fa = fb; fb = fc; fc = fa;
			}

			FT tol = FT(2.0)*FT(epsilon)*fabs(b) + tolerance;
			FT m = FT(0.5)*(c - b);
			
			if (fabs(m) <= tol || fb == FT(0.0)) return b;
			if (fabs(e) >= tol && fa > fb){
				//inverse quadratic
				FT s = fa / fb;
				FT p = 0.0, q = 0.0;
				if (a == c){
					p = FT(2.0)*m*s;
					q = FT(1.0) - s;
				}
				else{
					FT t = fa / fc;
					FT r = fb / fc;
					p = s * (FT(2.0)*t*(r - t)*m - (FT(1.0) - r)*(b - a));
					q = (t - FT(1.0))*(r - FT(1.0))*(s - FT(1.0));
				}
				if (p > 0.0) q = -q;
				else p = -p;

				if (FT(2.0)*p < min(FT(3.0)*q*m - fabs(tol*q), fabs(e*q))){
					e = d;
					d = p / q;
				}
				else{
					//bisection
					d = m;
					e = d;
				}
			}
			else{
				//bisection
				d = m;
				e = d;
			}
			a = b;
			fa = fb;
			if (fabs(d) > tol)
				b += d;
			else
				b += (m > 0 ? tol : -tol);
			fb = f(b);
		}
		Severe("max iteration in ODER::findRoot");
		return 0.0;
	}
}

#endif