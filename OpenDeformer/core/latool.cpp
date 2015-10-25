#include "stdafx.h"
#include "latool.h"
#include "memory.h"

namespace ODER{
	Mat4::Mat4(float mat[4][4]){
		memcpy(m, mat, 16 * sizeof(float));
	}

	Mat4::Mat4(float m00, float m01, float m02, float m03,
		float m10, float m11, float m12, float m13,
		float m20, float m21, float m22, float m23,
		float m30, float m31, float m32, float m33){
		m[0][0] = m00; m[0][1] = m01; m[0][2] = m02; m[0][3] = m03;
		m[1][0] = m10; m[1][1] = m11; m[1][2] = m12; m[1][3] = m13;
		m[2][0] = m20; m[2][1] = m21; m[2][2] = m22; m[2][3] = m23;
		m[3][0] = m30; m[3][1] = m31; m[3][2] = m32; m[3][3] = m33;
	}

	Mat4 Transpose(const Mat4& mat){
		Mat4 m;
		for (int i = 0; i < 4; i++){
			for (int j = 0; j < 4; j++){
				m.m[i][j] = mat.m[j][i];
			}
		}
		return m;
	}

	Mat4 Inverse(const Mat4& mat){
		Mat4 m;
		float scaleSquareInverse = 1.0f;
		float inverse_mat_ij = 0.0;
		for (int j = 0; j < 3; j++){
			scaleSquareInverse = 1.0f / (mat.m[0][j] * mat.m[0][j] + mat.m[1][j] * mat.m[1][j] + mat.m[2][j] * mat.m[2][j]);
			for (int i = 0; i < 3; i++){
				inverse_mat_ij = mat.m[i][j] * scaleSquareInverse;
				m.m[j][i] = inverse_mat_ij;//inverse rotation and scale
				m.m[j][3] -= inverse_mat_ij*mat.m[i][3];//inverse translate
			}
		}
		return m;
	}

	Mat4::Mat4(const Quaternion &q){
		float xx = q.v.x*q.v.x, yy = q.v.y*q.v.y, zz = q.v.z*q.v.z;
		float xy = q.v.x*q.v.y, xz = q.v.x*q.v.z, yz = q.v.y*q.v.z;
		float wx = q.w*q.v.x, wy = q.w*q.v.y, wz = q.w*q.v.z;

		m[0][0] = 1.0f - 2.0f*(yy + zz); m[0][1] = 2.0f*(xy + wz);      m[0][2] = 2.0f*(xz - wy); m[0][3] = 0.0f;
		m[1][0] = 2.0f*(xy - wz);      m[1][1] = 1.0f - 2.0f*(xx + zz); m[1][2] = 2.0f*(yz - wx); m[1][3] = 0.0f;
		m[2][0] = 2.0f*(xz + wy);      m[2][1] = 2.0f*(yz - wx);      m[2][2] = 1.0f - 2.0f*(xx + yy); m[2][3] = 0.0f;
		m[3][0] = 0.0f; m[3][1] = 0.0f; m[3][2] = 0.0f; m[3][3] = 1.0f;
	}

	Quaternion::Quaternion(const Mat4 &m){
		float trace = m.m[0][0] + m.m[0][1] + m.m[0][2];
		if (trace > 0.0f){
			float s = sqrtf(trace + 1.0f);
			w = s / 2.0f;
			s = 0.5f / s;
			v.x = (m.m[2][1] - m.m[1][2])*s;
			v.y = (m.m[0][2] - m.m[2][0])*s;
			v.z = (m.m[1][0] - m.m[0][1])*s;
		}
		else{
			int i = 0;
			if (m.m[1][1] > m.m[0][0]) i = 1;
			if (m.m[2][2] > m.m[i][i]) i = 2;
			float s = 0.f;
			switch (i){
			case 0:
				s = sqrtf(m.m[0][0] - m.m[1][1] - m.m[2][2] + 1.0f);
				v.x = s*0.5f;
				if (s != 0.0f) s = 0.5f / s;
				w = (m.m[1][2] - m.m[2][1])*s;
				v.y = (m.m[0][1] + m.m[1][0])*s;
				v.z = (m.m[2][0] + m.m[0][2])*s;
				break;

			case 1:
				s = sqrtf(m.m[1][1] - m.m[0][0] - m.m[2][2] + 1.0f);
				v.y = s*0.5f;
				if (s != 0.0f) s = 0.5f / s;
				w = (m.m[2][0] - m.m[0][2])*s;
				v.x = (m.m[0][1] + m.m[1][0])*s;
				v.z = (m.m[2][1] + m.m[1][2])*s;
				break;

			case 2:
				s = sqrtf(m.m[2][2] - m.m[1][1] - m.m[0][0] + 1.0f);
				v.z = s*0.5f;
				if (s != 0.0f) s = 0.5f / s;
				w = (m.m[0][1] - m.m[1][0])*s;
				v.y = (m.m[2][0] + m.m[0][2])*s;
				v.z = (m.m[1][2] + m.m[2][1])*s;
				break;
			}
		}
	}

}