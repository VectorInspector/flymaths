#include "mat4x4.h"
#include <cmath>
#include <algorithm>

template <typename X>
Mat4x4 <X>& Mat4x4 <X>::Transpose () {
	using std::swap;
	swap(u.y, v.x);
	swap(u.z, w.x);
	swap(u.w, x.x);
	swap(v.z, w.y);
	swap(v.w, x.y);
	swap(w.w, x.z);
	return *this;
}

// Matrix component addition
template <typename X>
Mat4x4 <X> Mat4x4 <X>::operator+ (const Mat4x4 <X>& m) const {
	return Mat4x4 <X> (u + m.u, v + m.v, w + m.w, x + m.x);
}

// Matrix x Matrix multiplication
template <typename X>
Mat4x4 <X> Mat4x4 <X>::operator* (const Mat4x4 <X>& m) const {
	return Mat4x4 <X> (
	Vec4 <X> (
		u.x * m.u.x + v.x * m.u.y + w.x * m.u.z + x.x * m.u.w,
		u.y * m.u.x + v.y * m.u.y + w.y * m.u.z + x.y * m.u.w,
		u.z * m.u.x + v.z * m.u.y + w.z * m.u.z + x.z * m.u.w,
		u.w * m.u.x + v.w * m.u.y + w.w * m.u.z + x.w * m.u.w),
	Vec4 <X> (
		u.x * m.v.x + v.x * m.v.y + w.x * m.v.z + x.x * m.v.w,
		u.y * m.v.x + v.y * m.v.y + w.y * m.v.z + x.y * m.v.w,
		u.z * m.v.x + v.z * m.v.y + w.z * m.v.z + x.z * m.v.w,
		u.w * m.v.x + v.w * m.v.y + w.w * m.v.z + x.w * m.v.w),
	Vec4 <X> (
		u.x * m.w.x + v.x * m.w.y + w.x * m.w.z + x.x * m.w.w,
		u.y * m.w.x + v.y * m.w.y + w.y * m.w.z + x.y * m.w.w,
		u.z * m.w.x + v.z * m.w.y + w.z * m.w.z + x.z * m.w.w,
		u.w * m.w.x + v.w * m.w.y + w.w * m.w.z + x.w * m.w.w),
	Vec4 <X> (
		u.x * m.x.x + v.x * m.x.y + w.x * m.x.z + x.x * m.x.w,
		u.y * m.x.x + v.y * m.x.y + w.y * m.x.z + x.y * m.x.w,
		u.z * m.x.x + v.z * m.x.y + w.z * m.x.z + x.z * m.x.w,
		u.w * m.x.x + v.w * m.x.y + w.w * m.x.z + x.w * m.x.w));
}

template <typename X>
Mat4x4 <X>& Mat4x4 <X>::operator+= (const Mat4x4 <X>& m) {
	u += m.u;
	v += m.v;
	w += m.w;
	x += m.x;
	return *this;
}

template <typename X>
Mat4x4 <X>& Mat4x4 <X>::operator-= (const Mat4x4 <X>& m) {
	u -= m.u;
	v -= m.v;
	w -= m.w;
	x -= m.x;
	return *this;
}

template <typename X>
Mat4x4 <X>& Mat4x4 <X>::operator*= (const Mat4x4 <X>& m) {
	*this = (*this) * m;
	return *this;
}

template <typename X>
Mat4x4 <X> Mat4x4 <X>::operator* (double f) const {
	return Mat4x4 <X> (f * u, f * v, f * w, f * x);
}

template <typename X>
Mat4x4 <X> Mat4x4 <X>::operator/ (double f) const {
	return Mat4x4 <X> (u / f, v / f, w / f, x / f);
}

// Component scaling
template <typename X>
Mat4x4 <X>& Mat4x4 <X>::operator*= (double f) {
	u *= f;
	v *= f;
	w *= f;
	x *= f;
	return *this;
}

// Component inverse scaling
template <typename X>
Mat4x4 <X>& Mat4x4 <X>::operator/ (double f) {
	u /= f;
	v /= f;
	w /= f;
	x /= f;
	return *this;
}

template <typename X>
X Mat4x4 <X>::Trace () const {
	return u.x + v.y + w.z + x.w;
}

template <typename X>
bool Mat4x4 <X>::DiagonalDominant () const {
	return
		std::abs(v.x) + std::abs(w.x) + std::abs(x.x) <= std::abs(u.x) &&
		std::abs(v.y) + std::abs(w.y) + std::abs(x.y) <= std::abs(v.y) &&
		std::abs(v.z) + std::abs(w.z) + std::abs(x.z) <= std::abs(w.z) &&
		std::abs(v.w) + std::abs(w.w) + std::abs(x.w) <= std::abs(x.w);
}

template <typename X>
bool Mat4x4 <X>::StrictDiagonalDominant () const {
	return
		std::abs(v.x) + std::abs(w.x) + std::abs(x.x) < std::abs(u.x) &&
		std::abs(v.y) + std::abs(w.y) + std::abs(x.y) < std::abs(v.y) &&
		std::abs(v.z) + std::abs(w.z) + std::abs(x.z) < std::abs(w.z) &&
		std::abs(v.w) + std::abs(w.w) + std::abs(x.w) < std::abs(x.w);
}

template <typename X>
Mat4x4 <X> operator* (double f, Mat4x4 <X> m) {
	return Mat4x4 <X> (f * m.u, f * m.v, f * m.w, f * m.x);
}

template <typename X>
Mat4x4 <X> operator/ (double f, Mat4x4 <X> m) {
	return Mat4x4 <X> (m.u / f, m.v / f, m.w / f, m.x / f);
}

template <typename X>
Mat4x4 <X> IdentityMat4x4 () {
	return Mat4x4 <X> (Vec4 <X> (1, 0, 0, 0), Vec4 <X> (0, 1, 0, 0), Vec4 <X> (0, 0, 1, 0), Vec4 <X> (0, 0, 0, 1));
}

template <typename X>
Mat4x4 <X> ZeroMat4x4 () {
	return Mat4x4 <X> (Vec4 <X> (0, 0, 0, 0), Vec4 <X> (0, 0, 0, 0), Vec4 <X> (0, 0, 0, 0), Vec4 <X> (0, 0, 0, 0));
}

template <typename X>
Mat4x4 <X> FilledMat4x4 (X e) {
	return Mat4x4 <X> (Vec4 <X> (e, e, e, e), Vec4 <X> (e, e, e, e), Vec4 <X> (e, e, e, e), Vec4 <X> (e, e, e, e));
}
