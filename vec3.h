#ifndef VEC3_H
#define VEC3_H

#include "space.h"
#include <cmath>
#include <algorithm>

// It's a drop-in header!
template <typename X>
struct Vec3 {
	Vec3 <X> Abs () const;
	X Max () const;
	X Min () const;
	X Dot (const Vec3 <X>& v) const;
	X LenSquared () const;
	double Len () const;
	Vec3& Unit0 ();
	Vec3 <X> Cross (const Vec3 <X>& v) const;
	Vec3 <X> operator- () const;
	Vec3 <X>& operator+= (const Vec3 <X>& v);
	Vec3 <X>& operator-= (const Vec3 <X>& v);
	Vec3 <X> operator+ (const Vec3 <X>& v) const;
	Vec3 <X> operator- (const Vec3 <X>& v) const;
	
	X x;
	X y;
	X z;
};

template <typename X>
Vec3 <X> RadVec3 (double yaw_rad, double pitch_rad);

template <typename X>
Vec3 <X> DegVec3 (double yaw_deg, double pitch_deg);

template <typename X>
Vec3 <X> operator* (double f, Vec3 <X> v);

template <typename X>
Vec3 <X> operator* (Vec3 <X> v, double f);

template <typename X>
Vec3 <X> operator/ (Vec3 <X> v, double f);

template <typename X>
Vec3 <X>& operator*= (Vec3 <X>& v, double f);

template <typename X>
Vec3 <X>& operator/= (Vec3 <X>& v, double f);

#include <cmath>
#include <algorithm>

template <typename X>
Vec3 <X> Vec3 <X>::Abs () const {
	return Vec3 <X> (std::abs(x), std::abs(y), std::abs(z));
}

template <typename X>
X Vec3 <X>::Max () const {
	return std::max(std::max(x, y), z);
}

template <typename X>
X Vec3 <X>::Min () const {
	return std::min(std::min(x, y), z);
}

template <typename X>
X Vec3 <X>::Dot (const Vec3 <X>& v) const {
	return x * v.x + y * v.y + z * v.z;
}

template <typename X>
X Vec3 <X>::LenSquared () const {
	return x * x + y * y + z * z;
}

template <typename X>
double Vec3 <X>::Len () const {
	return std::sqrt(x * x + y * y + z * z);
}

template <typename X>
Vec3 <X>& Vec3 <X>::Unit0 () {
	auto len = Len();
	
	if(0 < len) {
		x /= len;
		y /= len;
		z /= len;
	}
	
	else {
		x = 0;
		y = 0;
		z = 0;
	}
	
	return *this;
}

template <typename X>
Vec3 <X> Vec3 <X>::Cross (const Vec3 <X>& v) const {
	return Vec3 <X> (
		y * v.z - z * v.y,
		z * v.x - x * v.z,
		x * v.y - y * v.x);
}

template <typename X>
Vec3 <X> Vec3 <X>::operator- () const {
	return Vec3 <X> (-x, -y, -z);
}

template <typename X>
Vec3 <X>& Vec3 <X>::operator+= (const Vec3 <X>& v) {
	x += v.x;
	y += v.y;
	z += v.z;
	return *this;
}

template <typename X>
Vec3 <X>& Vec3 <X>::operator-= (const Vec3 <X>& v) {
	x -= v.x;
	y -= v.y;
	z -= v.z;
	return *this;
}

template <typename X>
Vec3 <X> Vec3 <X>::operator+ (const Vec3 <X>& v) const {
	return Vec3 <X> (x + v.x, y + v.y, z + v.z);
}

template <typename X>
Vec3 <X> Vec3 <X>::operator- (const Vec3 <X>& v) const {
	return Vec3 <X> (x - v.x, y - v.y, z - v.z);
}


template <typename X>
Vec3 <X> RadVec3 (double yaw_rad, double pitch_rad) {
	auto cos_yaw = std::cos(yaw_rad);
	auto sin_yaw = std::sin(yaw_rad);
	auto cos_pitch = std::cos(pitch_rad);
	auto sin_pitch = std::sin(pitch_rad);
	
	return Vec3 <X> (cos_pitch * cos_yaw, cos_pitch * sin_yaw, sin_pitch);
}

template <typename X>
Vec3 <X> DegVec3 (double yaw_deg, double pitch_deg) {
	auto yaw_rad = DegToRad(yaw_deg);
	auto pitch_rad = DegToRad(pitch_deg);
	auto cos_yaw = std::cos(yaw_rad);
	auto sin_yaw = std::sin(yaw_rad);
	auto cos_pitch = std::cos(pitch_rad);
	auto sin_pitch = std::sin(pitch_rad);
	
	return Vec3 <X> (cos_pitch * cos_yaw, cos_pitch * sin_yaw, sin_pitch);
}

template <typename X>
Vec3 <X> operator* (double f, Vec3 <X> v) {
	return Vec3 <X> (f * v.x, f * v.y, f * v.z);
}

template <typename X>
Vec3 <X> operator* (Vec3 <X> v, double f) {
	return f * v;
}

template <typename X>
Vec3 <X> operator/ (Vec3 <X> v, double f) {
	return Vec3 <X> (v.x / f, v.y / f, v.z / f);
}

template <typename X>
Vec3 <X>& operator*= (Vec3 <X>& v, double f) {
	v.x *= f;
	v.y *= f;
	v.z *= f;
	return v;
}

template <typename X>
Vec3 <X>& operator/= (Vec3 <X>& v, double f) {
	v.x /= f;
	v.y /= f;
	v.z /= f;
	return v;
}

using Vec3f = Vec3 <float>;
using Vec3d = Vec3 <double>;

#endif
