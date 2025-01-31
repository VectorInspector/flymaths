#include "vec4.h"
#include "space.h"

template <typename X>
Vec4 <X> Vec4 <X>::Abs () const {
	return Vec4 <X> (abs(x), abs(y), abs(z), abs(w));
}

template <typename X>
X Vec4 <X>::Max () const {
	return max(max(max(x, y), z), w);
}

template <typename X>
X Vec4 <X>::Min () const {
	return min(min(min(x, y), z), w);
}

template <typename X>
X Vec4 <X>::Dot (const Vec4 <X>& v) const {
	return x * v.x + y * v.y + z * v.z + w * v.w;
}

template <typename X>
X Vec4 <X>::LenSquared () const {
	return x * x + y * y + z * z + w * w;
}

template <typename X>
double Vec4 <X>::Len () const {
	return sqrt(x * x + y * y + z * z + w * w);
}

template <typename X>
Vec4 <X>& Vec4 <X>::Unit0 () {
	auto len = Len();
	
	if(0 < len) {
		x /= len;
		y /= len;
		z /= len;
		w /= len;
	}
	
	else {
		x = 0;
		y = 0;
		z = 0;
		w = 0;
	}
	
	return *this;
}

template <typename X>
Vec4 <X> Vec4 <X>::operator- () const {
	return Vec4 <X> (-x, -y, -z, -w);
}

template <typename X>
Vec4 <X>& Vec4 <X>::operator+= (const Vec4 <X>& v) {
	x += v.x;
	y += v.y;
	z += v.z;
	w += v.w;
	return *this;
}

template <typename X>
Vec4 <X>& Vec4 <X>::operator-= (const Vec4 <X>& v) {
	x -= v.x;
	y -= v.y;
	z -= v.z;
	w -= v.w;
	return *this;
}

template <typename X>
Vec4 <X> Vec4 <X>::operator+ (const Vec4 <X>& v) const {
	return Vec4 <X> (x + v.x, y + v.y, z + v.z, w + v.w);
}

template <typename X>
Vec4 <X> Vec4 <X>::operator- (const Vec4 <X>& v) const {
	return Vec4 <X> (x - v.x, y - v.y, z - v.z, w - v.w);
}


template <typename X>
Vec4 <X> RadVec4 (double yaw_rad, double pitch_rad) {
	auto cos_yaw = cos(yaw_rad);
	auto sin_yaw = sin(yaw_rad);
	auto cos_pitch = cos(pitch_rad);
	auto sin_pitch = sin(pitch_rad);
	
	return Vec4 <X> (cos_pitch * cos_yaw, cos_pitch * sin_yaw, sin_pitch);
}

template <typename X>
Vec4 <X> DegVec4 (double yaw_deg, double pitch_deg) {
	auto yaw_rad = DegToRad(yaw_deg);
	auto pitch_rad = DegToRad(pitch_deg);
	auto cos_yaw = cos(yaw_rad);
	auto sin_yaw = sin(yaw_rad);
	auto cos_pitch = cos(pitch_rad);
	auto sin_pitch = sin(pitch_rad);
	
	return Vec4 <X> (cos_pitch * cos_yaw, cos_pitch * sin_yaw, sin_pitch);
}

template <typename X>
Vec4 <X> operator* (double f, Vec4 <X> v) {
	return Vec4 <X> (f * v.x, f * v.y, f * v.z, f * v.w);
}

template <typename X>
Vec4 <X> operator* (Vec4 <X> v, double f) {
	return f * v;
}

template <typename X>
Vec4 <X> operator/ (Vec4 <X> v, double f) {
	return Vec4 <X> (v.x / f, v.y / f, v.z / f, v.w / f);
}

template <typename X>
Vec4 <X>& operator*= (Vec4 <X>& v, double f) {
	v.x *= f;
	v.y *= f;
	v.z *= f;
	v.w *= f;
	return v;
}

template <typename X>
Vec4 <X>& operator/= (Vec4 <X>& v, double f) {
	v.x /= f;
	v.y /= f;
	v.z /= f;
	v.w /= f;
	return v;
}

struct <float> Vec4;
struct <double> Vec4;
