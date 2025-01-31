#ifndef FLY_MATHS_H
#define FLY_MATHS_H

#include <cmath>
#include <algorithm>
#include <vector>
#include <string>
#include <sstream>

// It's a drop-in header!
template <typename X> struct Vec2;
template <typename X> struct Vec3;
template <typename X> struct Vec4;
template <typename X> struct Mat3x3;
template <typename X> struct Mat4x4;
template <typename X> struct Quaternion;
struct CosSinVals;
struct SpaceConst;

using Vec2f = Vec2 <float>;
using Vec2d = Vec2 <double>;
using Vec3f = Vec3 <float>;
using Vec3d = Vec3 <double>;
using Mat3x3f = Mat3x3 <float>;
using Mat3x3d = Mat3x3 <double>;
using Mat4x4f = Mat4x4 <float>;
using Mat4x4d = Mat4x4 <double>;
using Quatf = Quaternion <float>;
using Quatd = Quaternion <double>;

// Space and conversions
// Helper which stores pairs of trigonometric values for convenience.
struct CosSinVals {
	CosSinVals (double yaw, double pitch, double roll);
	CosSinVals (double yaw, double pitch);
	CosSinVals ();
	
	double ca;
	double cv;
	double cr;
	double sa;
	double sv;
	double sr;
};

struct SpaceConst {
	static constexpr double Pi () {
		return std::acos(-1);
	}
	
	static constexpr double PiHalf () {
		return std::acos(0);
	}
	
	static constexpr double RadToDegFactor () {
		return 180.0 / Pi();
	}
	
	static constexpr double DegToRadFactor () {
		return Pi() / 180;
	}
};

double DegToRad (double deg);
double RadToDeg (double rad);

// Conversions back to radiants
struct VectorToRad {
	
	static constexpr double sin_extreme = 0.5 - 1e-8;
	
	template <typename X>
	void FromQuat (Quaternion <X> r);
	double yaw_rad;
	double pitch_rad;
	double roll_rad;
};

// 3D vector.
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

// 4D vector
template <typename X>
struct Vec4 {

	Vec4 <X> Abs () const;
	X Max () const;
	X Min () const;
	X Dot (const Vec4 <X>& v) const;
	X LenSquared () const;
	double Len () const;
	Vec4 <X>& Unit0 ();
	Vec4 <X> operator- () const;
	Vec4 <X>& operator+= (const Vec4 <X>& v);
	Vec4 <X>& operator-= (const Vec4 <X>& v);
	Vec4 <X> operator+ (const Vec4 <X>& v) const;
	Vec4 <X> operator- (const Vec4 <X>& v) const;
	
	X x;
	X y;
	X z;
	X w;
};

// 3x3 Matrix
template <typename X>
struct Mat3x3 {
	Mat3x3& Transpose ();
	
	// Form the matrix of co-determinants.
	Mat3x3& Adjoint ();
	
	bool IsSingular () const;
	
	// Inverts a matrix if possible.
	Mat3x3& Invert ();
	
	// Matrix determinant
	X Det () const;
	
	// Matrix component addition
	Mat3x3 <X> operator+ (const Mat3x3 <X>& m) const;
	
	// Matrix x Matrix multiplication
	Mat3x3 <X> operator* (const Mat3x3 <X>& m) const;
	
	Mat3x3 <X>& operator+= (const Mat3x3 <X>& m);
	
	Mat3x3 <X>& operator-= (const Mat3x3 <X>& m);
	
	Mat3x3 <X>& operator*= (const Mat3x3 <X>& m);
	
	Mat3x3 <X> operator* (double f) const;
	
	Mat3x3 <X> operator/ (double f) const;
	
	// Component scaling
	Mat3x3 <X>& operator*= (double f);
	
	// Component inverse scaling
	Mat3x3 <X>& operator/ (double f);
	
	X Trace () const;
	
	bool DiagonalDominant () const;
	
	bool StrictDiagonalDominant () const;
	
	Vec3 <X> u;
	Vec3 <X> v;
	Vec3 <X> w;
};

template <typename X>
Mat3x3 <X> operator* (double f, Mat3x3 <X> m);

template <typename X>
Mat3x3 <X> operator/ (double f, Mat3x3 <X> m);

template <typename X>
Mat3x3 <X> IdentityMat3x3 ();

template <typename X>
Mat3x3 <X> ZeroMat3x3 ();

template <typename X>
Mat3x3 <X> FilledMat3x3 (X e);

// Create a rotation matrix. It is an orthonormal linear transformation, so the resulting 3 matrix
// columns (and rows) are pairwise orthogonal and unit size. (0,0,0) will give the identity matrix.
template <typename X>
Mat3x3 <X> RotationMat3x3 (double yaw, double pitch, double roll);

template <typename X>
Mat3x3 <X> RotationMat3x3noRoll (double yaw, double pitch);

// 4x4 Matrix

template <typename X>
struct Mat4x4 {
	Mat4x4 <X>& Transpose ();
	Mat4x4 <X> operator+ (const Mat4x4 <X>& m) const;
	Mat4x4 <X> operator* (const Mat4x4 <X>& m) const;
	Mat4x4 <X>& operator+= (const Mat4x4 <X>& m);
	Mat4x4 <X>& operator-= (const Mat4x4 <X>& m);
	Mat4x4 <X>& operator*= (const Mat4x4 <X>& m);
	Mat4x4 <X> operator* (double f) const;
	Mat4x4 <X> operator/ (double f) const;
	Mat4x4 <X>& operator*= (double f);
	Mat4x4 <X>& operator/ (double f);
	X Trace () const;
	bool DiagonalDominant () const;
	bool StrictDiagonalDominant () const;
	
	Vec4 <X> u;
	Vec4 <X> v;
	Vec4 <X> w;
	Vec4 <X> x;
};

// Quaternion 
template <typename X>
struct Quaternion {
	double Length () const;
	double LengthSquared () const;
	Quaternion <X>& Unit1 ();
	Quaternion <X>& Unit0 ();
	X Dot (const Quaternion <X>& r) const;
	Quaternion <X> operator* (const Quaternion <X>& r) const;
	Quaternion <X> operator* (double f) const;
	Quaternion <X> operator/ (double f) const;
	Quaternion <X> operator+ (const Quaternion <X>& r) const;
	Quaternion <X> operator- (const Quaternion <X>& r) const;
	Quaternion <X>& operator+= (const Quaternion <X>& r);
	Quaternion <X>& operator-= (const Quaternion <X>& r);
	Quaternion <X>& operator*= (const Quaternion <X>& r);
	
	// Storage in scalar + vector format.
	X q;
	Vec3 <X> v;
};

// This creates a unit quaternion which is a rotation quaternion, so you can use these to modify
// orientation quaternions. Combined with the Quaternion to radiant conversion, you can make 6DOF
// (six degrees of freedom) movement. Alternatively, you may create the adjoint view matrix from
// your view quaternion entirely!
template <typename X>
Quaternion <X> RadQuaternion (double yaw_rad, double pitch_rad, double roll_rad);

template <typename X>
Quaternion <X> operator* (double f, const Quaternion <X>& q);

template <typename X>
Quaternion <X> IdentityQuat ();

template <typename X>
Quaternion <X> ZeroQuat ();

template <typename X>
Quaternion <X> AxisRotQuat (Vec3 <X> axis, double rad);


// Vec3 implementation
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

// Quaternion implementation
template <typename X>
double Quaternion <X>::Length () const {
	return sqrt(q * q + v.x * v.x + v.y * v.y + v.z * v.z);
}

template <typename X>
double Quaternion <X>::LengthSquared () const {
	return q * q + v.x * v.x + v.y * v.y + v.z * v.z;
}

template <typename X>
Quaternion <X>& Quaternion <X>::Unit1 () {
	auto length = Length();
	
	if(0 < length) {
		q /= length;
		v /= length;
	}
	
	else {
		q = 1.0;
		v = Vec3 <X> (0.0, 0.0, 0.0);
	}
	
	return *this;
}

template <typename X>
Quaternion <X>& Quaternion <X>::Unit0 () {
	auto length = Length();
	
	if(0 < length) {
		q /= length;
		v /= length;
	}
	
	else {
		q = 0.0;
		v = Vec3 <X> (0.0, 0.0, 0.0);
	}
	
	return *this;
}

template <typename X>
X Quaternion <X>::Dot (const Quaternion <X>& r) const {
	return q * r.q + v.x * r.v.x + v.y * r.v.y + v.z * r.v.z;
}

template <typename X>
Quaternion <X> Quaternion <X>::operator* (const Quaternion <X>& r) const {
	return Quaternion <X> (q * r.q - v.Dot(r.v), q * r.v + r.q * v + v.Cross(r.v));
}

template <typename X>
Quaternion <X> Quaternion <X>::operator* (double f) const {
	return Quaternion <X> (f * q, f * v);
}

template <typename X>
Quaternion <X> Quaternion <X>::operator/ (double f) const {
	return Quaternion <X> (q / f, v / f);
}

template <typename X>
Quaternion <X> Quaternion <X>::operator+ (const Quaternion <X>& r) const {
	return Quaternion <X> (q + r.q, v + r.v);
}

template <typename X>
Quaternion <X> Quaternion <X>::operator- (const Quaternion <X>& r) const {
	return Quaternion <X> (q - r.q, v - r.v);
}

template <typename X>
Quaternion <X>& Quaternion <X>::operator+= (const Quaternion <X>& r) {
	q += r.q;
	v += r.v;
	return *this;
}

template <typename X>
Quaternion <X>& Quaternion <X>::operator-= (const Quaternion <X>& r) {
	q -= r.q;
	v -= r.v;
	return *this;
}

template <typename X>
Quaternion <X>& Quaternion <X>::operator*= (const Quaternion <X>& r) {
	*this = *this * r;
	return *this;
};

// This creates a unit quaternion which is a rotation quaternion, so you can use these to modify
// orientation quaternions. Combined with the Quaternion to radiant conversion, you can make 6DOF
// (six degrees of freedom) movement. Alternatively, you may create the adjoint view matrix from
// your view quaternion entirely!
template <typename X>
Quaternion <X> RadQuaternion (double yaw_rad, double pitch_rad, double roll_rad) {
	
	// Only half because the quaternion rotation resolves with an addition formula for sin and cos.
	yaw_rad		/= 2;
	pitch_rad	/= 2;
	roll_rad	/= 2;
	
	auto ca = cos(yaw_rad);
	auto cv = cos(pitch_rad);
	auto cr = cos(roll_rad);
	auto sa = sin(yaw_rad);
	auto sv = sin(pitch_rad);
	auto sr = sin(roll_rad);
	
	return Quaternion <X> (
		ca * cv * cr + sa * sv * sr, Vec3 <X> (
		ca * cv * sr - sa * sv * cr,
		ca * sv * cr + sa * cv * sr,
		sa * cv * cr - ca * sv * sr));
}

template <typename X>
Quaternion <X> operator* (double f, const Quaternion <X>& q) {
	return q * f;
}

template <typename X>
Quaternion <X> IdentityQuat () {
	return Quaternion <X> (1, Vec3 <X> (0, 0, 0));
}

template <typename X>
Quaternion <X> ZeroQuat () {
	return Quaternion <X> (0, Vec3 <X> (0, 0, 0));
}

template <typename X>
Quaternion <X> AxisRotQuat (Vec3 <X> axis, double rad) {
	axis = axis.Unit1();
	rad /= 2;
	return Quaternion <X> (cos(rad), sin(rad) * axis);
}

// 3x3 Matrix implementation
template <typename X>
Mat3x3 <X>& Mat3x3 <X>::Transpose () {
	using std::swap;
	swap(u.y, v.x);
	swap(u.z, w.x);
	swap(v.z, w.y);
	return *this;
}

// Form the matrix of co-determinants.
template <typename X>
Mat3x3 <X>& Mat3x3 <X>::Adjoint () {
	Vec3 <X> new_u (
		v.y * w.z - v.z * w.y,
		v.z * w.x - v.x * w.z,
		v.x * w.y - v.y * w.x);
	
	Vec3 <X> new_v (
		w.y * u.z - w.z * u.y,
		w.z * u.x - w.x * u.z,
		w.x * u.y - w.y * u.x);
		
	Vec3 <X> new_w (
		u.y * v.z - u.z * v.y,
		u.z * v.x - u.x * v.z,
		u.x * v.y - u.y * v.x);
	
	u = new_u;
	v = new_v;
	w = new_w;
	return *this;
}

template <typename X>
bool Mat3x3 <X>::IsSingular () const {
	return abs(Det()) < 1e-9;
}

// Inverts a matrix if possible.
template <typename X>
Mat3x3 <X>& Mat3x3 <X>::Invert () {
	auto det = Det();
	
	// This matrix is singular and cannot be inverted.
	if(abs(det) < 1e-9) {
		return *this;
	}
	
	Transpose();
	Mat3x3 adjoint = *this;
	adjoint.Adjoint();
	*this = adjoint / det;
	return *this;
}

// Matrix determinant
template <typename X>
X Mat3x3 <X>::Det () const {
	return (u.Cross(v)).Dot(w);
}

// Matrix component addition
template <typename X>
Mat3x3 <X> Mat3x3 <X>::operator+ (const Mat3x3 <X>& m) const {
	return Mat3x3 <X> (u + m.u, v + m.v, w + m.w);
}

// Matrix x Matrix multiplication
template <typename X>
Mat3x3 <X> Mat3x3 <X>::operator* (const Mat3x3 <X>& m) const {
	return Mat3x3 <X> (
	Vec3 <X> (
		u.x * m.u.x + v.x * m.u.y + w.x * m.u.z,
		u.y * m.u.x + v.y * m.u.y + w.y * m.u.z,
		u.z * m.u.x + v.z * m.u.y + w.z * m.u.z),
	Vec3 <X> (
		u.x * m.v.x + v.x * m.v.y + w.x * m.v.z,
		u.y * m.v.x + v.y * m.v.y + w.y * m.v.z,
		u.z * m.v.x + v.z * m.v.y + w.z * m.v.z),
	Vec3 <X> (
		u.x * m.w.x + v.x * m.w.y + w.x * m.w.z,
		u.y * m.w.x + v.y * m.w.y + w.y * m.w.z,
		u.z * m.w.x + v.z * m.w.y + w.z * m.w.z));
}

template <typename X>
Mat3x3 <X>& Mat3x3 <X>::operator+= (const Mat3x3 <X>& m) {
	u += m.u;
	v += m.v;
	w += m.w;
	return *this;
}

template <typename X>
Mat3x3 <X>& Mat3x3 <X>::operator-= (const Mat3x3 <X>& m) {
	u -= m.u;
	v -= m.v;
	w -= m.w;
	return *this;
}

template <typename X>
Mat3x3 <X>& Mat3x3 <X>::operator*= (const Mat3x3 <X>& m) {
	*this = (*this) * m;
	return *this;
}

template <typename X>
Mat3x3 <X> Mat3x3 <X>::operator* (double f) const {
	return Mat3x3 <X> (f * u, f * v, f * w);
}

template <typename X>
Mat3x3 <X> Mat3x3 <X>::operator/ (double f) const {
	return Mat3x3 <X> (u / f, v / f, w / f);
}

// Component scaling
template <typename X>
Mat3x3 <X>& Mat3x3 <X>::operator*= (double f) {
	u *= f;
	v *= f;
	w *= f;
	return *this;
}

// Component inverse scaling
template <typename X>
Mat3x3 <X>& Mat3x3 <X>::operator/ (double f) {
	u /= f;
	v /= f;
	w /= f;
	return *this;
}

template <typename X>
X Mat3x3 <X>::Trace () const {
	return u.x + v.y + w.z;
}

template <typename X>
bool Mat3x3 <X>::DiagonalDominant () const {
	return
		abs(v.x) + abs(w.x) <= abs(u.x) &&
		abs(v.y) + abs(w.y) <= abs(v.y) &&
		abs(v.z) + abs(w.z) <= abs(w.z);
}

template <typename X>
bool Mat3x3 <X>::StrictDiagonalDominant () const {
	return
		abs(v.x) + abs(w.x) < abs(u.x) &&
		abs(v.y) + abs(w.y) < abs(v.y) &&
		abs(v.z) + abs(w.z) < abs(w.z);
}

template <typename X>
Mat3x3 <X> operator* (double f, Mat3x3 <X> m) {
	return Mat3x3 <X> (f * m.u, f * m.v, f * m.w);
}

template <typename X>
Mat3x3 <X> operator/ (double f, Mat3x3 <X> m) {
	return Mat3x3 <X> (m.u / f, m.v / f, m.w / f);
}

template <typename X>
Mat3x3 <X> IdentityMat3x3 () {
	return Mat3x3 <X> (Vec3 <X> (1, 0, 0), Vec3 <X> (0, 1, 0), Vec3 <X> (0, 0, 1));
}

template <typename X>
Mat3x3 <X> ZeroMat3x3 () {
	return Mat3x3 <X> (Vec3 <X> (0, 0, 0), Vec3 <X> (0, 0, 0), Vec3 <X> (0, 0, 0));
}

template <typename X>
Mat3x3 <X> FilledMat3x3 (X e) {
	return Mat3x3 <X> (Vec3 <X> (e, e, e), Vec3 <X> (e, e, e), Vec3 <X> (e, e, e));
}

// Create a rotation matrix. It is an orthonormal linear transformation, so the resulting 3 matrix
// columns (and rows) are pairwise orthogonal and unit size. (0,0,0) will give the identity matrix.
template <typename X>
Mat3x3 <X> RotationMat3x3 (double yaw, double pitch, double roll) {
	CosSinVals m(yaw, pitch, roll);
	Vec3 <X> forw(m.cv * m.ca, m.cv * m.sa, m.sv);
	Vec3 <X> side(-m.sa, m.ca, 0);
	Vec3 <X> up(-m.sv * m.ca, -m.sv * m.sa, m.cv);
	return Mat3x3 <X> (forw, m.cr * side - m.sr * up, m.cr * up + m.sr * side);
}

template <typename X>
Mat3x3 <X> RotationMat3x3noRoll (double yaw, double pitch) {
	CosSinVals m(yaw, pitch);
	Vec3 <X> forw(m.cv * m.ca, m.cv * m.sa, m.sv);
	Vec3 <X> side(-m.sa, m.ca, 0);
	Vec3 <X> up(-m.sv * m.ca, -m.sv * m.sa, m.cv);
	return Mat3x3 <X> (forw, side, up);
}

// 4x4 Matrix implementation
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

// 4D vector implementation
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

// Further space stuff implementations
template <typename X>
void Vec3toString (std::string& s, const Vec3 <X>& v);

template <typename X>
void Vec2toString (std::string& s, const Vec2 <X>& v);

template <typename X>
void QuatToString (std::string& s, const Quaternion <X>& q);

template <typename X>
void Mat3x3toString (std::string& s, const Mat3x3 <X>& m);

CosSinVals::CosSinVals (double yaw, double pitch, double roll) {
	ca = cos(yaw);
	cv = cos(pitch);
	cr = cos(roll);
	sa = sin(yaw);
	sv = sin(pitch);
	sr = sin(roll);
}

CosSinVals::CosSinVals (double yaw, double pitch) {
	ca = cos(yaw);
	cv = cos(pitch);
	cr = 1;
	sa = sin(yaw);
	sv = sin(pitch);
	sr = 0;
}

CosSinVals::CosSinVals () {
	ca = 1;
	cv = 1;
	cr = 1;
	sa = 0;
	sv = 0;
	sr = 0;
}

double DegToRad (double deg) {
	return SpaceConst::DegToRadFactor() * deg;
}

double RadToDeg (double rad) {
	return SpaceConst::RadToDegFactor() * rad;
}

template <typename X>
void VectorToRad::FromQuat (Quaternion <X> r) {
	yaw_rad		= 0;
	pitch_rad	= 0;
	roll_rad	= 0;
	
	auto q = r.q;
	auto v = r.v;
	auto sinvang_2 = v.x * v.z - q * v.y;
	
	if(sin_extreme < sinvang_2) {
		pitch_rad	= -SpaceConst::PiHalf();
		roll_rad	= 0;
		yaw_rad		= atan2(
			q * q - v.x * v.x + v.y * v.y - v.z * v.z,
			2 * (v.x * v.y - q * v.z)) - SpaceConst::PiHalf();
	}
	
	else if(sinvang_2 < -sin_extreme) {
		roll_rad	= 0;
		pitch_rad	= SpaceConst::PiHalf();
		yaw_rad		= atan2(
			q * q - v.x * v.x + v.y * v.y - v.z * v.z,
			2 * (v.x * v.y - q * v.z)) - SpaceConst::PiHalf();
	}
	
	else {
		auto q2_m_y2 = q * q - v.y * v.y;
		auto z2_m_x2 = v.z * v.z - v.x * v.y;
		
		pitch_rad	= asin(2 * sinvang_2);
		roll_rad	= -atan2(
			2 * (q * v.x + v.y * v.z),
			1 - 2 * (v.x * v.x + v.y * v.y));
		yaw_rad		= atan2(
			2 * (q * v.z + v.x * v.y),
			1 - 2 * (v.y * v.y + v.z * v.z));
	}
}


#include <iostream>
#include <string>
#include <sstream>

template <typename X>
void Vec3toString (std::string& s, const Vec3 <X>& v) {
	std::stringstream ss("");
	ss << "[" << v.x << ", " << v.y << ", " << v.z << "]";
	
	// Output.
	s = ss.str();
}

template <typename X>
void Vec2toString (std::string& s, const Vec2 <X>& v) {
	std::stringstream ss("");
	ss << "[" << v.x << ", " << v.y << "]";
	
	// Output.
	s = ss.str();
}

template <typename X>
void QuatToString (std::string& s, const Quaternion <X>& q) {
	std::stringstream ss("");
	ss << "[" << q.q << ", " << q.v.x << ", " << q.v.y << ", " << q.v.z << "]";
	
	// Output.
	s = ss.str();
}

template <typename X>
void Mat3x3toString (std::string& s, const Mat3x3 <X>& m) {
	std::stringstream ss("");
	ss << "[ ";
	
	auto insert_vec3 = [&ss] (const Vec3 <X>& v) {
		ss << "[";
		ss << v.x << ", " << v.y << ", " << v.z;
		ss << "] ";
	};
	
	insert_vec3(m.u);
	insert_vec3(m.v);
	insert_vec3(m.w);
	
	ss << "]";
	
	// Output.
	s = ss.str();
}

#endif
