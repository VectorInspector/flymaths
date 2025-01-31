#ifndef SPACE_H
#define SPACE_H

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

#include "quaternion.h"
#include "mat3x3.h"
#include "vec3.h"
#include "vec2.h"
#include <string>

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
