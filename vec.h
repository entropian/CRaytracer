#ifndef VEC_H
#define VEC_H

#include <math.h>
#include <assert.h>
#include <algorithm>

// TODO: figure out this stuff
static const double PI = 3.14159265358979323846264338327950288;
static const double EPS = 1e-8;
static const double EPS2 = EPS*EPS;
static const double EPS3 = EPS*EPS*EPS;

template <int n>
class Vec
{
	float f_[n];

public:
	Vec()
	{
		for (int i = 0; i < n; i++)
			f_[i] = 0;
	}

	Vec(float t)
	{
		for (int i = 0; i < n; i++)
			f_[i] = t;
	}


	Vec(float a, float b)
	{
		assert(n == 2);
		f_[0] = a;
		f_[1] = b;
	}

	Vec(float a, float b, float c)
	{
		assert(n == 3);
		f_[0] = a;
		f_[1] = b;
		f_[2] = c;
	}

	Vec(float a, float b, float c, float d)
	{
		assert(n == 4);
		f_[0] = a;
		f_[1] = b;
		f_[2] = c;
		f_[3] = d;
	}

	// truncate if m < n, or extend the vector with extendValue
	template<int m>
	explicit Vec(const Vec<m>& a, const float extendValue = 0.0)
	{
		for (int i = 0; i < std::min(m, n); i++)
			f_[i] = a[i];

		for (int i = std::min(m, n); i < n; i++)
			f_[i] = extendValue;
	}

	float& operator [](const int i)
	{
		return f_[i];
	}

	const float& operator[](const int i) const
	{
		return f_[i];
	}

	Vec operator - () const
	{
		return Vec(*this) *= -1;
	}

	Vec& operator +=(const Vec& v)
	{
		for (int i = 0; i < n; i++)
			f_[i] += v[i];
		return *this;
	}

	Vec& operator -=(const Vec& v)
	{
		for (int i = 0; i < n; i++)
			f_[i] -= v[i];
		return *this;
	}

	Vec& operator *=(const float a)
	{
		for (int i = 0; i < n; i++)
			f_[i] *= a;
		return *this;
	}

	Vec& operator /=(const float a)
	{
		const float inva = 1 / a;
		for (int i = 0; i < n; i++)
			f_[i] *= inva;
		return *this;
	}

	Vec operator +(const Vec& v) const
	{
		return Vec(*this) += v;
	}

	Vec operator -(const Vec& v) const
	{
		return Vec(*this) -= v;
	}

	Vec operator *(const float a) const
	{
		return Vec(*this) *= a;
	}

	Vec operator /(const float a) const
	{
		return Vec(*this) /= a;
	}

	Vec& normalize()
	{
		assert((double)dot(*this, *this) > EPS2);
		return *this /= sqrt(dot(*this, *this));
	}
};

typedef Vec<2> Vec2;
typedef Vec<3> Vec3;
typedef Vec<4> Vec4;

inline Vec3 cross(const Vec3& a, const Vec3& b)
{
	return Vec3((a[1] * b[2]) - (a[2] * b[1]), (a[2] * b[0]) - (a[0] * b[2]), (a[0] * b[1]) - (a[1] * b[0]));
}

template<int n>
inline float dot(const Vec<n>& a, const Vec<n>& b)
{
	float r = 0;
	for (int i = 0; i < n; i++)
		r += a[i] * b[i];
	return r;
}

template<int n>
inline float norm2(const Vec<n>& v)
{
	return dot(v, v);
}

template<int n>
inline float norm(const Vec<n>& v)
{
	return sqrt(dot(v, v));
}

// Return a normalized vector without modifying the input
template<int n>
inline Vec<n> normalize(const Vec<n>& v)
{
    /*
    assert((double)(dot(v, v)) > EPS2);
    return v / norm(v);
    */

    if((double)(dot(v, v)) < EPS2)
        return Vec3(0.0f, 0.0f, 0.0f);
    else
        return v / norm(v);

}


#endif
