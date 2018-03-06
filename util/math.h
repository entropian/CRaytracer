#pragma once

#include <stdio.h>
#include <assert.h>
#include <emmintrin.h>
#include <xmmintrin.h>
#include "vec.h"
#include "constants.h"
#include "mat.h"
#include "ray.h"
#include "util.h"

inline float cosTheta(const vec3 w) { return w[2]; }
inline float cos2Theta(const vec3 w) { return w[2] * w[2]; }
inline float absCosTheta(const vec3 w) { return fabs(w[2]); }
inline float sin2Theta(const vec3 w) {
    return max((float)0, (float)1 - cos2Theta(w));
}

inline float sinTheta(const vec3 w) { return sqrt(sin2Theta(w)); }

inline float tanTheta(const vec3 w) { return sinTheta(w) / cosTheta(w); }

inline float tan2Theta(const vec3 w) {
    return sin2Theta(w) / cos2Theta(w);
}

inline float cosPhi(const vec3 w) {
    float sin_theta = sinTheta(w);
    return (sin_theta == 0) ? 1 : clamp(w[0] / sin_theta, -1.0f, 1.0f);
}

inline float sinPhi(const vec3 w) {
    float sin_theta = sinTheta(w);
    return (sin_theta == 0) ? 0 : clamp(w[1] / sin_theta, -1.0f, 1.0f);
}

inline float cos2Phi(const vec3 w) { return cosPhi(w) * cosPhi(w); }

inline float sin2Phi(const vec3 w) { return sinPhi(w) * sinPhi(w); }

inline float cosDPhi(const vec3 wa, const vec3 wb) {
    return clamp(
        (wa[0] * wb[0] + wa[1] * wb[1]) / sqrt((wa[0] * wa[0] + wa[1] * wa[1]) *
                                                (wb[0] * wb[0] + wb[1] * wb[1])),
        -1, 1);
}

inline bool sameHemisphere(const vec3 a, const vec3 b)
{
    return a[2] * b[2] > 0.0f;
}


void orthoNormalTransform(vec3 r, const vec3 u, const vec3 v, const vec3 w, const vec3 a);
void transposeTransform(vec3 r, const vec3 u, const vec3 v, const vec3 w, const vec3 a);

void getVec3InLocalBasis(vec3 r, const vec3 a, const vec3 normal);

void transformRay(Ray* dest_ray, const mat4 mat, const Ray src_ray);

void defaultInvTransform(mat4 r, const vec3 scaling, const vec3 axis, const float theta,
                         const vec3 translation);

inline float degToRad(const float degree)
{
    return degree / 180.0f * (float)PI;
}

void eulerAngToMat3(mat3 r, const vec3 euler_ang);

void eulerAngToMat4(mat4 r, const vec3 euler_ang);

inline void printMat3(mat3 a)
{
    for(int i = 0; i < 3; i++)
    {
        printf("%f %f %f\n", a[i][0], a[i][1], a[i][2]);
    }
}

float distanceSquared(const vec3 a, const vec3 b);


inline float lerp(const float x, const float a, const float b)
{
    assert(x >= 0.0f && x <= 1.0f);
    return a + (b - a) * x;
}

inline void sphericalDirection(vec3 out, float sin_theta, float cos_theta, float phi) {
    vec3_assign(out, sin_theta * cosf(phi), sin_theta * sinf(phi), cos_theta);
}

inline void cartesianToSpherical(vec2 out, const vec3 in)
{
    out[0] = atan2f(in[2], in[0]);
    out[0] += PI;
    out[1] = acosf(in[1]);
}

inline void sphericalToUV(vec2 uv, const vec2 spherical)
{
    uv[0] = spherical[0] / (2.0f * PI);
    uv[1] = 1.0f - spherical[1] / PI;
}

__m128 fourKnotSplineSSE(__m128* x, __m128* k0, __m128* k1, __m128* k2, __m128* k3);


// Assumes knots[] has four elements
inline float fourKnotSpline(const float x, const float knots[])
{
    float c3 = -0.5f * knots[0] + 1.5f * knots[1] - 1.5f * knots[2] + 0.5f * knots[3];
    float c2 = knots[0] - 2.5f * knots[1] + 2.0f * knots[2] - 0.5f * knots[3];
    float c1 = 0.5f * (-knots[0] + knots[2]);
    float c0 = knots[1];
    return ((c3*x + c2)*x + c1)*x + c0;
}

/*
  Utility functions to find cubic and quartic roots.
  Copied from code for Ray Tracing from the Ground up.
  Author: Jochen Schwarze (schwarze@isa.de)
 */

int solveQuadric(double c[3], double s[2]);
int solveCubic(double c[4], double s[3]);
int solveQuartic(double c[5], double s[4]);
