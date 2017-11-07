#pragma once

#include <stdio.h>
#include <assert.h>
#include <emmintrin.h>
#include <xmmintrin.h>
#include "vec.h"
#include "constants.h"
#include "mat.h"
#include "ray.h"


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
