#pragma once

#include "../util/vec.h"
#include "../util/ray.h"
#include "../util/shaderec.h"
#include "../util/constants.h"
#include "../materials.h"

typedef struct Triangle
{
    bool shadow;
    vec3 v0, v1, v2;
    vec3 normal;
    Material mat;
} Triangle;

void calcTriangleNormal(Triangle* triangle)
{
    // Assume the vertices are ordered counterclock wise
    vec3 tmp1, tmp2, tmp3;
    vec3_sub(tmp1, triangle->v1, triangle->v0);
    vec3_sub(tmp2, triangle->v2, triangle->v0);
    vec3_cross(tmp3, tmp1, tmp2);
    vec3_normalize(triangle->normal, tmp3);    
}

float calcTriangleIntersect(const Triangle* tri, const Ray ray)
{
    // o + td = v0 + v(v1-v0) + w(v2-v0)
    // v(v0-v1) + w(v0-v2) + td = v0 - o

    float a = tri->v0[0] - tri->v1[0], b = tri->v0[0] - tri->v2[0], c = ray.direction[0], d = tri->v0[0] - ray.origin[0];
    float e = tri->v0[1] - tri->v1[1], f = tri->v0[1] - tri->v2[1], g = ray.direction[1], h = tri->v0[1] - ray.origin[1];
    float i = tri->v0[2] - tri->v1[2], j = tri->v0[2] - tri->v2[2], k = ray.direction[2], l = tri->v0[2] - ray.origin[2];

    float m = f*k - g*j, n = h*k - g*l, p = f*l - h*j;
    float q = g*i - e*k, s = e*j - f*i;

    float inv_denom = 1.0f / (a*m + b*q + c*s);

    float e1 = d*m - b*n - c*p;
    float beta = e1 * inv_denom;

    if(beta < 0.0f)
    {
        return TMAX;
    }

    float r = e*l - h*i;
    float e2 = a*n + d*q + c*r;
    float gamma = e2 * inv_denom;

    if(gamma < 0.0f)
    {
        return TMAX;
    }

    if(beta + gamma > 1.0f)
    {
        return TMAX;
    }

    float e3 = a*p - b*r + d*s;
    float t = e3 * inv_denom;

    if(t < K_EPSILON)
    {
        return TMAX;
    }
    return t;
}

// NOTE: perhaps rewrite this so that it makes more sense to me
float rayIntersectTriangle(ShadeRec* sr, Triangle* tri, const Ray ray)
{
    float t = calcTriangleIntersect(tri, ray);

    vec3_copy(sr->normal, tri->normal);
    vec3 displacement;
    vec3_scale(displacement, ray.direction, t);
    vec3_add(sr->hit_point, ray.origin, displacement);
    vec3_negate(sr->wo, ray.direction);    
    sr->mat = &(tri->mat);
    return t;
}

float shadowRayIntersectTriangle(const Triangle* tri, const Ray ray)
{
    return calcTriangleIntersect(tri, ray);
}
