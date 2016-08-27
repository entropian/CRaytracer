#pragma once

#include <emmintrin.h>
#include "../util/vec.h"
#include "../util/ray.h"
#include "../util/shaderec.h"
#include "../util/constants.h"
#include "../materials.h"
#include "../mesh.h"
#include "../texture.h"
#include "../util/mat.h"
#include "../util/math.h"


typedef struct Triangle
{
    bool shadow;
    vec3 v0, v1, v2;
    vec3 normal;
    Material* mat;
}Triangle;

typedef struct FlatTriangle_s
{
    bool shadow;
    int i0, i1, i2;
    vec3 v0, v1, v2;
    vec3 normal;
    Mesh* mesh_ptr;
    Material* mat;
}FlatTriangle;

typedef struct SmoothTriangle_s
{
    bool shadow;
    int i0, i1, i2;
    vec3 v0, v1, v2;
    mat3* normal_mat;
    Mesh* mesh_ptr;
    Material* mat;
}SmoothTriangle;

void calcTriangleNormal(vec3 r, const vec3 v0, const vec3 v1, const vec3 v2);

float calcTriangleIntersect(float* beta_out, float* gamma_out,
                            const vec3 v0, const vec3 v1, const vec3 v2, const Ray ray);

// NOTE: perhaps rewrite this so that it makes more sense to me
float rayIntersectTriangle(ShadeRec* sr, Triangle* tri, const Ray ray);

void interpTexcoord(vec2 uv_out,const float beta, const float gamma,
                    const Mesh* mesh_ptr, const int i0, const int i1, const int i2);

float rayIntersectFlatTriangle(ShadeRec* sr, FlatTriangle* tri, const Ray ray);

void interpTriangleVec3(vec3 out, const float beta, const float gamma,
                        const vec3 v0, const vec3 v1, const vec3 v2);

float rayIntersectSmoothTriangle(ShadeRec* sr, SmoothTriangle* tri, const Ray ray);

inline float shadowRayIntersectTriangle(const Triangle* tri, const Ray ray)
{
    float beta, gamma; // Unused
    return calcTriangleIntersect(&beta, &gamma, tri->v0, tri->v1, tri->v2, ray);
}

inline float shadowRayIntersectFlatTriangle(const FlatTriangle* tri, const Ray ray)
{
    float beta, gamma;    // Unused
    return calcTriangleIntersect(&beta, &gamma, tri->v0, tri->v1, tri->v2, ray);
}

inline float shadowRayIntersectSmoothTriangle(const SmoothTriangle* tri, const Ray ray)
{
    float beta, gamma;  // Unused
    return calcTriangleIntersect(&beta, &gamma, tri->v0, tri->v1, tri->v2, ray);
}


