#pragma once

#include "../util/vec.h"
#include "../util/ray.h"
#include "../util/shaderec.h"
#include "../util/constants.h"
#include "../materials.h"
#include "../objloader/objloader.h"

typedef struct Triangle
{
    bool shadow;
    vec3 v0, v1, v2;
    vec3 normal;
    Material* mat;
}Triangle;

typedef struct 
{
    bool shadow;
    int i0, i1, i2;
    vec3 normal;
    OBJShape* mesh_ptr;
    Material* mat;
}MeshTriangle;

void getMeshTriangleVertPos(vec3 v0, vec3 v1, vec3 v2, const MeshTriangle* tri)
{
    OBJShape* mesh = tri->mesh_ptr;
    int index = tri->i0 * 3;
    vec3_assign(v0, mesh->positions[index], mesh->positions[index+1], mesh->positions[index+2]);
    index = tri->i1 * 3;    
    vec3_assign(v1, mesh->positions[index], mesh->positions[index+1], mesh->positions[index+2]);
    index = tri->i2 * 3;
    vec3_assign(v2, mesh->positions[index], mesh->positions[index+1], mesh->positions[index+2]);
}

//void calcTriangleNormal(Triangle* triangle)
void calcTriangleNormal(vec3 r, const vec3 v0, const vec3 v1, const vec3 v2)
{
    // Assume the vertices are ordered counterclock wise
    vec3 tmp1, tmp2, tmp3;
    vec3_sub(tmp1, v1, v0);
    vec3_sub(tmp2, v2, v0);
    vec3_cross(tmp3, tmp1, tmp2);
    vec3_normalize(r, tmp3);    
}

//float calcTriangleIntersect(const Triangle* tri, const Ray ray)
float calcTriangleIntersect(const vec3 v0, const vec3 v1, const vec3 v2, const Ray ray)
{
    // o + td = v0 + v(v1-v0) + w(v2-v0)
    // v(v0-v1) + w(v0-v2) + td = v0 - o
    /*
    float a = tri->v0[0] - tri->v1[0], b = tri->v0[0] - tri->v2[0], c = ray.direction[0], d = tri->v0[0] - ray.origin[0];
    float e = tri->v0[1] - tri->v1[1], f = tri->v0[1] - tri->v2[1], g = ray.direction[1], h = tri->v0[1] - ray.origin[1];
    float i = tri->v0[2] - tri->v1[2], j = tri->v0[2] - tri->v2[2], k = ray.direction[2], l = tri->v0[2] - ray.origin[2];
    */
    float a = v0[0] - v1[0], b = v0[0] - v2[0], c = ray.direction[0], d = v0[0] - ray.origin[0];
    float e = v0[1] - v1[1], f = v0[1] - v2[1], g = ray.direction[1], h = v0[1] - ray.origin[1];
    float i = v0[2] - v1[2], j = v0[2] - v2[2], k = ray.direction[2], l = v0[2] - ray.origin[2];

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
    float t = calcTriangleIntersect(tri->v0, tri->v1, tri->v2, ray);

    vec3_copy(sr->normal, tri->normal);
    getPointOnRay(sr->hit_point, ray, t);
    vec3_negate(sr->wo, ray.direction);    
    sr->mat = tri->mat;
    return t;
}

float rayIntersectMeshTriangle(ShadeRec* sr, MeshTriangle* tri, const Ray ray)
{
    vec3 v0, v1, v2;
    getMeshTriangleVertPos(v0, v1, v2, tri);    
    float t = calcTriangleIntersect(v0, v1, v2, ray);

    vec3_copy(sr->normal, tri->normal);
    getPointOnRay(sr->hit_point, ray, t);
    vec3_negate(sr->wo, ray.direction);    
    sr->mat = tri->mat;
    return t;
}

float shadowRayIntersectTriangle(const Triangle* tri, const Ray ray)
{
    return calcTriangleIntersect(tri->v0, tri->v1, tri->v2, ray);
}

float shadowRayIntersectMeshTriangle(const MeshTriangle* tri, const Ray ray)
{
    /*
    OBJShape* mesh = tri->mesh_ptr;
    vec3 v0, v1, v2;    
    int index = tri->i0 * 3;
    vec3_assign(v0, mesh->positions[index], mesh->positions[index+1], mesh->positions[index+2]);
    index = tri->i1 * 3;    
    vec3_assign(v1, mesh->positions[index], mesh->positions[index+1], mesh->positions[index+2]);
    index = tri->i2 * 3;
    vec3_assign(v2, mesh->positions[index], mesh->positions[index+1], mesh->positions[index+2]);
    */
    vec3 v0, v1, v2;
    getMeshTriangleVertPos(v0, v1, v2, tri);
    return calcTriangleIntersect(v0, v1, v2, ray);
}
