#pragma once

#include <emmintrin.h>
#include "../util/vec.h"
#include "../util/ray.h"
#include "../util/shaderec.h"
#include "../util/constants.h"
#include "../materials.h"
#include "../mesh.h"
#include "../texture.h"

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
    vec3 n0, n1, n2;
    vec3 tan0, tan1, tan2;     // Tentative
    Mesh* mesh_ptr;
    Material* mat;
}SmoothTriangle;

void calcTriangleNormal(vec3 r, const vec3 v0, const vec3 v1, const vec3 v2)
{
    // Assume the vertices are ordered counterclock wise
    vec3 tmp1, tmp2, tmp3;
    vec3_sub(tmp1, v1, v0);
    vec3_sub(tmp2, v2, v0);
    vec3_cross(tmp3, tmp1, tmp2);
    vec3_normalize(r, tmp3);    
}

float calcTriangleIntersect(float* beta_out, float* gamma_out,
                            const vec3 v0, const vec3 v1, const vec3 v2, const Ray ray)
{
    // o + td = v0 + v(v1-v0) + w(v2-v0)
    // v(v0-v1) + w(v0-v2) + td = v0 - o

    float a = v0[0] - v1[0], b = v0[0] - v2[0], c = ray.direction[0], d = v0[0] - ray.origin[0];      
    float e = v0[1] - v1[1], f = v0[1] - v2[1], g = ray.direction[1], h = v0[1] - ray.origin[1];
    float i = v0[2] - v1[2], j = v0[2] - v2[2], k = ray.direction[2], l = v0[2] - ray.origin[2];
    
    float m = f*k - g*j, n = h*k - g*l, p = f*l - h*j;
    float q = g*i - e*k, s = e*j - f*i;

    float inv_denom = 1.0f / (a*m + b*q + c*s);

    float e1 = d*m - b*n - c*p;
    float beta = e1 * inv_denom;

    float r = e*l - h*i;
    float e2 = a*n + d*q + c*r;
    float gamma = e2 * inv_denom;
    *beta_out = beta;
    *gamma_out = gamma;

    if(beta < 0.0f)
    {
        return TMAX;
    }    

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
    float gamma, beta; // For smooth triangles. Unused here.
    float t = calcTriangleIntersect(&beta, &gamma, tri->v0, tri->v1, tri->v2, ray);

    vec3_copy(sr->normal, tri->normal);
    getPointOnRay(sr->hit_point, ray, t);
    vec3_negate(sr->wo, ray.direction);    
    sr->mat = tri->mat;
    return t;
}

void interpTexcoord(vec2 uv_out,const float beta, const float gamma,
                    const Mesh* mesh_ptr, const int i0, const int i1, const int i2)
{
    vec2 uv0 = {mesh_ptr->texcoords[i0*2], mesh_ptr->texcoords[i0*2 + 1]};
    vec2 uv1 = {mesh_ptr->texcoords[i1*2], mesh_ptr->texcoords[i1*2 + 1]};
    vec2 uv2 = {mesh_ptr->texcoords[i2*2], mesh_ptr->texcoords[i2*2 + 1]};

    vec2 tmp;
    vec2_scale(uv_out, uv0, 1.0f - beta - gamma);
    vec2_scale(tmp, uv1, beta);
    vec2_add(uv_out, uv_out, tmp);
    vec2_scale(tmp, uv2, gamma);
    vec2_add(uv_out, uv_out, tmp);
}

float rayIntersectFlatTriangle(ShadeRec* sr, FlatTriangle* tri, const Ray ray)
{
    float gamma, beta; // For smooth triangles. Unused here.    
    float t = calcTriangleIntersect(&beta, &gamma, tri->v0, tri->v1, tri->v2, ray);

    if(tri->mesh_ptr->num_texcoords > 0 && tri->mat->tex_flags != NO_TEXTURE)
    {
        interpTexcoord(sr->uv, beta, gamma, tri->mesh_ptr, tri->i0, tri->i1, tri->i2);
    }    
    vec3_copy(sr->normal, tri->normal);
    getPointOnRay(sr->hit_point, ray, t);
    vec3_negate(sr->wo, ray.direction);    
    sr->mat = tri->mat;
    return t;
}

void interpTriangleVec3(vec3 out, const float beta, const float gamma,
                        const vec3 v0, const vec3 v1, const vec3 v2)
{
    vec3 tmp;
    vec3_scale(out, v0, 1.0f - beta - gamma);
    vec3_scale(tmp, v1, beta);
    vec3_add(out, out, tmp);
    vec3_scale(tmp, v2, gamma);    
    vec3_add(out, out, tmp);
    vec3_normalize(out, out);
}

float rayIntersectSmoothTriangle(ShadeRec* sr, SmoothTriangle* tri, const Ray ray)
{
    float gamma, beta;  
    float t = calcTriangleIntersect(&beta, &gamma, tri->v0, tri->v1, tri->v2, ray);
    if(t == TMAX){return t;}

    interpTriangleVec3(sr->normal, beta, gamma, tri->n0, tri->n1, tri->n2);    
    if(tri->mesh_ptr->num_texcoords > 0 && tri->mat->tex_flags != NO_TEXTURE)
    {
        interpTexcoord(sr->uv, beta, gamma, tri->mesh_ptr, tri->i0, tri->i1, tri->i2);
        if(tri->mat->tex_flags & NORMAL)
        {
            vec3 tangent, binormal, tex_normal, normal, tmp;
            interpTriangleVec3(tangent, beta, gamma, tri->tan0, tri->tan1, tri->tan2);
            vec3_cross(binormal, sr->normal, tangent);
            vec3_normalize(binormal, binormal);
            
            getTexColor(tex_normal, tri->mat->tex_array[1], sr->uv);
            orthoNormalTransform(normal, tangent, binormal, sr->normal, tex_normal);
            vec3_normalize(normal, normal);
            float tmp_float = vec3_dot(normal, sr->normal);
            vec3_assign(tmp, tmp_float, tmp_float, tmp_float);
            vec3_copy(sr->normal, normal);
        }
    }
    getPointOnRay(sr->hit_point, ray, t);
    vec3_negate(sr->wo, ray.direction);    
    sr->mat = tri->mat;
    return t;
}

/*
float rayIntersectSmoothTriangle(ShadeRec* sr, SmoothTriangle* tri, const Ray ray)
{
    float gamma, beta;  
    float t = calcTriangleIntersect(&beta, &gamma, tri->v0, tri->v1, tri->v2, ray);
    if(t == TMAX){return t;}

    interpTriangleVec3(sr->normal, beta, gamma, tri->n0, tri->n1, tri->n2);    
    if(tri->mesh_ptr->num_texcoords > 0 && tri->mat->tex_flags != NO_TEXTURE)
    {
        interpTexcoord(sr->uv, beta, gamma, tri->mesh_ptr, tri->i0, tri->i1, tri->i2);
        if(tri->mat->tex_flags & NORMAL)
        {
            vec3 tangent, binormal, tex_normal, normal;
            interpTriangleVec3(tangent, beta, gamma, tri->tan0, tri->tan1, tri->tan2);
            vec3_cross(binormal, sr->normal, binormal);
            getTexColor(tex_normal, tri->mat->tex_array[1], sr->uv);
            orthoNormalTransform(normal, tangent, binormal, sr->normal, tex_normal);
            // Check if we're mirrored
            if(vec3_dot(normal, sr->normal) < 0.0f)
            {
                vec3_negate(binormal, binormal);
                orthoNormalTransform(normal, tangent, binormal, sr->normal, tex_normal);                
            }
            vec3_normalize(normal, normal);            
            vec3_copy(sr->normal, normal);
        }
    }
    getPointOnRay(sr->hit_point, ray, t);
    vec3_negate(sr->wo, ray.direction);    
    sr->mat = tri->mat;
    return t;
}
*/


float shadowRayIntersectTriangle(const Triangle* tri, const Ray ray)
{
    float beta, gamma; // Unused
    return calcTriangleIntersect(&beta, &gamma, tri->v0, tri->v1, tri->v2, ray);
}

float shadowRayIntersectFlatTriangle(const FlatTriangle* tri, const Ray ray)
{
    float beta, gamma;    // Unused
    return calcTriangleIntersect(&beta, &gamma, tri->v0, tri->v1, tri->v2, ray);
}

float shadowRayIntersectSmoothTriangle(const SmoothTriangle* tri, const Ray ray)
{
    float beta, gamma;  // Unused
    return calcTriangleIntersect(&beta, &gamma, tri->v0, tri->v1, tri->v2, ray);
}


