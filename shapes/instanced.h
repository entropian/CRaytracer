#pragma once

#include "../util/mat.h"
#include "shapes.h"

typedef struct InstancedShape
{
    ObjectType obj_type;
    void* obj_ptr;
    mat4 inv_transform;
} InstancedShape;

float rayIntersectInstanced(ShadeRec* sr, InstancedShape* is, const Ray src_ray)
{
    // NOTE: ignore recursive instanced object for now
    Ray ray;    
    transformRay(&ray, is->inv_transform, src_ray);
    float t = TMAX;
    switch(is->obj_type)
    {
        case SPHERE:
            t = rayIntersectSphere(sr, (Sphere*)is->obj_ptr, ray);            
            break;
        case PLANE:
            t = rayIntersectPlane(sr, (Plane*)is->obj_ptr, ray);                        
            break;
        case RECTANGLE:
            t = rayIntersectRect(sr, (Rectangle*)is->obj_ptr, ray);
            break;
        case AABOX:
            t = rayIntersectAABox(sr, (AABox*)is->obj_ptr, ray);
            break;
        case TRIANGLE:
            t = rayIntersectTriangle(sr, (Triangle*)is->obj_ptr, ray);
            break;
        case OPENCYLINDER:
            t = rayIntersectOpenCylinder(sr, (OpenCylinder*)is->obj_ptr, ray);            
            break;
        case DISK:
            t = rayIntersectDisk(sr, (Disk*)is->obj_ptr, ray);
            break;
        case TORUS:
            t = rayIntersectTorus(sr, (Torus*)is->obj_ptr, ray);
            break;
        case SOLIDCYLINDER:
            t = rayIntersectSolidCylinder(sr, (SolidCylinder*)is->obj_ptr, ray);
            break;
    }
    if(t < TMAX)
    {
        mat4 normal_mat;
        mat4_transpose(normal_mat, is->inv_transform);
        vec4 normal, transformed_normal;
        vec4FromVec3(normal, sr->normal, 0.0f);
        mat4_mult_vec4(transformed_normal, normal_mat, normal);
        vec3FromVec4(sr->normal, transformed_normal);
        getPointOnRay(sr->hit_point, src_ray, t);
    }
    return t;
}

float shadowRayIntersectInstanced(InstancedShape* is, const Ray src_ray)
{
    // NOTE: ignore recursive instanced object for now
    Ray shadow_ray;
    transformRay(&shadow_ray, is->inv_transform, src_ray);    
    float t = TMAX;
    switch(is->obj_type)
    {
        case SPHERE:
            t = shadowRayIntersectSphere((Sphere*)is->obj_ptr, shadow_ray);            
            break;
        case PLANE:
            t = shadowRayIntersectPlane((Plane*)is->obj_ptr, shadow_ray);                        
            break;
        case RECTANGLE:
            t = shadowRayIntersectRect((Rectangle*)is->obj_ptr, shadow_ray);
            break;
        case AABOX:
            t = shadowRayIntersectAABox((AABox*)is->obj_ptr, shadow_ray);
            break;
        case TRIANGLE:
            t = shadowRayIntersectTriangle((Triangle*)is->obj_ptr, shadow_ray);
            break;
        case OPENCYLINDER:
            t = shadowRayIntersectOpenCylinder((OpenCylinder*)is->obj_ptr, shadow_ray);            
            break;
        case DISK:
            t = shadowRayIntersectDisk((Disk*)is->obj_ptr, shadow_ray);
            break;
        case TORUS:
            t = shadowRayIntersectTorus((Torus*)is->obj_ptr, shadow_ray);
            break;
        case SOLIDCYLINDER:
            t = shadowRayIntersectSolidCylinder((SolidCylinder*)is->obj_ptr, shadow_ray);
            break;
    }
    return t;
}
