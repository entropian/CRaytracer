#pragma once

#include "../util/mat.h"
#include "objecttype.h"
#include "../aabb.h"

typedef struct CompoundObject
{
    bool shadow;
    int num_obj;
    ObjectType* obj_types;
    void** obj_ptrs;
    AABB aabb;    
} CompoundObject;

float rayIntersectCompound(ShadeRec* sr, CompoundObject* co, const Ray ray);
float shadowRayIntersectCompound(CompoundObject* co, const Ray shadow_ray);

typedef struct InstancedShape
{
    ObjectType obj_type;
    void* obj_ptr;
    Material* mat;
    mat4 inv_transform;
} InstancedShape;

void freeInstancedShape(InstancedShape* is);

void freeCompoundObject(CompoundObject* co)
{
    for(int i = 0; i < co->num_obj; i++)
    {
        if(co->obj_types[i] == INSTANCED)
        {
            freeInstancedShape((InstancedShape*)(co->obj_ptrs[i]));
        }else if(co->obj_types[i] == COMPOUND)
        {
            freeCompoundObject((CompoundObject*)(co->obj_ptrs[i]));
        }else
        {
            free(co->obj_ptrs[i]);
        }
    }
    free(co->obj_ptrs);
    free(co->obj_types);
}

void freeInstancedShape(InstancedShape* is)
{
    if(is->obj_type == COMPOUND)
    {
        freeCompoundObject((CompoundObject*)is->obj_ptr);
    }else
    {
        free(is->obj_ptr);
    }
}

float rayIntersectInstanced(ShadeRec* sr, InstancedShape* is, const Ray src_ray)
{
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
    case INSTANCED:
        t = rayIntersectInstanced(sr, (InstancedShape*)is->obj_ptr, ray);
        break;
    case COMPOUND:
        t = rayIntersectCompound(sr, (CompoundObject*)is->obj_ptr, ray);
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
        is->mat = sr->mat;
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
        case INSTANCED:
            t = shadowRayIntersectInstanced((InstancedShape*)is->obj_ptr, shadow_ray);
            break;
        case COMPOUND:
            t = shadowRayIntersectCompound((CompoundObject*)is->obj_ptr, shadow_ray);
            break;                        
    }
    return t;
}

float rayIntersectCompound(ShadeRec* sr, CompoundObject* co, const Ray ray)
{
    if(rayIntersectAABB(&(co->aabb), ray) == TMAX)
    {
        return TMAX;
    }       
    float min_t = TMAX, tmp_t = TMAX;
    ShadeRec min_sr, tmp_sr;
    int num_obj = co->num_obj;
    for(int i = 0; i < num_obj; i++)
    {
        void* obj_ptr = co->obj_ptrs[i];
        switch(co->obj_types[i])
        {
        case SPHERE:
            tmp_t = rayIntersectSphere(&tmp_sr, (Sphere*)obj_ptr, ray);            
            break;
        case PLANE:
            tmp_t = rayIntersectPlane(&tmp_sr, (Plane*)obj_ptr, ray);                        
            break;
        case RECTANGLE:
            tmp_t = rayIntersectRect(&tmp_sr, (Rectangle*)obj_ptr, ray);
            break;
        case AABOX:
            tmp_t = rayIntersectAABox(&tmp_sr, (AABox*)obj_ptr, ray);
            break;
        case TRIANGLE:
            tmp_t = rayIntersectTriangle(&tmp_sr, (Triangle*)obj_ptr, ray);
            break;
        case OPENCYLINDER:
            tmp_t = rayIntersectOpenCylinder(&tmp_sr, (OpenCylinder*)obj_ptr, ray);            
            break;
        case DISK:
            tmp_t = rayIntersectDisk(&tmp_sr, (Disk*)obj_ptr, ray);
            break;
        case TORUS:
            tmp_t = rayIntersectTorus(&tmp_sr, (Torus*)obj_ptr, ray);
            break;
        case INSTANCED:
            tmp_t = rayIntersectInstanced(&tmp_sr, (InstancedShape*)obj_ptr, ray);
            break;
        }
        if(tmp_t < min_t)
        {
            min_t = tmp_t;
            min_sr = tmp_sr;
        }
    }
    *sr = min_sr;
    return min_t;
}

float shadowRayIntersectCompound(CompoundObject* co, const Ray shadow_ray)
{
    if(rayIntersectAABB(&(co->aabb), shadow_ray) == TMAX)
    {
        return TMAX;
    }               
    float t = TMAX;
    int num_obj = co->num_obj;
    for(int i = 0; i < num_obj; i++)
    {
        void* obj_ptr = co->obj_ptrs[i];
        switch(co->obj_types[i])
        {
        case SPHERE:
            t = shadowRayIntersectSphere((Sphere*)obj_ptr, shadow_ray);  
            break;
        case PLANE:
            t = shadowRayIntersectPlane((Plane*)obj_ptr, shadow_ray);                        
            break;
        case RECTANGLE:
            t = shadowRayIntersectRect((Rectangle*)obj_ptr, shadow_ray);
            break;
        case AABOX:
            t = shadowRayIntersectAABox((AABox*)obj_ptr, shadow_ray);
            break;
        case TRIANGLE:
            t = shadowRayIntersectTriangle((Triangle*)obj_ptr, shadow_ray);
            break;
        case OPENCYLINDER:
            t = shadowRayIntersectOpenCylinder((OpenCylinder*)obj_ptr, shadow_ray);            
            break;
        case DISK:
            t = shadowRayIntersectDisk((Disk*)obj_ptr, shadow_ray);
            break;
        case TORUS:
            t = shadowRayIntersectTorus((Torus*)obj_ptr, shadow_ray);
            break;
        case INSTANCED:
            t = shadowRayIntersectInstanced((InstancedShape*)obj_ptr, shadow_ray);
            break;
        }
        if(t < TMAX)
        {
            return t;
        }
    }
    return t;
}

typedef CompoundObject SolidCylinder;

void calcAABBSolidCylinder(AABB* aabb, const Disk* top, const Disk* bottom, const OpenCylinder* tube)
{
    aabb->min[0] = -(tube->radius);
    aabb->min[1] = -(tube->half_height);
    aabb->min[2] = -(tube->radius);

    aabb->max[0] = tube->radius;
    aabb->max[1] = tube->half_height;
    aabb->max[2] = (tube->radius);    
}

SolidCylinder* initSolidCylinder(const float radius, const float half_height, const float phi,
                    Material* mat, const bool shadow)
{
    int num_obj = 3;
    SolidCylinder* sc = (SolidCylinder*)malloc(sizeof(SolidCylinder));
    sc->obj_ptrs = (void**)malloc(sizeof(void*) * num_obj);
    sc->obj_types = (ObjectType*)malloc(sizeof(ObjectType) * num_obj);
    sc->shadow = shadow;
    Disk* top = (Disk*)malloc(sizeof(Disk));
    Disk* bottom = (Disk*)malloc(sizeof(Disk));
    OpenCylinder* tube = (OpenCylinder*)malloc(sizeof(OpenCylinder));
    top->radius = bottom->radius = tube->radius = radius;
    tube->half_height = half_height;
    vec3_assign(top->center, 0.0f, half_height, 0.0f);
    vec3_assign(bottom->center, 0.0f, -half_height, 0.0f);
    vec3_assign(top->normal, 0.0f, 1.0f, 0.0f);
    vec3_assign(bottom->normal, 0.0f, -1.0f, 0.0f);
    tube->phi = phi;
    tube->normal_type = OPEN;
    top->mat = bottom->mat = tube->mat = mat;

    sc->obj_ptrs[0] = top;
    sc->obj_types[0] = DISK;
    sc->obj_ptrs[1] = bottom;
    sc->obj_types[1] = DISK;
    sc->obj_ptrs[2] = tube;
    sc->obj_types[2] = OPENCYLINDER;
    sc->num_obj = num_obj;
    calcAABBSolidCylinder(&(sc->aabb), top, bottom, tube);
    return sc;
}
