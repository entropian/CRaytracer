#pragma once

#include "../util/mat.h"
#include "objecttype.h"
#include "../aabb.h"
#include "generic.h"

typedef struct CompoundObject
{
    bool shadow;
    int num_obj;
    Object_t* objects;
    AABB aabb;    
} CompoundObject;

float rayIntersectCompound(ShadeRec* sr, CompoundObject* co, const Ray ray);
float shadowRayIntersectCompound(CompoundObject* co, const Ray shadow_ray);

typedef struct InstancedShape
{
    Object_t obj;
    mat4 inv_transform;
    Material* mat;        
} InstancedShape;

void freeInstancedShape(InstancedShape* is);

void freeCompoundObject(CompoundObject* co)
{
    for(int i = 0; i < co->num_obj; i++)
    {
        if(co->objects[i].type == INSTANCED)
        {
            freeInstancedShape((InstancedShape*)(co->objects[i].ptr));
        }else if(co->objects[i].type == COMPOUND)
        {
            freeCompoundObject((CompoundObject*)(co->objects[i].ptr));
        }else
        {
            free(co->objects[i].ptr);
        }
    }
    free(co->objects);
}

void freeInstancedShape(InstancedShape* is)
{
    if(is->obj.type == COMPOUND)
    {
        freeCompoundObject((CompoundObject*)is->obj.ptr);
    }else
    {
        free(is->obj.ptr);
    }
}

float rayIntersectInstanced(ShadeRec* sr, InstancedShape* is, const Ray src_ray)
{
    Ray ray;    
    transformRay(&ray, is->inv_transform, src_ray);
    float t = TMAX;
    switch(is->obj.type)
    {
    case SPHERE:
        t = rayIntersectSphere(sr, (Sphere*)is->obj.ptr, ray);            
        break;
    case PLANE:
        t = rayIntersectPlane(sr, (Plane*)is->obj.ptr, ray);                        
        break;
    case RECTANGLE:
        t = rayIntersectRect(sr, (Rectangle*)is->obj.ptr, ray);
        break;
    case AABOX:
        t = rayIntersectAABox(sr, (AABox*)is->obj.ptr, ray);
        break;
    case TRIANGLE:
        t = rayIntersectTriangle(sr, (Triangle*)is->obj.ptr, ray);
        break;
    case MESH_TRIANGLE:
      t = rayIntersectMeshTriangle(sr, (MeshTriangle*)is->obj.ptr, ray);
      break;
    case GENERICOPENCYLINDER:
        t = rayIntersectGenericOpenCylinder(sr, (GenericOpenCylinder*)is->obj.ptr, ray);            
        break;
    case DISK:
        t = rayIntersectDisk(sr, (Disk*)is->obj.ptr, ray);
        break;
    case GENERICTORUS:
        t = rayIntersectGenericTorus(sr, (GenericTorus*)is->obj.ptr, ray);
        break;
    case INSTANCED:
        t = rayIntersectInstanced(sr, (InstancedShape*)is->obj.ptr, ray);
        break;
    case COMPOUND:
        t = rayIntersectCompound(sr, (CompoundObject*)is->obj.ptr, ray);
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
    switch(is->obj.type)
    {
        case SPHERE:
            t = shadowRayIntersectSphere((Sphere*)is->obj.ptr, shadow_ray);            
            break;
        case PLANE:
            t = shadowRayIntersectPlane((Plane*)is->obj.ptr, shadow_ray);                        
            break;
        case RECTANGLE:
            t = shadowRayIntersectRect((Rectangle*)is->obj.ptr, shadow_ray);
            break;
        case AABOX:
            t = shadowRayIntersectAABox((AABox*)is->obj.ptr, shadow_ray);
            break;
        case TRIANGLE:
            t = shadowRayIntersectTriangle((Triangle*)is->obj.ptr, shadow_ray);
            break;
        case MESH_TRIANGLE:
            t = shadowRayIntersectMeshTriangle((MeshTriangle*)is->obj.ptr, shadow_ray);
            break;
        case GENERICOPENCYLINDER:
            t = shadowRayIntersectGenericOpenCylinder((GenericOpenCylinder*)is->obj.ptr, shadow_ray);            
            break;
        case DISK:
            t = shadowRayIntersectDisk((Disk*)is->obj.ptr, shadow_ray);
            break;
        case GENERICTORUS:
            t = shadowRayIntersectGenericTorus((GenericTorus*)is->obj.ptr, shadow_ray);
            break;
        case INSTANCED:
            t = shadowRayIntersectInstanced((InstancedShape*)is->obj.ptr, shadow_ray);
            break;
        case COMPOUND:
            t = shadowRayIntersectCompound((CompoundObject*)is->obj.ptr, shadow_ray);
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
        void* obj_ptr = co->objects[i].ptr;
        switch(co->objects[i].type)
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
        case MESH_TRIANGLE:
            tmp_t = rayIntersectMeshTriangle(&tmp_sr, (MeshTriangle*)obj_ptr, ray);
            break;
        case GENERICOPENCYLINDER:
            tmp_t = rayIntersectGenericOpenCylinder(&tmp_sr, (GenericOpenCylinder*)obj_ptr, ray);            
            break;
        case DISK:
            tmp_t = rayIntersectDisk(&tmp_sr, (Disk*)obj_ptr, ray);
            break;
        case GENERICTORUS:
            tmp_t = rayIntersectGenericTorus(&tmp_sr, (GenericTorus*)obj_ptr, ray);
            break;
        case INSTANCED:
            tmp_t = rayIntersectInstanced(&tmp_sr, (InstancedShape*)obj_ptr, ray);
            break;
	case COMPOUND:
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
        void* obj_ptr = co->objects[i].ptr;
        switch(co->objects[i].type)
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
        case MESH_TRIANGLE:
            t = shadowRayIntersectMeshTriangle((MeshTriangle*)obj_ptr, shadow_ray);
            break;
        case GENERICOPENCYLINDER:
            t = shadowRayIntersectGenericOpenCylinder((GenericOpenCylinder*)obj_ptr, shadow_ray);            
            break;
        case DISK:
            t = shadowRayIntersectDisk((Disk*)obj_ptr, shadow_ray);
            break;
        case GENERICTORUS:
            t = shadowRayIntersectGenericTorus((GenericTorus*)obj_ptr, shadow_ray);
            break;
        case INSTANCED:
            t = shadowRayIntersectInstanced((InstancedShape*)obj_ptr, shadow_ray);
            break;
	case COMPOUND:
	  break;
        }
        if(t < TMAX)
        {
            return t;
        }
    }
    return t;
}


