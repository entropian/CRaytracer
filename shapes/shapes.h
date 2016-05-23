#pragma once

#include "objecttype.h"
#include "aabox.h"
#include "plane.h"
#include "rect.h"
#include "sphere.h"
#include "triangle.h"
#include "cylinder.h"
#include "disk.h"
#include "torus.h"
#include "instanced.h"
#include "../aabb.h"

float rayIntersectObject(ShadeRec* sr, const void* obj_ptr, const ObjectType obj_type, const Ray ray)
{
    float t = TMAX;
    switch(obj_type)
    {
    case SPHERE:
        t = rayIntersectSphere(sr, (Sphere*)obj_ptr, ray);
        break;
    case PLANE:
        t = rayIntersectPlane(sr, (Plane*)obj_ptr, ray);
        break;
    case RECTANGLE:
        t = rayIntersectRect(sr, (Rectangle*)obj_ptr, ray);
        break;
    case AABOX:
        t = rayIntersectAABox(sr, (AABox*)obj_ptr, ray);
        break;
    case TRIANGLE:
        t = rayIntersectTriangle(sr, (Triangle*)obj_ptr, ray);
        break;
    case OPENCYLINDER:
        t = rayIntersectOpenCylinder(sr, (OpenCylinder*)obj_ptr, ray);
        break;
    case DISK:
        t = rayIntersectDisk(sr, (Disk*)obj_ptr, ray);
        break;
    case TORUS:
        t = rayIntersectTorus(sr, (Torus*)obj_ptr, ray);
        break;
    case INSTANCED:
        t = rayIntersectInstanced(sr, (InstancedShape*)obj_ptr, ray);
        break;
    case COMPOUND:
        t = rayIntersectCompound(sr, (CompoundObject*)obj_ptr, ray);
        break;
    }
    return t;
}

float shadowRayIntersectObject(const void* obj_ptr, const ObjectType obj_type, const Ray ray)
{
    float t = TMAX;
    switch(obj_type)
    {
    case SPHERE:
        t = shadowRayIntersectSphere((Sphere*)obj_ptr, ray);
        break;
    case PLANE:
        t = shadowRayIntersectPlane((Plane*)obj_ptr, ray);
        break;
    case RECTANGLE:
        t = shadowRayIntersectRect((Rectangle*)obj_ptr, ray);
        break;
    case AABOX:
        t = shadowRayIntersectAABox((AABox*)obj_ptr, ray);
        break;
    case TRIANGLE:
        t = shadowRayIntersectTriangle((Triangle*)obj_ptr, ray);
        break;
    case OPENCYLINDER:
        t = shadowRayIntersectOpenCylinder((OpenCylinder*)obj_ptr, ray);
        break;
    case DISK:
        t = shadowRayIntersectDisk((Disk*)obj_ptr, ray);
        break;
    case TORUS:
        t = shadowRayIntersectTorus((Torus*)obj_ptr, ray);
        break;
    case INSTANCED:
        t = shadowRayIntersectInstanced((InstancedShape*)obj_ptr, ray);
        break;
    case COMPOUND:
        t = shadowRayIntersectCompound((CompoundObject*)obj_ptr, ray);
        break;
    }
    return t;
}   



// TODO: add return value
void getObjectAABB(AABB* aabb, const void* obj_ptr, const ObjectType obj_type)
{
    switch(obj_type)
    {
    case SPHERE:
    {
        Sphere* sphere = (Sphere*)obj_ptr;
        aabb->min[0] = sphere->center[0] - sphere->radius;
        aabb->min[1] = sphere->center[1] - sphere->radius;
        aabb->min[2] = sphere->center[2] - sphere->radius;

        aabb->max[0] = sphere->center[0] + sphere->radius;
        aabb->max[1] = sphere->center[1] + sphere->radius;
        aabb->max[2] = sphere->center[2] + sphere->radius;
    } break;
    case RECTANGLE:
    {
        // NOTE: shit code
        Rectangle* rect = (Rectangle*)obj_ptr;
        vec3 pts[4];
        vec3 pt1, pt2, pt3, pt4;
        vec3_copy(pts[0], rect->point);
        vec3_add(pts[1], pts[0], rect->width);
        vec3_add(pts[2], pts[0], rect->height);
        vec3 tmp;
        vec3_add(tmp, rect->width, rect->height);
        vec3_add(pts[3], pts[0], tmp);
            
        vec3_copy(aabb->min, pts[0]);
        vec3_copy(aabb->max, pts[0]);

        for(int i = 1; i < 4; i++)
        {
            for(int j = 0; j < 3; j++)
            {
                if(pts[i][j] < aabb->min[j]){aabb->min[j] = pts[i][j] - K_EPSILON;}
                if(pts[i][j] > aabb->max[j]){aabb->max[j] = pts[i][j] + K_EPSILON;}
            }
        }
    } break;
    case AABOX:
    {
        AABox* aabox = (AABox*)obj_ptr;
        vec3_copy(aabb->min, aabox->min);
        vec3_copy(aabb->max, aabox->max);            
    } break;
    case TRIANGLE:
    {
        Triangle* triangle = (Triangle*)obj_ptr;
        vec3_copy(aabb->min, triangle->v0);
        vec3_copy(aabb->max, triangle->v0);

        for(int i = 0; i < 3; i++)
        {
            if(triangle->v1[i] < aabb->min[i]){aabb->min[i] = triangle->v1[i];}
            if(triangle->v1[i] > aabb->max[i]){aabb->max[i] = triangle->v1[i];}
            if(triangle->v2[i] < aabb->min[i]){aabb->min[i] = triangle->v2[i];}
            if(triangle->v2[i] > aabb->max[i]){aabb->max[i] = triangle->v2[i];}
        }
        for(int i = 0; i < 3; i++)
        {
            aabb->min[i] -= K_FLAT_AABB;
            aabb->max[i] += K_FLAT_AABB;
        }
    } break;
    case OPENCYLINDER:
    {
        OpenCylinder* oc = (OpenCylinder*)obj_ptr;
        aabb->min[0] = -(oc->radius);
        aabb->min[1] = -(oc->half_height);
        aabb->min[2] = -(oc->radius);
        aabb->max[0] = oc->radius;
        aabb->max[1] = oc->half_height;
        aabb->max[2] = oc->radius;
    } break;
    case DISK:
    {
        // TODO: improve this crap
        Disk* disk = (Disk*)obj_ptr;
        aabb->min[0] = disk->center[0] - disk->radius;
        aabb->min[1] = disk->center[1] - disk->radius;
        aabb->min[2] = disk->center[2] - disk->radius;

        aabb->max[0] = disk->center[0] + disk->radius;
        aabb->max[1] = disk->center[1] + disk->radius;
        aabb->max[2] = disk->center[2] + disk->radius;        
    } break;
    case TORUS:
    {
        Torus* torus = (Torus*)obj_ptr;
        *aabb = torus->aabb;
    } break;
    case INSTANCED:
    {
        InstancedShape* is = (InstancedShape*)obj_ptr;
        getObjectAABB(aabb, is->obj_ptr, is->obj_type);
    } break;
    case COMPOUND:
    {
        CompoundObject* co = (CompoundObject*)obj_ptr;
        *aabb = co->aabb;
    } break;
    }
}


