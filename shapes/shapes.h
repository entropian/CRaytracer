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






