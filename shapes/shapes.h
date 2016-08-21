#pragma once

#include "objecttype.h"
#include "aabox.h"
#include "plane.h"
#include "rect.h"
#include "sphere.h"
#include "triangle.h"
#include "disk.h"
#include "instanced.h"
#include "generic.h"
#include "cylinder.h"
#include "torus.h"
#include "../aabb.h"


float rayIntersectObject(ShadeRec* sr, const Object_t obj, const Ray ray)
{
    float t = TMAX;
    switch(obj.type)
    {
    case SPHERE:
        t = rayIntersectSphere(sr, (Sphere*)obj.ptr, ray);
        break;
    case PLANE:
        t = rayIntersectPlane(sr, (Plane*)obj.ptr, ray);
        break;
    case RECTANGLE:
        t = rayIntersectRect(sr, (Rectangle*)obj.ptr, ray);
        break;
    case AABOX:
        t = rayIntersectAABox(sr, (AABox*)obj.ptr, ray);
        break;
    case TRIANGLE:
        t = rayIntersectTriangle(sr, (Triangle*)obj.ptr, ray);
        break;
    case GENERICOPENCYLINDER:
        t = rayIntersectGenericOpenCylinder(sr, (GenericOpenCylinder*)obj.ptr, ray);
        break;
    case DISK:
        t = rayIntersectDisk(sr, (Disk*)obj.ptr, ray);
        break;
    case GENERICTORUS:
        t = rayIntersectGenericTorus(sr, (GenericTorus*)obj.ptr, ray);
        break;
    case INSTANCED:
        t = rayIntersectInstanced(sr, (InstancedShape*)obj.ptr, ray);
        break;
    case COMPOUND:
        t = rayIntersectCompound(sr, (CompoundObject*)obj.ptr, ray);
        break;
    case FLAT_TRIANGLE:
        t = rayIntersectFlatTriangle(sr, (FlatTriangle*)obj.ptr, ray);
        break;
    case SMOOTH_TRIANGLE:
        t = rayIntersectSmoothTriangle(sr, (SmoothTriangle*)obj.ptr, ray);
        break;                
    }
    return t;
}

float shadowRayIntersectObject(const Object_t obj, const Ray ray)
{
    float t = TMAX;
    switch(obj.type)
    {
    case SPHERE:
        t = shadowRayIntersectSphere((Sphere*)obj.ptr, ray);
        break;
    case PLANE:
        t = shadowRayIntersectPlane((Plane*)obj.ptr, ray);
        break;
    case RECTANGLE:
        t = shadowRayIntersectRect((Rectangle*)obj.ptr, ray);
        break;
    case AABOX:
        t = shadowRayIntersectAABox((AABox*)obj.ptr, ray);
        break;
    case TRIANGLE:
        t = shadowRayIntersectTriangle((Triangle*)obj.ptr, ray);
        break;
    case GENERICOPENCYLINDER:
        t = shadowRayIntersectGenericOpenCylinder((GenericOpenCylinder*)obj.ptr, ray);
        break;
    case DISK:
        t = shadowRayIntersectDisk((Disk*)obj.ptr, ray);
        break;
    case GENERICTORUS:
        t = shadowRayIntersectGenericTorus((GenericTorus*)obj.ptr, ray);
        break;
    case INSTANCED:
        t = shadowRayIntersectInstanced((InstancedShape*)obj.ptr, ray);
        break;
    case COMPOUND:
        t = shadowRayIntersectCompound((CompoundObject*)obj.ptr, ray);
        break;
    case FLAT_TRIANGLE:
        t = shadowRayIntersectFlatTriangle((FlatTriangle*)obj.ptr, ray);
        break;
    case SMOOTH_TRIANGLE:
        t = shadowRayIntersectSmoothTriangle((SmoothTriangle*)obj.ptr, ray);
        break;                        
    }
    return t;
}   

bool isGridObjType(const Object_t obj)
{
    switch(obj.type)
    {
    case SPHERE:
    case RECTANGLE:
    case AABOX:
    case TRIANGLE:
    case GENERICOPENCYLINDER:
    case DISK:
    case GENERICTORUS:
    case FLAT_TRIANGLE:
    case SMOOTH_TRIANGLE:        
        return true;
        break;
    case INSTANCED:
    {
        InstancedShape* is = (InstancedShape*)(obj.ptr);
        return isGridObjType(is->obj);
    } break;
    }
    return false;
}

// TODO: add return value
void getObjectAABB(AABB* aabb, const Object_t obj)
{
    switch(obj.type)
    {
    case PLANE:
      break;
    case SPHERE:
    {
        Sphere* sphere = (Sphere*)obj.ptr;
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
        Rectangle* rect = (Rectangle*)obj.ptr;
        vec3 pts[4];
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
        AABox* aabox = (AABox*)obj.ptr;
        float gap = 0.00001f;
        vec3_assign(aabb->min, aabox->min[0] - gap, aabox->min[1] - gap, aabox->min[2] - gap);
        vec3_assign(aabb->max, aabox->max[0] + gap, aabox->max[1] + gap, aabox->max[2] + gap);        
    } break;
    case TRIANGLE:
    {
        Triangle* triangle = (Triangle*)obj.ptr;
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
    case FLAT_TRIANGLE:
    {
        FlatTriangle* mesh_tri = (FlatTriangle*)obj.ptr;
        vec3_copy(aabb->min, mesh_tri->v0);
        vec3_copy(aabb->max, mesh_tri->v0);

        for(int i = 0; i < 3; i++)
        {
            if(mesh_tri->v1[i] < aabb->min[i]){aabb->min[i] = mesh_tri->v1[i];}
            if(mesh_tri->v1[i] > aabb->max[i]){aabb->max[i] = mesh_tri->v1[i];}
            if(mesh_tri->v2[i] < aabb->min[i]){aabb->min[i] = mesh_tri->v2[i];}
            if(mesh_tri->v2[i] > aabb->max[i]){aabb->max[i] = mesh_tri->v2[i];}
        }
        for(int i = 0; i < 3; i++)
        {
            aabb->min[i] -= K_FLAT_AABB;
            aabb->max[i] += K_FLAT_AABB;
        }
    } break;
    case SMOOTH_TRIANGLE:
    {
        SmoothTriangle* mesh_tri = (SmoothTriangle*)obj.ptr;
        vec3_copy(aabb->min, mesh_tri->v0);
        vec3_copy(aabb->max, mesh_tri->v0);

        for(int i = 0; i < 3; i++)
        {
            if(mesh_tri->v1[i] < aabb->min[i]){aabb->min[i] = mesh_tri->v1[i];}
            if(mesh_tri->v1[i] > aabb->max[i]){aabb->max[i] = mesh_tri->v1[i];}
            if(mesh_tri->v2[i] < aabb->min[i]){aabb->min[i] = mesh_tri->v2[i];}
            if(mesh_tri->v2[i] > aabb->max[i]){aabb->max[i] = mesh_tri->v2[i];}
        }
        for(int i = 0; i < 3; i++)
        {
            aabb->min[i] -= K_FLAT_AABB;
            aabb->max[i] += K_FLAT_AABB;
        }
    } break;        
    case GENERICOPENCYLINDER:
    {
        GenericOpenCylinder* oc = (GenericOpenCylinder*)obj.ptr;
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
        Disk* disk = (Disk*)obj.ptr;
        aabb->min[0] = disk->center[0] - disk->radius;
        aabb->min[1] = disk->center[1] - disk->radius;
        aabb->min[2] = disk->center[2] - disk->radius;

        aabb->max[0] = disk->center[0] + disk->radius;
        aabb->max[1] = disk->center[1] + disk->radius;
        aabb->max[2] = disk->center[2] + disk->radius;        
    } break;
    case GENERICTORUS:
    {
        GenericTorus* torus = (GenericTorus*)obj.ptr;
        //*aabb = torus->aabb;
        calcAABBGenericTorus(aabb, torus);
    } break;
    case INSTANCED:
    {
        InstancedShape* is = (InstancedShape*)obj.ptr;
        getObjectAABB(aabb, is->obj);
        mat4 transform;
        affine_inverse(transform, is->inv_transform);
        vec3 tmp1, tmp2;
        vec4_assign(tmp1, aabb->min[0], aabb->min[1], aabb->min[2], 1.0f);
        mat4_mult_vec4(tmp2, transform, tmp1);
        vec3_assign(aabb->min, tmp2[0], tmp2[1], tmp2[2]);
        vec4_assign(tmp1, aabb->max[0], aabb->max[1], aabb->max[2], 1.0f);
        mat4_mult_vec4(tmp2, transform, tmp1);
        vec3_assign(aabb->max, tmp2[0], tmp2[1], tmp2[2]);        
    } break;
    case COMPOUND:
    {
        CompoundObject* co = (CompoundObject*)obj.ptr;
        *aabb = co->aabb;
    } break;
    }
}
