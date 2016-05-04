#pragma once

#include "../util/util.h"
#include "../util/vec.h"
#include "../util/ray.h"
#include "../util/shaderec.h"
#include "../util/constants.h"
#include "../materials.h"
#include "../aabb.h"
#include "disk.h"

enum NormalType
{
    OPEN,
    CONVEX,
    CONCAVE
};

typedef struct OpenCylinder
{
    bool shadow;
    NormalType normal_type;
    float half_height;
    float radius;
    float phi; // for setting how much the cylinder is visible
    Material mat;    
} OpenCylinder;

void fillShadeRecOpenCylinder(ShadeRec *sr, OpenCylinder* oc, const Ray ray, const vec3 hit_point)
{
    sr->hit_status = true;
    vec3_negate(sr->wo, ray.direction);
    vec3_copy(sr->hit_point, hit_point);
    switch(oc->normal_type)
    {
    case OPEN:
    {
        vec3_assign(sr->normal, hit_point[0]/oc->radius, 0.0f, hit_point[2]/oc->radius);
        if(vec3_dot(sr->wo, sr->normal) < 0)
        {
            vec3_negate(sr->normal, sr->normal);
        }
    } break;
    case CONVEX:
    {
        vec3_assign(sr->normal, hit_point[0]/oc->radius, 0.0f, hit_point[2]/oc->radius);        
    } break;
    case CONCAVE:
    {
        vec3_assign(sr->normal, (-hit_point[0])/oc->radius, 0.0f, (-hit_point[2])/oc->radius);
    } break;
    }
    sr->mat = &(oc->mat);
}

float rayIntersectOpenCylinder(ShadeRec* sr, OpenCylinder* oc, const Ray ray)
{
    // intersection equation: at^2 + bt + c = 0 quadratic
    float a = ray.direction[0]*ray.direction[0] + ray.direction[2]*ray.direction[2];
    float b = 2*(ray.origin[0]*ray.direction[0] + ray.origin[2]*ray.direction[2]);
    float c = ray.origin[0]*ray.origin[0] + ray.origin[2]*ray.origin[2] - oc->radius*oc->radius;
    float disc = b*b - 4*a*c;

    if(disc < 0.0f)
    {
        return TMAX;
    }else
    {
        float e = sqrt(disc);
        float denom = 2*a;
        float t = (-b - e)/denom;
        vec3 hit_point;        
        if(t > K_EPSILON)
        {
            getPointOnRay(hit_point, ray, t);
            if(abs(hit_point[1]) <= oc->half_height)
            {
                float phi = (float)atan2(hit_point[0], hit_point[2]);
                if(abs(phi) <= oc->phi)
                {
                    fillShadeRecOpenCylinder(sr, oc, ray, hit_point);
                    return t;                                        
                }
            }
        }
        t = (-b + e)/denom;
        if(t > K_EPSILON)
        {
            getPointOnRay(hit_point, ray, t);
            if(abs(hit_point[1]) <= oc->half_height)
            {
                float phi = (float)atan2(hit_point[0], hit_point[2]);
                if(abs(phi) <= oc->phi)
                {
                    fillShadeRecOpenCylinder(sr, oc, ray, hit_point);
                    return t;                    
                }
            }
        }
    }
    return TMAX;
}

float shadowRayIntersectOpenCylinder(const OpenCylinder* oc, const Ray ray)
{
    float a = ray.direction[0]*ray.direction[0] + ray.direction[2]*ray.direction[2];
    float b = 2*(ray.origin[0]*ray.direction[0] + ray.origin[2]*ray.direction[2]);
    float c = ray.origin[0]*ray.origin[0] + ray.origin[2]*ray.origin[2] - oc->radius*oc->radius;
    float disc = b*b - 4*a*c;

    if(disc < 0.0f)
    {
        return TMAX;
    }else
    {
        float e = sqrt(disc);
        float denom = 2*a;
        float t = (-b - e)/denom;
        vec3 hit_point;        
        if(t > K_EPSILON)
        {
            getPointOnRay(hit_point, ray, t);            
            if(abs(hit_point[1]) <= oc->half_height)
            {
                float phi = (float)atan2(hit_point[0], hit_point[2]);
                if(abs(phi) <= oc->phi)
                {
                    return t;
                }
            }
        }
        t = (-b + e)/denom;
        if(t > K_EPSILON)
        {
            getPointOnRay(hit_point, ray, t);            
            if(abs(hit_point[1]) <= oc->half_height)
            {
                float phi = (float)atan2(hit_point[0], hit_point[2]);
                if(abs(phi) <= oc->phi)
                {
                    return t;
                }
            }
        }
    }
    return TMAX;
}

typedef struct SolidCylinder
{
    bool shadow;
    Disk* top;
    Disk* bottom;
    OpenCylinder* tube;
    AABB aabb;
} SolidCylinder;

void calcAABBSolidCylinder(SolidCylinder* sc)
{
    AABB* aabb = &(sc->aabb);
    aabb->min[0] = -(sc->tube->radius);
    aabb->min[1] = -(sc->tube->half_height);
    aabb->min[2] = sc->tube->radius;

    aabb->max[0] = sc->tube->radius;
    aabb->max[1] = sc->tube->half_height;
    aabb->max[2] = -(sc->tube->radius);
}

SolidCylinder* initSolidCylinder(const float radius, const float half_height,
                       const float phi, const Material* mat, bool shadow)
{
    SolidCylinder* sc = (SolidCylinder*)malloc(sizeof(SolidCylinder));
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
    top->mat = bottom->mat = tube->mat = *mat;
    sc->top = top;
    sc->bottom = bottom;
    sc->tube = tube;
    calcAABBSolidCylinder(sc);
    return sc;
}

void freeSolidCylinder(SolidCylinder* sc)
{
    if(sc != NULL)
    {
        if(sc->top){free(sc->top);}
        if(sc->bottom){free(sc->bottom);}
        if(sc->tube){free(sc->tube);}
        free(sc);
    }
}

float rayIntersectSolidCylinder(ShadeRec* sr, SolidCylinder* sc, const Ray ray)
{
    if(!rayIntersectAABB(&(sc->aabb), ray))
    {
        return TMAX;
    }
    float min_t = TMAX, tmp_t = TMAX;
    ShadeRec min_sr, tmp_sr;

    min_t = rayIntersectDisk(&min_sr, sc->top, ray);
    tmp_t = rayIntersectDisk(&tmp_sr, sc->bottom, ray);
    if(tmp_t < min_t)
    {
        min_t = tmp_t;
        min_sr = tmp_sr;
    }
    tmp_t = rayIntersectOpenCylinder(&tmp_sr, sc->tube, ray);
    if(tmp_t < min_t)
    {
        min_t = tmp_t;
        min_sr = tmp_sr;
    }
    *sr = min_sr;
    return min_t;
}

float shadowRayIntersectSolidCylinder(const SolidCylinder* sc, const Ray ray)
{
    if(!rayIntersectAABB(&(sc->aabb), ray))
    {
        return TMAX;
    }    
    float min_t = TMAX, tmp_t = TMAX;

    min_t = shadowRayIntersectDisk(sc->top, ray);
    tmp_t = shadowRayIntersectDisk(sc->bottom, ray);
    if(tmp_t < min_t)
    {
        min_t = tmp_t;
    }
    tmp_t = shadowRayIntersectOpenCylinder(sc->tube, ray);
    if(tmp_t < min_t)
    {
        min_t = tmp_t;
    }
    return min_t;
}
