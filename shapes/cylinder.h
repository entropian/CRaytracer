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
    Material* mat;    
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
    sr->mat = oc->mat;
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
            if(abs((float)hit_point[1]) <= oc->half_height)
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
