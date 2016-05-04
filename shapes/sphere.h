#pragma once

#include "../util/util.h"
#include "../util/vec.h"
#include "../util/ray.h"
#include "../util/shaderec.h"
#include "../util/constants.h"
#include "../materials.h"

typedef struct Sphere
{
    bool shadow;
    bool partial;        // NOTE: Maybe worth it, maybe not
    float min_theta, max_theta, phi; // For setting how much of the sphere is visible if partial is true
    float radius;    
    vec3 center;
    Material mat;
} Sphere;

void fillShadeRecSphere(ShadeRec *sr, Sphere *sphere, const vec3 hit_point, const Ray ray, const float t)
{
    sr->hit_status = true;
    vec3_copy(sr->hit_point, hit_point);
    vec3 hit_point_to_center;
    vec3_sub(hit_point_to_center, sr->hit_point, sphere->center);
    vec3_normalize(sr->normal, hit_point_to_center);
    vec3_negate(sr->wo, ray.direction);
    sr->mat = &(sphere->mat);
}

float rayIntersectSphere(ShadeRec *sr, Sphere *sphere, const Ray ray)
{
    // The analytic solution is so much better than the shitty loop from before
    // though I probably should have used smaller increment for the loop
    float a = vec3_dot(ray.direction, ray.direction);
    vec3 origin_to_center, tmp;
    vec3_sub(origin_to_center, ray.origin, sphere->center);
    vec3_scale(tmp, origin_to_center, 2);
    float b = vec3_dot(tmp, ray.direction);
    float c = vec3_dot(origin_to_center, origin_to_center) - sphere->radius*sphere->radius;
    float disc = b*b - 4*a*c;

    if(disc < 0.0f)
    {
        return TMAX;
    }else
    {
        float e = sqrt(disc);
        float denom = 2*a;
        float t = (-b - e)/denom;
        if(t > K_EPSILON)
        {
            vec3 hit_point;            
            getPointOnRay(hit_point, ray, t);                        
            float phi = (float)atan2(hit_point[0] - sphere->center[0], hit_point[2] - sphere->center[2]);
            float theta = (float)acos(((hit_point[1] - sphere->center[1]) / sphere->radius));
            if(abs(phi) <= sphere->phi && theta >= sphere->min_theta&&
               theta <= sphere->max_theta)            
            {
                fillShadeRecSphere(sr, sphere, hit_point, ray, t);            
                return t;
            }
        }

        t = (-b + e)/denom;
        if(t > K_EPSILON)
        {
            vec3 hit_point;            
            getPointOnRay(hit_point, ray, t);
            float phi = (float)atan2(hit_point[0] - sphere->center[0], hit_point[2] - sphere->center[2]);
            float theta = (float)acos(((hit_point[1] - sphere->center[1]) / sphere->radius));            
            if(abs(phi) <= sphere->phi && theta >= sphere->min_theta&&
               theta <= sphere->max_theta)                            
            {
                fillShadeRecSphere(sr, sphere, hit_point, ray, t);
                return t;
            }
        }
    }
    return TMAX;
}

float shadowRayIntersectSphere(Sphere *sphere, const Ray ray)
{
    if(!sphere->shadow)
    {
        return TMAX;
    }
    float a = vec3_dot(ray.direction, ray.direction);
    vec3 origin_to_center, tmp;
    vec3_sub(origin_to_center, ray.origin, sphere->center);
    vec3_scale(tmp, origin_to_center, 2);
    float b = vec3_dot(tmp, ray.direction);
    float c = vec3_dot(origin_to_center, origin_to_center) - sphere->radius*sphere->radius;
    float disc = b*b - 4*a*c;

    if(disc < 0.0f)
    {
        return TMAX;
    }else
    {
        float e = sqrt(disc);
        float denom = 2*a;
        float t = (-b - e)/denom;
        if(t > K_EPSILON)
        {
            vec3 hit_point;            
            getPointOnRay(hit_point, ray, t);
            float phi = (float)atan2(hit_point[0] - sphere->center[0], hit_point[2] - sphere->center[2]);
            float theta = (float)acos(((hit_point[1] - sphere->center[1]) / sphere->radius));
            if(abs(phi) <= sphere->phi && theta >= sphere->min_theta&&
               theta <= sphere->max_theta)                        
            {
                return t;
            }
        }

        t = (-b + e)/denom;
        if(t > K_EPSILON)
        {
            vec3 hit_point;            
            getPointOnRay(hit_point, ray, t);
            float phi = (float)atan2(hit_point[0] - sphere->center[0], hit_point[2] - sphere->center[2]);
            float theta = (float)acos(((hit_point[1] - sphere->center[1]) / sphere->radius));
            if(abs(phi) <= sphere->phi && theta >= sphere->min_theta&&
               theta <= sphere->max_theta)                        
            {
                return t;
            }            
        }
    }
    return TMAX;
}
