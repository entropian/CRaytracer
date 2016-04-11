#pragma once

#include "../vec.h"
#include "../materials.h"
#include "../constants.h"
#include "../ray.h"
#include "../shaderec.h"

typedef struct Plane
{
    bool shadow;
    vec3 point;
    vec3 normal;
    Material mat;
} Plane;

float rayIntersectPlane(ShadeRec *sr, Plane *plane, const Ray ray)
{
    vec3 displacement;
    vec3_sub(displacement, plane->point, ray.origin);
    float t = vec3_dot(displacement, plane->normal) / vec3_dot(ray.direction, plane->normal);
    if(t > k_epsilon)
    {
        vec3_copy(sr->normal, plane->normal);
        vec3_scale(displacement, ray.direction, t);
        vec3_add(sr->hit_point, ray.origin, displacement);
        vec3_negate(sr->wo, ray.direction);
        sr->mat = &(plane->mat);
        return t;
    }
    return TMAX;
}

float shadowRayIntersectPlane(Plane *plane, const Ray ray)
{
    if(!plane->shadow)
    {
        return TMAX;
    }
    vec3 displacement;
    vec3_sub(displacement, plane->point, ray.origin);
    float t = vec3_dot(displacement, plane->normal) / vec3_dot(ray.direction, plane->normal);
    if(t > k_epsilon)
    {
        return t;
    }
    return TMAX;
}
