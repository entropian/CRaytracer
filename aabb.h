#pragma once

#include "util/constants.h"
#include "util/vec.h"

typedef struct AABB
{
    vec3 min;
    vec3 max;
} AABB;

bool rayIntersectAABB(const AABB* aabb, const Ray ray)
{
    float tx_min, ty_min, tz_min;
    float tx_max, ty_max, tz_max;

    float a = 1.0f/ray.direction[0];
    if(a >= 0)
    {
        tx_min = (aabb->min[0] - ray.origin[0]) * a;
        tx_max = (aabb->max[0] - ray.origin[0]) * a;
    }else
    {
        tx_min = (aabb->max[0] - ray.origin[0]) * a;
        tx_max = (aabb->min[0] - ray.origin[0]) * a;        
    }

    float b = 1.0f/ray.direction[1];
    if(b >= 0)
    {
        ty_min = (aabb->min[1] - ray.origin[1]) * b;
        ty_max = (aabb->max[1] - ray.origin[1]) * b;        
    }else
    {
        ty_min = (aabb->max[1] - ray.origin[1]) * b;
        ty_max = (aabb->min[1] - ray.origin[1]) * b;        
    }

    float c = 1.0f/-ray.direction[2];
    if(c >= 0)
    {
        tz_min = -(aabb->min[2] - ray.origin[2]) * c;
        tz_max = -(aabb->max[2] - ray.origin[2]) * c;
    }else
    {
        tz_min = -(aabb->max[2] - ray.origin[2]) * c;
        tz_max = -(aabb->min[2] - ray.origin[2]) * c;
    }
    
    float t0, t1;
    if(tx_min > ty_min)
    {
        t0 = tx_min;
    }else
    {
        t0 = ty_min;
    }

    if(tz_min > t0)
    {
        t0 = tz_min;
    }

    if(tx_max < ty_max)
    {
        t1 = tx_max;
    }else
    {
        t1 = ty_max;
    }

    if(tz_max < t1)
    {
        t1 = tz_max;
    }

    return (t0 < t1 && t1 > K_EPSILON);
}
