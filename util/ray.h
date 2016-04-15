#pragma once

#include "vec.h"

typedef struct Ray
{
    vec3 origin;
    vec3 direction;
} Ray;

void getPointOnRay(vec3 r, const Ray ray, const float t)
{
    vec3 displacement;
    vec3_scale(displacement, ray.direction, t);
    vec3_add(r, ray.origin, displacement);
}
