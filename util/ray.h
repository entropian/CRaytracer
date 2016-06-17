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

void calcReflectRayDir(vec3 reflect_dir, const vec3 normal, const vec3 primary_dir)
{
    float magnitude = -2.0f * vec3_dot(primary_dir, normal);
    vec3 adjusted_normal;
    vec3_scale(adjusted_normal, normal, magnitude);
    vec3_add(reflect_dir, primary_dir, adjusted_normal);
}
