#pragma once

#include "../util/vec.h"
#include "../util/ray.h"
#include "../util/shaderec.h"
#include "../materials.h"

typedef struct Plane
{
    bool shadow;
    vec3 point;
    vec3 normal;
    Material* mat;
} Plane;

float rayIntersectPlane(ShadeRec *sr, Plane *plane, const Ray ray);
float shadowRayIntersectPlane(Plane *plane, const Ray ray);
