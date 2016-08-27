#pragma once

#include "../util/util.h"
#include "../util/vec.h"
#include "../util/ray.h"
#include "../util/shaderec.h"
#include "../util/constants.h"
#include "../materials.h"

typedef struct Sphere_s
{
    bool shadow;
    bool partial;        // NOTE: Maybe worth it, maybe not
    float min_theta, max_theta, phi; // For setting how much of the sphere is visible if partial is true
    float radius;    
    vec3 center;
    Material* mat;
} Sphere;

void fillShadeRecSphere(ShadeRec *sr, Sphere *sphere, const vec3 hit_point, const Ray ray, const float t);
float rayIntersectSphere(ShadeRec *sr, Sphere *sphere, const Ray ray);
float shadowRayIntersectSphere(Sphere *sphere, const Ray ray);
