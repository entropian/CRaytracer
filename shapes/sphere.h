#pragma once

#include "../util/util.h"
#include "../util/vec.h"
#include "../util/ray.h"
#include "../util/shaderec.h"
#include "../util/constants.h"
#include "../materials.h"

typedef struct Sphere_s
{
    float min_theta, max_theta, phi; // For setting how much of the sphere is visible if partial is true
    float radius;        
    Material* mat;
    //mat4 world_to_object, object_to_world;
    vec3 center;
} Sphere;

void fillShadeRecSphere(ShadeRec *sr, Sphere *sphere, const vec3 hit_point, const Ray ray, const float t,
                        const float theta, const float phi);
float rayIntersectSphere(ShadeRec *sr, Sphere *sphere, const Ray ray);
float shadowRayIntersectSphere(Sphere *sphere, const Ray ray);
