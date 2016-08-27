#pragma once

#include "../util/util.h"
#include "../util/vec.h"
#include "../util/ray.h"
#include "../util/shaderec.h"
#include "../util/constants.h"
#include "../util/math.h"
#include "../materials.h"
#include "../aabb.h"
#include "disk.h"

enum NormalType
{
    OPEN,
    CONVEX,
    CONCAVE
};

typedef struct GenericOpenCylinder
{
    bool shadow;
    NormalType normal_type;
    float half_height;
    float radius;
    float phi; // for setting how much the cylinder is visible
    Material* mat;    
} GenericOpenCylinder;

void fillShadeRecGenericOpenCylinder(ShadeRec *sr, GenericOpenCylinder* oc, const Ray ray, const vec3 hit_point);
float rayIntersectGenericOpenCylinder(ShadeRec* sr, GenericOpenCylinder* oc, const Ray ray);
float shadowRayIntersectGenericOpenCylinder(const GenericOpenCylinder* oc, const Ray ray);

typedef struct GenericTorus
{
    bool shadow;
    float swept_radius;
    float tube_radius;
    float phi;
    Material* mat;
} GenericTorus;

void computeGenericTorusNormal(vec3 r, const GenericTorus* torus, const vec3 hit_point);
void calcAABBGenericTorus(AABB* aabb, GenericTorus* torus);
float rayIntersectGenericTorus(ShadeRec* sr, GenericTorus* torus, const Ray ray);
float shadowRayIntersectGenericTorus(const GenericTorus* torus, const Ray ray);
