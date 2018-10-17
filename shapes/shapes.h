#pragma once

#include "objecttype.h"
#include "plane.h"
#include "rect.h"
#include "sphere.h"
#include "triangle.h"
#include "disk.h"
#include "instanced.h"
#include "generic.h"
#include "cylinder.h"
#include "torus.h"
#include "box.h"
#include "../aabb.h"


float rayIntersectObject(ShadeRec* sr, const Object_t obj, const Ray ray);
float shadowRayIntersectObject(const Object_t obj, const Ray ray);
bool isGridObjType(const Object_t obj);
bool getObjectAABB(AABB* aabb, const Object_t obj);
Material* getObjectMatPtr(const Object_t obj);
bool calcBoundingSphere(vec3 center, float *radius, const Object_t obj);
int testAABBPlane(AABB* aabb, vec3 plane_normal, float plane_d);
int testTriangleAABB(vec3 tv0, vec3 tv1, vec3 tv2, AABB* aabb, int* has_zero_vector, int is_debugging);
int triangleAABBIntersect(Object_t obj, AABB* aabb);
float shapePdf(const Object_t obj, const ShadeRec* sr, const vec3 wi);

 
