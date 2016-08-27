#pragma once

#include "objecttype.h"
#include "aabox.h"
#include "plane.h"
#include "rect.h"
#include "sphere.h"
#include "triangle.h"
#include "disk.h"
#include "instanced.h"
#include "generic.h"
#include "cylinder.h"
#include "torus.h"
#include "../aabb.h"


float rayIntersectObject(ShadeRec* sr, const Object_t obj, const Ray ray);
float shadowRayIntersectObject(const Object_t obj, const Ray ray);
bool isGridObjType(const Object_t obj);
void getObjectAABB(AABB* aabb, const Object_t obj);
