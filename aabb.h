#pragma once

#include "util/constants.h"
#include "util/vec.h"
#include "util/ray.h"

typedef struct AABB
{
    vec3 min;
    vec3 max;
} AABB;

float rayIntersectAABB(const AABB* aabb, const Ray ray);
bool isInsideAABB(const AABB* aabb, const vec3 point);
void addToAABB(AABB* r, const AABB* a);
int testAABB(AABB a, AABB b);
