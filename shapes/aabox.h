#pragma once

#include "../util/vec.h"
#include "../util/ray.h"
#include "../util/shaderec.h"
#include "../util/constants.h"
#include "../materials.h"

typedef struct AABox
{
    bool shadow;
    vec3 min, max;
    Material* mat;
} AABox;

void getAABoxNormal(vec3 r, const int face_hit);
float rayIntersectAABox(ShadeRec* sr, AABox* aabox, const Ray ray);
float shadowRayIntersectAABox(const AABox* aabox, const Ray ray);
