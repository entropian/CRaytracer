#pragma once

#include "../util/vec.h"
#include "../util/ray.h"
#include "../util/shaderec.h"
#include "../util/constants.h"
#include "../materials.h"

typedef struct Rectangle
{
    bool shadow;
    vec3 point;                // Bottom left of the rectangle
    vec3 width, height;
    vec3 normal;
    Material* mat;
} Rectangle;

float rayIntersectRect(ShadeRec *sr, Rectangle *rect, const Ray ray);
float shadowRayIntersectRect(Rectangle *rect, const Ray ray);
