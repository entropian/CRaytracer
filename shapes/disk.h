#pragma once

#include "../util/vec.h"
#include "../util/ray.h"
#include "../util/shaderec.h"
#include "../util/constants.h"
#include "../materials.h"

typedef struct Disk
{
    bool shadow;
    float radius;
    vec3 center;
    vec3 normal;
    Material* mat;
} Disk;

float rayIntersectDisk(ShadeRec* sr, Disk* disk, const Ray ray);
float shadowRayIntersectDisk(const Disk* disk, const Ray ray);
