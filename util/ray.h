#pragma once

#include "vec.h"

typedef struct Ray
{
    vec3 origin;
    vec3 direction;
} Ray;

void getPointOnRay(vec3 r, const Ray ray, const float t);
void calcReflectRayDir(vec3 reflect_dir, const vec3 normal, const vec3 incident_dir);
float calcTransmitDir(vec3 transmit_dir, const vec3 normal, const vec3 wo, const float ior_in,
                      const float ior_out);
