#pragma once

#include "util/constants.h"
#include "util/vec.h"
#include "util/util.h"
#include "materials.h"
#include "shapes/shapes.h"
#include "lights.h"
#include "scene/scenedata.h"
#include "sampling.h"
#include "intersect.h"
#include "texture.h"
#include "noise.h"


float AOTest(const vec3 h_sample, const SceneObjects *so, const ShadeRec* sr)
{
    //vec3 h_sample;
    //getNextSample3D(h_sample, samples);
    Ray shadow_ray;
    getVec3InLocalBasis(shadow_ray.direction, h_sample, sr->normal);
    vec3_copy(shadow_ray.origin, sr->hit_point);
    return shadowIntersectTest(so, shadow_ray, TMAX);
}

//extern void maxToOne(vec3, const vec3);
inline float GammaCorrect(float value) {
    if (value <= 0.0031308f) return 12.92f * value;
    return 1.055f * powf(value, (1.f / 2.4f)) - 0.055f;
}

// divide vec3 a by its max component if max component > 1
void toneMap(vec3 r, const vec3 a)
{
    /*
    float max;
    max = (a[0] > a[1]) ? a[0] : a[1];
    max = (max > a[2]) ? max : a[2];
    if(max > 1.0f)
    {
        vec3_scale(r, a, 1.0f/max);
    }
    */

    float exposure = -2.0f;
    r[0] = 1.0f - expf(a[0] * exposure);
    r[1] = 1.0f - expf(a[1] * exposure);
    r[2] = 1.0f - expf(a[2] * exposure);
    vec3_pow(r, r, 1.0f / 2.2f);
    /*
    vec3 denom;
    denom[0] = 1.0f / (a[0] + 1.0f);
    denom[1] = 1.0f / (a[1] + 1.0f);
    denom[2] = 1.0f / (a[2] + 1.0f);    
    vec3_mult(r, a, denom);
    vec3_pow(r, r, 1.0f / 2.2f);
    */
    /*
    r[0] = clamp(GammaCorrect(a[0]) + 0.5f, 0.f, 1.f);
    r[1] = clamp(GammaCorrect(a[1]) + 0.5f, 0.f, 1.f);
    r[2] = clamp(GammaCorrect(a[2]) + 0.5f, 0.f, 1.f);
    */
}
