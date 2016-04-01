#pragma once

#include "constants.h"
#include "vec.h"
#include "materials.h"
#include "shapes.h"
#include "lights.h"
#include "sceneobj.h"
#include "sampling.h"
 
void diffuseShading(vec3 radiance, const float ndotwi, const vec3 inc_radiance_cos,
                    const ShadeRec* sr)
{
    vec3 reflectance;
    vec3_scale(reflectance, sr->mat->cd, sr->mat->kd);
    vec3 f;
    vec3_scale(f, reflectance, 1.0f/(float)PI);
    vec3 tmp;
    vec3_mult(tmp, inc_radiance_cos, f);
    vec3_add(radiance, radiance, tmp);
}

void specularShading(vec3 radiance, const vec3 wo, const vec3 light_dir,
                     const vec3 inc_radiance_cos, const ShadeRec* sr)
{
    // ks*cs * cos(theta)^exp * inc_radiance_cos
    // theta = angle between light and reflected view vector
    // reflect vector: 2*dot(n, v)*n - v    
    float ndotwo = vec3_dot(sr->normal, wo);
    vec3 tmp, wo_neg, reflect_dir;
    vec3_negate(wo_neg, wo);
    vec3_scale(tmp, sr->normal, ndotwo * 2);
    vec3_add(reflect_dir, wo_neg, tmp);    
    vec3_normalize(reflect_dir, reflect_dir);
    float rdotwi = vec3_dot(reflect_dir, light_dir);
    vec3_scale(tmp, sr->mat->cs, sr->mat->ks * pow(rdotwi, sr->mat->exp));
    vec3_mult(tmp, inc_radiance_cos, tmp);
    vec3_add(radiance, radiance, tmp);
}
/*
float AOTest(Samples* samples, const SceneObjects *so, const ShadeRec* sr)
{
    vec3 h_sample;
    getNextHemisphereSample(h_sample, samples);
    vec3 u, v, w;
    vec3_copy(w, sr->normal);
    vec3_cross(v, w, JITTERED_UP);
    vec3_normalize(v, v);
    vec3_cross(u, v, w);
    Ray shadow_ray;
    orthoNormalTransform(shadow_ray.direction, u, v, w, h_sample);
    vec3_copy(shadow_ray.origin, sr->hit_point);
    return shadowIntersectTest(so, shadow_ray);
}
*/
float AOTest(Samples3D* samples, const SceneObjects *so, const ShadeRec* sr)
{
    vec3 h_sample;
    getNextSample3D(h_sample, samples);
    vec3 u, v, w;
    vec3_copy(w, sr->normal);
    vec3_cross(v, w, JITTERED_UP);
    vec3_normalize(v, v);
    vec3_cross(u, v, w);
    Ray shadow_ray;
    orthoNormalTransform(shadow_ray.direction, u, v, w, h_sample);
    vec3_copy(shadow_ray.origin, sr->hit_point);
    return shadowIntersectTest(so, shadow_ray);
}

