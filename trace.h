#pragma once

#include "sampling.h"
#include "util/vec.h"
#include "util/ray.h"
#include "util/shaderec.h"
#include "util/constants.h"
#include "scene/scenedata.h"
#include "shading.h"
#include "intersect.h"

extern int MAX_DEPTH;
#define SEPARATE_DIRECT_INDIRECT

enum TraceType
{
    RAYCAST,
    WHITTED,
    PATHTRACE
};

typedef struct TraceArgs_s
{
    // Data that are constant throughout a single sample
    int sample_index;
    SceneObjects *objects;
    SceneLights *lights;
    vec3 h_sample;

    // Data that may get modified during the tracing of a sample
    Material* medium_mat;
}TraceArgs;

/*
  float raycast(vec3, int, const vec3, const Ray,
  const SceneObjects*, const SceneLights*, const int, const unsignedint);
  float whittedTrace(vec3, int, const vec3, const Ray,
  const SceneObjects*, const SceneLights*, const int, const unsigned int);
  float pathTrace(vec3, int, const vec3, const Ray,
  const SceneObjects*, const SceneLights*, const int, cosnt unsigned int);
*/

float raycast(vec3, int, const Ray, TraceArgs trace_args);
float whittedTrace(vec3, int, const Ray, TraceArgs trace_args);
float pathTrace(vec3, int, const Ray, TraceArgs trace_args);
float raycastMedium(vec3, int, const Ray, TraceArgs trace_args);
float whittedTraceMedium(vec3, int, const Ray, TraceArgs trace_args);

/*
typedef float (*traceFunc)(vec3, int, const vec3, const Ray,
 const SceneObjects*, const SceneLights*, const int, const unsigned int);
*/
typedef float (*traceFunc)(vec3, int, const Ray, TraceArgs trace_args);

traceFunc getTraceFunc(const TraceType trace_type)
{
    traceFunc func;
    switch(trace_type)
    {
    case RAYCAST:
        func = &raycast;
        break;
    case WHITTED:
        func = &whittedTrace;
        break;
    case PATHTRACE:
        func = &pathTrace;
        break;
    default:
        func = &raycast;
    }
    return func;
}

float raycast(vec3 radiance, int depth, const Ray ray, TraceArgs trace_args)
{
    const SceneObjects *so = trace_args.objects;
    const SceneLights *sl = trace_args.lights;
    const int sample_index = trace_args.sample_index;
    vec3 h_sample;
    vec3_copy(h_sample, trace_args.h_sample);

    vec3_copy(radiance, BLACK);
    ShadeRec min_sr;
    float min_t = intersectTest(&min_sr, so, ray);

    // Shading
    if(min_t < TMAX)
    {
        if(min_sr.mat->mat_type == EMISSIVE)
        {
            vec3_scale(radiance, min_sr.mat->ce, min_sr.mat->ke/1.0f);
            maxToOne(radiance, radiance);
        }else if(min_sr.mat->mat_type == PARTICIPATING)
        {
            TraceArgs new_trace_args = trace_args;
            new_trace_args.medium_mat = min_sr.mat;
            Ray new_ray;
            getPointOnRay(new_ray.origin, ray, min_t);
            vec3_copy(new_ray.direction, ray.direction);
            raycastMedium(radiance, depth-1, new_ray, new_trace_args);
        }else
        {
            // Add ambient component to radiance
            ambientShading(radiance, sl->amb_light, h_sample, so, &min_sr);

            // Direct illumination
            for(int i = 0; i < sl->num_lights; i++)
            {
                vec3 light_dir;
                getLightDir(light_dir, sl->light_types[i], sl->light_ptrs[i], &min_sr, sample_index);
                float ndotwi = vec3_dot(light_dir, min_sr.normal);
                if(ndotwi > 0)
                {
                    bool in_shadow = shadowTest(i, sl, so, light_dir, &min_sr);
                    if(!in_shadow)
                    {
                        directIllumShading(radiance, ndotwi, light_dir, sl->light_ptrs[i],
                                           sl->light_types[i], &min_sr);
                    }
                }
            }
        }
    }else
    {
        vec3_copy(radiance, sl->bg_color);
    }
    return min_t;
}

// Calculate specular incident radiance and brdf
float calcSpecRefRadiance(vec3 spec_ref_radiance,
                          const int depth,  const Ray ray, const ShadeRec* sr,
                          TraceArgs trace_args)
{
    const int sample_index = trace_args.sample_index;
    vec3 reflect_dir, normal;
    vec3_copy(normal, sr->normal);
    calcReflectRayDir(reflect_dir, normal, ray.direction);
    if(vec3_dot(normal, reflect_dir) < 0.0f && sr->mat->mat_type == TRANSPARENT)
    {
        vec3_negate(normal, normal);
    }
    vec3 new_sample;
    getSample3D(new_sample, sr->mat->h_samples, sample_index);
    Ray sample_ray;
    vec3_copy(sample_ray.origin, sr->hit_point);
    getVec3InLocalBasis(sample_ray.direction, new_sample, reflect_dir);
    if(vec3_dot(sample_ray.direction, normal) < 0.0f)
    {
        vec3 reflected_sample = {-new_sample[0], -new_sample[1], new_sample[2]};
        getVec3InLocalBasis(sample_ray.direction, reflected_sample, reflect_dir);
    }
    float raydotr = vec3_dot(sample_ray.direction, reflect_dir);
    float phong_lobe = pow(vec3_dot(sample_ray.direction, reflect_dir), sr->mat->exp);
    float pdf = phong_lobe * (vec3_dot(normal, sample_ray.direction));
    vec3 brdf;
    vec3_scale(brdf, sr->mat->cs, phong_lobe);   // Deferring reflectance scaling
    float t = whittedTrace(spec_ref_radiance, depth-1, sample_ray, trace_args);
    vec3_mult(spec_ref_radiance, spec_ref_radiance, brdf);
    vec3_scale(spec_ref_radiance, spec_ref_radiance, vec3_dot(normal, sample_ray.direction) / pdf);
    return t;
}

bool totalInternalReflection(const ShadeRec* sr)
{
    float cos_theta_i = vec3_dot(sr->normal, sr->wo);
    float eta = sr->mat->ior_in / sr->mat->ior_out;

    if(cos_theta_i < 0.0f)
    {
        eta = 1.0f / eta;
    }
    return (1.0f - (1.0f - cos_theta_i * cos_theta_i) / (eta * eta) < 0.0f);
}

float calcTransmitDir(vec3 transmit_dir, const ShadeRec* sr)
{
    vec3 n;
    vec3_copy(n, sr->normal);
    float cos_theta_i = vec3_dot(n, sr->wo);
    float eta = sr->mat->ior_in / sr->mat->ior_out;    
    if(cos_theta_i < 0.0f)
    {
        cos_theta_i = -cos_theta_i;
        vec3_negate(n, n);
        eta = 1.0f / eta;                    
    }

    float tmp = 1.0f - (1.0f - cos_theta_i * cos_theta_i) / (eta * eta);
    float cos_theta2 = (float)sqrt(tmp);
    vec3 tmp_vec3_1, tmp_vec3_2;
    vec3_scale(tmp_vec3_1, sr->wo, -1.0f / eta);
    vec3_scale(tmp_vec3_2, n, cos_theta2 - cos_theta_i / eta);
    vec3_sub(transmit_dir, tmp_vec3_1, tmp_vec3_2);
    return eta;
}

float calcFresnelReflectance(const ShadeRec* sr)
{
    vec3 n;
    vec3_copy(n, sr->normal);
    // NOTE: not sure about the next line
    float cos_theta_i = vec3_dot(n, sr->wo);
    float eta = sr->mat->ior_in / sr->mat->ior_out;    
    if(cos_theta_i < 0.0f)
    {
        cos_theta_i = -cos_theta_i;
        vec3_negate(n, n);
        eta = 1.0f / eta;
    }
    float tmp = 1.0f - (1.0f - cos_theta_i * cos_theta_i) / (eta * eta);
    if(tmp < 0.0f)
    {
        return 1.0f;
    }
    float cos_theta_t = (float)sqrt(tmp);
    float r_parallel = (eta * cos_theta_i - cos_theta_t) / (eta * cos_theta_i + cos_theta_t);
    float r_perpendicular = (cos_theta_i - eta * cos_theta_t) / (cos_theta_i + eta * cos_theta_t);
    float fresnel_reflectance = 0.5f * (r_parallel * r_parallel + r_perpendicular * r_perpendicular);
    return fresnel_reflectance;
}

float whittedTrace(vec3 radiance, int depth, const Ray ray, TraceArgs trace_args)
{
    const SceneObjects *so = trace_args.objects;
    const SceneLights *sl = trace_args.lights;
    const int sample_index = trace_args.sample_index;
    vec3 h_sample;
    vec3_copy(h_sample, trace_args.h_sample);

    // debug variables
    bool intersected = false, indirect = false, not_tir = false, tir = false;

    float reflect_t = 0, transmit_t = 0;
    float kr = 0.0f, kt = 0.0f;
    vec3 reflected_illum = {0.0f, 0.0f, 0.0f};
    vec3 transmitted_illum = {0.0f, 0.0f, 0.0f};
    vec3_copy(radiance, ORIGIN);

    float min_t = TMAX;
    ShadeRec min_sr;
    min_t = intersectTest(&min_sr, so, ray);
    // Shading
    if(min_t < TMAX)
    {
        intersected = true;
        if(min_sr.mat->mat_type == EMISSIVE)
        {
            vec3_scale(radiance, min_sr.mat->ce, min_sr.mat->ke/1.0f);
            //maxToOne(radiance, radiance);
        }else
        {
            ambientShading(radiance, sl->amb_light, h_sample, so, &min_sr);
            // Direct illumination

            for(int i = 0; i < sl->num_lights; i++)
            {
                vec3 light_dir;
                getLightDir(light_dir, sl->light_types[i], sl->light_ptrs[i], &min_sr, sample_index);
                float ndotwi = vec3_dot(light_dir, min_sr.normal);
                if(ndotwi > 0)
                {
                    bool in_shadow = shadowTest(i, sl, so, light_dir, &min_sr);
                    if(!in_shadow)
                    {
                        directIllumShading(radiance, ndotwi, light_dir, sl->light_ptrs[i],
                                           sl->light_types[i], &min_sr);
                    }
                }
            }

            // Indirect illumination
            if(depth > 0 && (min_sr.mat->mat_type == REFLECTIVE || min_sr.mat->mat_type == TRANSPARENT))
            {
                indirect = true;
                //float reflect_t;
                //vec3 reflected_illum = {0.0f, 0.0f, 0.0f};
                reflect_t = calcSpecRefRadiance(reflected_illum, depth, ray, &min_sr, trace_args);

                if(min_sr.mat->mat_type == REFLECTIVE)
                {
                    vec3_scale(reflected_illum, reflected_illum, min_sr.mat->ks);
                    vec3_add(radiance, radiance, reflected_illum);
                    return min_t;
                }
                float ndotwo = vec3_dot(min_sr.normal, min_sr.wo);
                if(!totalInternalReflection(&min_sr))
                {
                    not_tir = true;
                    // Transmitted radiance
                    vec3 transmit_dir = {0, 0, 0};
                    float eta = calcTransmitDir(transmit_dir, &min_sr);
                    float ndotwt = fabs(vec3_dot(min_sr.normal, transmit_dir));
                    //float kr = calcFresnelReflectance(&min_sr);
                    //float kt = 1.0f - kr;
                    kr = calcFresnelReflectance(&min_sr);
                    kt = 1.0f - kr;

                    vec3 btdf;
                    vec3_scale(btdf, WHITE, kt / (eta*eta) / ndotwt);
                    Ray transmitted_ray;
                    vec3_copy(transmitted_ray.origin, min_sr.hit_point);
                    vec3_copy(transmitted_ray.direction, transmit_dir);

                    //float transmit_t;
                    //vec3 transmitted_illum = {0.0f, 0.0f, 0.0f};
                    transmit_t = whittedTrace(transmitted_illum, depth-1, transmitted_ray, trace_args);
                    vec3_scale(transmitted_illum, transmitted_illum, ndotwt);
                    vec3_mult(transmitted_illum, transmitted_illum, btdf);
                    // Scaling reflection since there's no total internal reflection
                    vec3_scale(reflected_illum, reflected_illum, kr);

                    vec3 color_filter_ref, color_filter_trans;
                    if(ndotwo > 0.0f)
                    {
                        vec3_pow(color_filter_ref, min_sr.mat->cf_out, reflect_t);
                        vec3_pow(color_filter_trans, min_sr.mat->cf_in, transmit_t);
                    }else
                    {
                        vec3_pow(color_filter_ref, min_sr.mat->cf_in, reflect_t);
                        vec3_pow(color_filter_trans, min_sr.mat->cf_out, transmit_t);
                    }
                    vec3_mult(reflected_illum, reflected_illum, color_filter_ref);
                    vec3_mult(transmitted_illum, transmitted_illum, color_filter_trans);
                    vec3_add(radiance, radiance, transmitted_illum);
                }else
                {

                    tir = true;
                    vec3 color_filter;
                    if(ndotwo > 0.0f)
                    {
                        vec3_pow(color_filter, min_sr.mat->cf_out, reflect_t);
                    }else
                    {
                        vec3_pow(color_filter, min_sr.mat->cf_in, reflect_t);
                    }
                    vec3_mult(reflected_illum, reflected_illum, color_filter);
                }
                vec3_add(radiance, radiance, reflected_illum);
            }
        }
    }else
    {
        vec3_copy(radiance, sl->bg_color);
    }
    //vec3_copy(radiance, min_sr.normal);
    /*
    float tmpf = (min_t - 800.0f) / 600.0f;
    vec3 depth_f = {tmpf, tmpf, tmpf};
    vec3_copy(radiance, depth_f);
    */
    return min_t;
}

float calcSpecRadiancePT(vec4 ref_radiance, const Ray ray, const ShadeRec* sr,
                         const int depth, TraceArgs trace_args)
{
    const int sample_index = trace_args.sample_index;
    vec3 reflect_dir;
    calcReflectRayDir(reflect_dir, sr->normal, ray.direction);
    vec3 new_sample;
    getSample3D(new_sample, sr->mat->h_samples, sample_index + (MAX_DEPTH - depth));
    Ray sample_ray;
    vec3_copy(sample_ray.origin, sr->hit_point);
    getVec3InLocalBasis(sample_ray.direction, new_sample, reflect_dir);
    if(vec3_dot(sample_ray.direction, sr->normal) < 0)
    {
        vec3 reflected_sample = {-new_sample[0], -new_sample[1], new_sample[2]};
        getVec3InLocalBasis(sample_ray.direction, reflected_sample, reflect_dir);
    }
    float phong_lobe = pow(vec3_dot(sample_ray.direction, reflect_dir), sr->mat->exp);
    float pdf = phong_lobe * (vec3_dot(sr->normal, sample_ray.direction));
    vec3 brdf;
    vec3_scale(brdf, sr->mat->cs, phong_lobe); // Defer reflectance scaling
    float t;
    t = pathTrace(ref_radiance, depth-1, sample_ray, trace_args);
    vec3_scale(ref_radiance, ref_radiance, vec3_dot(sr->normal, sample_ray.direction) / pdf);
    vec3_mult(ref_radiance, ref_radiance, brdf);
    return t;
}

float pathTrace(vec3 radiance, int depth, const Ray ray, TraceArgs trace_args)
{
    const SceneObjects *so = trace_args.objects;
    const SceneLights *sl = trace_args.lights;
    const int sample_index = trace_args.sample_index;
    vec3 h_sample;
    vec3_copy(h_sample, trace_args.h_sample);

    vec3_copy(radiance, ORIGIN);
    float min_t = TMAX;
    ShadeRec min_sr;
    min_t = intersectTest(&min_sr, so, ray);
    // Shading
    if(min_t < TMAX)
    {
        if(min_sr.mat->mat_type == EMISSIVE)
        {
#ifdef SEPARATE_DIRECT_INDIRECT
            if(depth == MAX_DEPTH - 1)
            {
                vec3_copy(radiance, BLACK);
            }else
            {
                vec3_scale(radiance, min_sr.mat->ce, min_sr.mat->ke);
            }
#else
            vec3_scale(radiance, min_sr.mat->ce, min_sr.mat->ke);
#endif
            //maxToOne(radiance, radiance);
        }else
        {
            if(depth > 0)
            {
#ifdef SEPARATE_DIRECT_INDIRECT
                // Direct illumination
                if(min_sr.mat->mat_type & (MATTE | PHONG) && depth == MAX_DEPTH)
                {
                    for(int i = 0; i < sl->num_lights; i++)
                    {
                        if(sl->light_types[i] == AREALIGHT)
                        {
                            vec3 light_dir;
                            getLightDir(light_dir, sl->light_types[i], sl->light_ptrs[i], &min_sr, sample_index);
                            float ndotwi = vec3_dot(light_dir, min_sr.normal);
                            if(ndotwi > 0)
                            {
                                bool in_shadow = shadowTest(i, sl, so, light_dir, &min_sr);
                                if(!in_shadow)
                                {
                                    directIllumShading(radiance, ndotwi, light_dir, sl->light_ptrs[i],
                                                       sl->light_types[i], &min_sr);
                                }
                            }
                        }
                    }
                }
#endif
                if(min_sr.mat->mat_type == MATTE)
                {
                    vec3 new_sample;
                    getSample3D(new_sample, min_sr.mat->h_samples, sample_index + (MAX_DEPTH - depth));
                    Ray sample_ray;
                    getVec3InLocalBasis(sample_ray.direction, new_sample, min_sr.normal);
                    vec3_copy(sample_ray.origin, min_sr.hit_point);

                    vec3 inc_radiance;
                    pathTrace(inc_radiance, depth-1, sample_ray, trace_args);
                    vec3 brdf;
                    vec3_scale(brdf, min_sr.mat->cd, min_sr.mat->kd / (float)PI);
                    float ndotwi = vec3_dot(min_sr.normal, sample_ray.direction);
                    float pdf = ndotwi / (float)PI;

                    vec3 tmp;
                    vec3_mult(tmp, inc_radiance, brdf);
                    vec3_scale(tmp, tmp, ndotwi / pdf);
                    vec3_add(radiance, radiance, tmp);
                }

                if(min_sr.mat->mat_type == REFLECTIVE || min_sr.mat->mat_type == PHONG)
                {
                    vec3 reflected_illum;
                    // NOTE: depth - 1 is a shitty hack
                    calcSpecRadiancePT(reflected_illum, ray, &min_sr, depth-1, trace_args);
                    vec3_scale(reflected_illum, reflected_illum, min_sr.mat->ks);
                    vec3_add(radiance, radiance, reflected_illum);
                }

                if(min_sr.mat->mat_type == TRANSPARENT)
                {
                    float reflect_t = TMAX;
                    vec3 reflected_illum = {0.0f, 0.0f, 0.0f};
                    float ndotwo = vec3_dot(min_sr.normal, min_sr.wo);
                    float kr = calcFresnelReflectance(&min_sr);
                    //if(!totalInternalReflection(&min_sr))

                    float transmit_t = TMAX;
                    vec3 transmitted_illum = {0.0f, 0.0f, 0.0f};

                    float rand_float = (float)rand() / (float)RAND_MAX;
                    if(rand_float <= kr) // Reflection
                    {
                        reflect_t = calcSpecRadiancePT(reflected_illum, ray, &min_sr, depth, trace_args);
                    }else // Transmission
                    {
                        vec3 transmit_dir = {0, 0, 0};
                        float eta = calcTransmitDir(transmit_dir, &min_sr);
                        float ndotwt = fabs(vec3_dot(min_sr.normal, transmit_dir));
                        float kt = 1.0f - kr;

                        vec3 btdf;
                        vec3_scale(btdf, WHITE, kt / (eta*eta) / ndotwt);
                        Ray transmitted_ray;
                        vec3_copy(transmitted_ray.origin, min_sr.hit_point);
                        vec3_copy(transmitted_ray.direction, transmit_dir);

                        transmit_t = pathTrace(transmitted_illum, depth-1, transmitted_ray, trace_args);
                        vec3_scale(transmitted_illum, transmitted_illum, ndotwt);
                        vec3_mult(transmitted_illum, transmitted_illum, btdf);
                    }

                    vec3 color_filter_ref, color_filter_trans;
                    if(ndotwo > 0.0f)
                    {
                        vec3_pow(color_filter_ref, min_sr.mat->cf_out, reflect_t);
                        vec3_pow(color_filter_trans, min_sr.mat->cf_in, transmit_t);
                    }else
                    {
                        vec3_pow(color_filter_ref, min_sr.mat->cf_in, reflect_t);
                        vec3_pow(color_filter_trans, min_sr.mat->cf_out, transmit_t);
                    }
                    vec3_mult(reflected_illum, reflected_illum, color_filter_ref);
                    vec3_mult(transmitted_illum, transmitted_illum, color_filter_trans);
                    vec3_add(radiance, radiance, transmitted_illum);
                    vec3_add(radiance, radiance, reflected_illum);
                }
            }else
            {
                vec3_copy(radiance, sl->bg_color);
            }
        }
    }else
    {
        vec3_copy(radiance, sl->bg_color);
    }
    return min_t;
}

