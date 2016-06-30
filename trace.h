#pragma once

#include "sampling.h"
#include "util/vec.h"
#include "util/ray.h"
#include "util/shaderec.h"
#include "util/constants.h"
#include "sceneobj.h"
#include "shading.h"
#include "intersect.h"

enum TraceType
{
    RAYCAST,
    WHITTED,
    PATHTRACE
};

void raycast(vec3, int, const vec3, const Ray, const SceneObjects*, const SceneLights*);
void whittedTrace(vec3, int, const vec3, const Ray, const SceneObjects*, const SceneLights*);
void pathTrace(vec3, int, const vec3, const Ray, const SceneObjects*, const SceneLights*);

typedef void (*traceFunc)(vec3, int, const vec3, const Ray, const SceneObjects*, const SceneLights*);

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

void raycast(vec3 radiance, int depth, const vec3 sample, const Ray ray,
             const SceneObjects* so, const SceneLights* sl)
{
    vec3_copy(radiance, ORIGIN);
    float min_t = TMAX;
    ShadeRec min_sr;
    min_t = intersectTest(&min_sr, so, ray);
            
    // Shading
    if(min_t < TMAX)
    {
        if(min_sr.mat->mat_type == EMISSIVE)
        {
            vec3_scale(radiance, min_sr.mat->ce, min_sr.mat->ke/1.0f);
            maxToOne(radiance, radiance);
        }else
        {
            // Add ambient component to radiance
            ambientShading(radiance, sl->amb_light, sample, so, &min_sr);

            // Direct illumination
            for(int i = 0; i < sl->num_lights; i++)
            {
                vec3 light_dir;
                getLightDir(light_dir, sl->light_types[i], sl->light_ptrs[i], &min_sr);
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
}

// Calculate specular incident radiance and brdf
void calcSpecRefRadiance(vec3 spec_ref_radiance, const int depth,  const vec3 h_sample,
                         const Ray ray, const ShadeRec* sr, const SceneObjects* so, const SceneLights*  sl)
{
    vec3 reflect_dir;
    calcReflectRayDir(reflect_dir, sr->normal, ray.direction);
    vec3 new_sample;
    getNextSample3D(new_sample, sr->mat->h_samples);
    Ray sample_ray;
    vec3_copy(sample_ray.origin, sr->hit_point);
    getVec3InLocalBasis(sample_ray.direction, new_sample, reflect_dir);
    if(vec3_dot(sample_ray.direction, sr->normal) < 0)
    {
        vec3 reflected_sample = {-new_sample[0], -new_sample[1], new_sample[2]};
        getVec3InLocalBasis(sample_ray.direction, reflected_sample, reflect_dir);
    }
    double phong_lobe = pow(vec3_dot(sample_ray.direction, reflect_dir), sr->mat->exp);
    float pdf = phong_lobe * (vec3_dot(sr->normal, sample_ray.direction));
    vec3 brdf;
    vec3_scale(brdf, sr->mat->cs, phong_lobe);   // Deferring reflection constant scaling
    whittedTrace(spec_ref_radiance, depth-1, h_sample, sample_ray, so, sl);
    vec3_mult(spec_ref_radiance, spec_ref_radiance, brdf);
    vec3_scale(spec_ref_radiance, spec_ref_radiance, vec3_dot(sr->normal, sample_ray.direction) / pdf);
}

bool totalInternalReflection(const ShadeRec* sr)
{
    float cos_theta_i = vec3_dot(sr->normal, sr->wo);
    float eta = sr->mat->ior;

    if(cos_theta_i < 0.0f)
    {
        eta = 1.0f / eta;
    }
    return (1.0f - (1.0f - cos_theta_i * cos_theta_i) / (eta * eta) < 0.0f);
}

float calcTransmitDir(vec3 transmit_dir, const vec3 normal, const vec3 wo, const float ior)
{
    vec3 n;
    vec3_copy(n, normal);
    float cos_theta_i = vec3_dot(n, wo);
    float eta = ior;
    if(cos_theta_i < 0.0f)
    {
        cos_theta_i = -cos_theta_i;
        vec3_negate(n, n);
        eta = 1.0f / eta;                    
    }

    float tmp = 1.0f - (1.0f - cos_theta_i * cos_theta_i) / (eta * eta);
    float cos_theta2 = (float)sqrt(tmp);
    vec3 tmp_vec3_1, tmp_vec3_2;
    vec3_scale(tmp_vec3_1, wo, -1.0f / eta);
    vec3_scale(tmp_vec3_2, n, cos_theta2 - cos_theta_i / eta);
    vec3_sub(transmit_dir, tmp_vec3_1, tmp_vec3_2);
    return eta;
}

void whittedTrace(vec3 radiance, int depth, const vec3 h_sample,
                     const Ray ray, const SceneObjects* so, const SceneLights* sl)
{
    vec3_copy(radiance, ORIGIN);
    float min_t = TMAX;
    ShadeRec min_sr;
    min_t = intersectTest(&min_sr, so, ray);            
    // Shading
    if(min_t < TMAX)
    {
        if(min_sr.mat->mat_type == EMISSIVE)
        {
            vec3_scale(radiance, min_sr.mat->ce, min_sr.mat->ke/1.0f);
            maxToOne(radiance, radiance);
        }else
        {
            ambientShading(radiance, sl->amb_light, h_sample, so, &min_sr);
            // Direct illumination
            for(int i = 0; i < sl->num_lights; i++)
            {
                vec3 light_dir;
                getLightDir(light_dir, sl->light_types[i], sl->light_ptrs[i], &min_sr);
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
                vec3 reflected_illum = {0.0f, 0.0f, 0.0f};
                calcSpecRefRadiance(reflected_illum, depth, h_sample, ray, &min_sr, so, sl);

                if(min_sr.mat->mat_type == REFLECTIVE)
                {
                    vec3_add(radiance, radiance, reflected_illum);
                    return;
                }

                // Transmitted radiance
                vec3 transmit_dir = {0, 0, 0};
                float eta = calcTransmitDir(transmit_dir, min_sr.normal, min_sr.wo, min_sr.mat->ior);
                float ndotwt = fabs(vec3_dot(min_sr.normal, transmit_dir));
                if(!totalInternalReflection(&min_sr))
                {
                    vec3 btdf;
                    vec3_scale(btdf, WHITE, min_sr.mat->kt / (eta*eta) / ndotwt);
                    Ray transmitted_ray;
                    vec3_copy(transmitted_ray.origin, min_sr.hit_point);
                    vec3_copy(transmitted_ray.direction, transmit_dir);

                    vec3 transmitted_illum = {0.0f, 0.0f, 0.0f};
                    whittedTrace(transmitted_illum, depth-1, h_sample, transmitted_ray, so, sl);
                    vec3_scale(transmitted_illum, transmitted_illum, ndotwt);
                    vec3_mult(transmitted_illum, transmitted_illum, btdf);
                    vec3_add(radiance, radiance, transmitted_illum);
                    // Scaling reflection since there's no total internal reflection
                    vec3_scale(reflected_illum, reflected_illum, min_sr.mat->kr);
                }
                vec3_add(radiance, radiance, reflected_illum);                    
            }
        }        
    }else
    {
        vec3_copy(radiance, sl->bg_color);
    }
}

void calcSpecRadiancePT(vec4 ref_radiance, const Ray ray, const ShadeRec* sr, const vec3 h_sample,
                           const int depth, const SceneObjects* so, const SceneLights* sl)
{
    vec3 reflect_dir;
    calcReflectRayDir(reflect_dir, sr->normal, ray.direction);
    vec3 new_sample;
    getNextSample3D(new_sample, sr->mat->h_samples);
    Ray sample_ray;
    vec3_copy(sample_ray.origin, sr->hit_point);
    getVec3InLocalBasis(sample_ray.direction, new_sample, reflect_dir);
    if(vec3_dot(sample_ray.direction, sr->normal) < 0)
    {
        vec3 reflected_sample = {-new_sample[0], -new_sample[1], new_sample[2]};
        getVec3InLocalBasis(sample_ray.direction, reflected_sample, reflect_dir);
    }
    double phong_lobe = pow(vec3_dot(sample_ray.direction, reflect_dir), sr->mat->exp);
    float pdf = phong_lobe * (vec3_dot(sr->normal, sample_ray.direction));
    vec3 brdf;
    vec3_scale(brdf, sr->mat->cs, phong_lobe); // Defer reflection scaling
    pathTrace(ref_radiance, depth-1, h_sample, sample_ray, so, sl);                
    vec3_scale(ref_radiance, ref_radiance, vec3_dot(sr->normal, sample_ray.direction) / pdf);
    vec3_mult(ref_radiance, ref_radiance, brdf);
}

void pathTrace(vec3 radiance, int depth, const vec3 h_sample, const Ray ray, const SceneObjects* so,
    const SceneLights* sl)
{
    vec3_copy(radiance, ORIGIN);
    float min_t = TMAX;
    ShadeRec min_sr;
    min_t = intersectTest(&min_sr, so, ray);
    // Shading
    if(min_t < TMAX)
    {    
        if(min_sr.mat->mat_type == EMISSIVE)
        {
            vec3_scale(radiance, min_sr.mat->ce, min_sr.mat->ke);
            maxToOne(radiance, radiance);                
        }else
        {
            if(depth > 0)
            {
                if(min_sr.mat->mat_type == MATTE)
                {
                    vec3 new_sample;
                    getNextSample3D(new_sample, min_sr.mat->h_samples);
                    Ray sample_ray;
                    getVec3InLocalBasis(sample_ray.direction, new_sample, min_sr.normal);
                    vec3_copy(sample_ray.origin, min_sr.hit_point);

                    vec3 inc_radiance;
                    pathTrace(inc_radiance, depth-1, h_sample, sample_ray, so, sl);
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
                    calcSpecRadiancePT(reflected_illum, ray, &min_sr, h_sample, depth, so, sl);
                    vec3_scale(reflected_illum, reflected_illum, min_sr.mat->ks);
                    vec3_add(radiance, radiance, reflected_illum);
                }

                if(min_sr.mat->mat_type == TRANSPARENT)
                {
                    vec3 reflected_illum;
                    calcSpecRadiancePT(reflected_illum, ray, &min_sr, h_sample, depth, so, sl);

                    vec3 transmit_dir = {0, 0, 0};
                    float eta = calcTransmitDir(transmit_dir, min_sr.normal, min_sr.wo, min_sr.mat->ior);
                    float ndotwt = fabs(vec3_dot(min_sr.normal, transmit_dir));
                    if(!totalInternalReflection(&min_sr))
                    {
                        vec3 btdf;
                        vec3_scale(btdf, WHITE, min_sr.mat->kt / (eta*eta) / ndotwt);
                        Ray transmitted_ray;
                        vec3_copy(transmitted_ray.origin, min_sr.hit_point);
                        vec3_copy(transmitted_ray.direction, transmit_dir);
                        
                        vec3 transmitted_illum = {0.0f, 0.0f, 0.0f};
                        pathTrace(transmitted_illum, depth-1, h_sample, transmitted_ray, so, sl);
                        vec3_scale(transmitted_illum, transmitted_illum, ndotwt);
                        vec3_mult(transmitted_illum, transmitted_illum, btdf);
                        vec3_add(radiance, radiance, transmitted_illum);
                        // Scaling reflection since there's no total internal reflection
                        vec3_scale(reflected_illum, reflected_illum, min_sr.mat->kr);
                    }                    
                    vec3_add(radiance, radiance, reflected_illum);                    
                }
            }
        }
    }else
    {
        vec3_copy(radiance, sl->bg_color);
    }
}
