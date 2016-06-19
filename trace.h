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

void raycast(vec3 radiance, int depth, const vec3 sample, const Ray ray, const SceneObjects* so, const SceneLights* sl)
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
            if(depth > 0 && min_sr.mat->mat_type == PHONG)
            {
                vec3 reflect_dir;
                calcReflectRayDir(reflect_dir, min_sr.normal, ray.direction);

                vec3 new_sample;
                getNextSample3D(new_sample, min_sr.mat->h_samples);
                Ray sample_ray;
                vec3_copy(sample_ray.origin, min_sr.hit_point);
                getVec3InLocalBasis(sample_ray.direction, new_sample, reflect_dir);
                if(vec3_dot(sample_ray.direction, min_sr.normal) < 0)
                {
                    vec3 reflected_sample = {-new_sample[0], -new_sample[1], new_sample[2]};
                    getVec3InLocalBasis(sample_ray.direction, reflected_sample, reflect_dir);
                }
                double phong_lobe = pow(vec3_dot(sample_ray.direction, reflect_dir), min_sr.mat->exp);

                float pdf = phong_lobe * (vec3_dot(min_sr.normal, sample_ray.direction));

                vec3 reflectance;
                vec3_scale(reflectance, min_sr.mat->cs, min_sr.mat->ks);                
                vec3 fr;
                vec3_scale(fr, reflectance, phong_lobe);

                vec3 indirect_illum = {0.0f, 0.0f, 0.0f};
                whittedTrace(indirect_illum, depth-1, h_sample, sample_ray, so, sl);
                
                vec3_scale(indirect_illum, indirect_illum, vec3_dot(min_sr.normal, sample_ray.direction) / pdf);
                vec3_mult(indirect_illum, indirect_illum, fr);                
                vec3_add(radiance, radiance, indirect_illum);
            }
        }        
    }else
    {
        vec3_copy(radiance, sl->bg_color);
    }
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

                if(min_sr.mat->mat_type == PHONG)
                {
                    vec3 reflect_dir;
                    calcReflectRayDir(reflect_dir, min_sr.normal, ray.direction);

                    vec3 new_sample;
                    getNextSample3D(new_sample, min_sr.mat->h_samples);
                    Ray sample_ray;
                    vec3_copy(sample_ray.origin, min_sr.hit_point);
                    getVec3InLocalBasis(sample_ray.direction, new_sample, reflect_dir);
                    if(vec3_dot(sample_ray.direction, min_sr.normal) < 0)
                    {
                        vec3 reflected_sample = {-new_sample[0], -new_sample[1], new_sample[2]};
                        getVec3InLocalBasis(sample_ray.direction, reflected_sample, reflect_dir);
                    }
                    double phong_lobe = pow(vec3_dot(sample_ray.direction, reflect_dir), min_sr.mat->exp);

                    float pdf = phong_lobe * (vec3_dot(min_sr.normal, sample_ray.direction));

                    vec3 reflectance;
                    vec3_scale(reflectance, min_sr.mat->cs, min_sr.mat->ks);                
                    vec3 fr;
                    vec3_scale(fr, reflectance, phong_lobe);

                    vec3 indirect_illum = {0.0f, 0.0f, 0.0f};
                    pathTrace(indirect_illum, depth-1, h_sample, sample_ray, so, sl);
                
                    vec3_scale(indirect_illum, indirect_illum, vec3_dot(min_sr.normal, sample_ray.direction) / pdf);
                    vec3_mult(indirect_illum, indirect_illum, fr);                
                    vec3_add(radiance, radiance, indirect_illum);
                }
            }
        }
    }else
    {
        vec3_copy(radiance, sl->bg_color);
    }
}
