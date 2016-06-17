#pragma once

#include "sampling.h"
#include "util/vec.h"
#include "util/ray.h"
#include "util/shaderec.h"
#include "util/constants.h"
#include "sceneobj.h"
#include "shading.h"
#include "intersect.h"

Samples3D other_samples = Samples3D_default;

void initOtherSampler(const int num_samples, const int num_sets)
{
    Samples2D unit_square_samples = Samples2D_default;
    Samples2D disk_samples = Samples2D_default;
    genMultijitteredSamples(&unit_square_samples, num_samples, num_sets);
    mapSamplesToDisk(&disk_samples, &unit_square_samples);
    mapSamplesToHemisphere(&other_samples, &disk_samples, 1);
    freeSamples2D(&unit_square_samples);
    freeSamples2D(&disk_samples);    
}

void traceRay(vec3 radiance, vec3 sample, const Ray ray, const SceneObjects* so, const SceneLights* sl)
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

void whittedTraceRay(vec3 radiance, int depth, const vec3 h_sample,
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
            // Add ambient component to radiance
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
            // TODO: experiment with branch length

            if(depth > 0 && min_sr.mat->mat_type == PHONG)
            {
                vec3 reflectance;
                vec3_scale(reflectance, min_sr.mat->cs, min_sr.mat->ks);
                vec3 reflect_dir;
                calcReflectRayDir(reflect_dir, min_sr.normal, ray.direction);

                vec3 new_sample;
                getNextSample3D(new_sample, &other_samples);
                Ray sample_ray;
                vec3_copy(sample_ray.origin, min_sr.hit_point);
                getVec3InLocalBasis(sample_ray.direction, new_sample, reflect_dir);
                if(vec3_dot(sample_ray.direction, min_sr.normal) < 0)
                {
                    vec3 reflected_sample = {-new_sample[0], -new_sample[1], new_sample[2]};
                    getVec3InLocalBasis(sample_ray.direction, reflected_sample, reflect_dir);
                }
                // Probblem: division by zero
                double phong_lobe = pow(vec3_dot(sample_ray.direction, reflect_dir), min_sr.mat->exp);
                //printf("phong lobe %f\n", phong_lobe);

                float pdf = phong_lobe * (vec3_dot(min_sr.normal, sample_ray.direction));

                vec3 fr;
                vec3_scale(fr, reflectance, phong_lobe);

                vec3 indirect_illum = {0.0f, 0.0f, 0.0f};
                whittedTraceRay(indirect_illum, depth-1, h_sample, sample_ray, so, sl);
                
                vec3 tmp;
                //vec3_mult(tmp, fr, indirect_illum);
                if(vec3_dot(min_sr.normal, sample_ray.direction) < 0.0f)
                {
                    printf("dot negative\n");
                }
                vec3_scale(tmp, indirect_illum, vec3_dot(min_sr.normal, sample_ray.direction) / pdf);
                vec3_mult(tmp, tmp, fr);
                /*
                if(tmp[0] <= 0.0f || tmp[1] <= 0.0f || tmp[2] <= 0.0f)
                {
                    printf("tmp negative\n");
                }
                */
                if(phong_lobe == 0.0f)
                {
                    printf("phong lobe is zero\n");
                }
                //vec3_add(radiance, radiance, indirect_illum);
                vec3_add(radiance, radiance, tmp);

            }
        }        
    }else
    {
        vec3_copy(radiance, sl->bg_color);
    }
}
