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

float raycast(vec3, int, const vec3, const Ray, const SceneObjects*, const SceneLights*, const int);
float whittedTrace(vec3, int, const vec3, const Ray, const SceneObjects*, const SceneLights*, const int);
float pathTrace(vec3, int, const vec3, const Ray, const SceneObjects*, const SceneLights*, const int);

typedef float (*traceFunc)(vec3, int, const vec3, const Ray, const SceneObjects*, const SceneLights*, const int);

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

float raycast(vec3 radiance, int depth, const vec3 sample, const Ray ray,
              const SceneObjects* so, const SceneLights* sl, const int sample_index)
{
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
        }else
        {
            // Add ambient component to radiance
            ambientShading(radiance, sl->amb_light, sample, so, &min_sr);

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
                          const int depth,  const vec3 h_sample, const Ray ray, const ShadeRec* sr,
                          const SceneObjects* so, const SceneLights*  sl, const int sample_index)
{
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
    float t = whittedTrace(spec_ref_radiance, depth-1, h_sample, sample_ray, so, sl, sample_index);
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

float whittedTrace(vec3 radiance, int depth, const vec3 h_sample,
                   const Ray ray, const SceneObjects* so, const SceneLights* sl, const int sample_index)
{
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
                reflect_t = calcSpecRefRadiance(reflected_illum, depth, h_sample, ray, &min_sr, so, sl, sample_index);

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
                    transmit_t = whittedTrace(transmitted_illum, depth-1, h_sample, transmitted_ray, so, sl, sample_index);
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

float calcSpecRadiancePT(vec4 ref_radiance, const Ray ray, const ShadeRec* sr, const vec3 h_sample,
                         const int depth, const SceneObjects* so, const SceneLights* sl, const int sample_index)
{
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
    t = pathTrace(ref_radiance, depth-1, h_sample, sample_ray, so, sl, sample_index);
    vec3_scale(ref_radiance, ref_radiance, vec3_dot(sr->normal, sample_ray.direction) / pdf);
    vec3_mult(ref_radiance, ref_radiance, brdf);
    return t;
}

float pathTrace(vec3 radiance, int depth, const vec3 h_sample, const Ray ray, const SceneObjects* so,
                const SceneLights* sl, const int sample_index)
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
                    pathTrace(inc_radiance, depth-1, h_sample, sample_ray, so, sl, sample_index);
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
                    calcSpecRadiancePT(reflected_illum, ray, &min_sr, h_sample, depth-1, so, sl, sample_index);
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
                        reflect_t = calcSpecRadiancePT(reflected_illum, ray, &min_sr, h_sample, depth, so, sl, sample_index);
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


                        transmit_t = pathTrace(transmitted_illum, depth-1, h_sample, transmitted_ray, so, sl, sample_index);
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

float schlickPhaseFunc(const float cos_theta, const float k)
{
    float numerator = 1.0f - k * k;
    float a = 1.0f + k * cos_theta;
    float denom = 4.0f * (float)PI * (a*a);
    return numerator / denom;
}

void directIllumInScatter(vec3 radiance, const vec3 point, const float extinct_coeff, const float scatter_coeff,
                          const vec3 view_dir,
                          const SceneObjects *so, const SceneLights *sl, const int sample_index)
{
    ShadeRec tmp_sr;
    vec3_copy(tmp_sr.hit_point, point);
    Ray ray;
    vec3_copy(ray.origin, point);
    for(int j = 0; j < sl->num_lights; j++)
    {
        // 4a construct a ray towards a light point
        getLightDir(ray.direction, sl->light_types[j], sl->light_ptrs[j], &tmp_sr, sample_index);
        if(shadowTest(j, sl, so, ray.direction, &tmp_sr)){continue;}
        ShadeRec light_sr;
        // 4b find exiting t value for that ray
        // TODO: handle cases where light is occluded
        float t_light_exit = intersectTest(&light_sr, so, ray);
        // 4c get radiance from light directly
        vec3 light_rad = {0.0f, 0.0f, 0.0f};
        getIncRadiance(light_rad, sl->light_types[j], sl->light_ptrs[j], light_sr.hit_point);
        // 4c adjust light radiance by distance traveled in medium
        float rad_dec = powf((float)K_E, -extinct_coeff * t_light_exit);
        vec3_scale(light_rad, light_rad, rad_dec);
        // 5. Calculate how much direct illum contributes via phase function
        //vec3_scale(light_rad, light_rad, phase_func * scatter_coeff);
        float phase = schlickPhaseFunc(vec3_dot(ray.direction, view_dir), 0.5f);
        vec3_scale(light_rad, light_rad, phase * scatter_coeff);
        vec3_add(radiance, radiance, light_rad);
    }
}

void raymarch(vec3 radiance, const Ray ray, const vec3 h_sample, const SceneObjects *so,
              const SceneLights *sl, const int sample_index)
{
    /*
      TODO:
      light inside medium
      other phase functions
      different scattering and absorption coefficients
      multiple scatterings
     */
    vec3 view_dir;
    vec3_negate(view_dir, ray.direction);
    const float extinct_coeff = 0.01f;
    const float phase_func = 1.0f / (float)(4.0 * PI);
    const float scatter_coeff = extinct_coeff * 0.5f;
    float t_seg = 5.0f;    
    ShadeRec in_sr;
    // 1. FInd exiting t value and location
    float t_exit = intersectTest(&in_sr, so, ray);
    if(t_exit == TMAX){return;}
    if(in_sr.mat->mat_type == PARTICIPATING)
    {
        // 2. Calculate initial radiance
        Ray exit_ray = ray;
        vec3_copy(exit_ray.origin, in_sr.hit_point);
        vec3 init_rad = {0.0f, 0.0f, 0.0f};
        raycast(init_rad, 0, h_sample, exit_ray, so, sl, sample_index);
        // 3. Ray march back to front
        for(float i = t_exit; i > K_EPSILON; i -= t_seg)
        {
            // 4. At each segment, calculate direction illumination
            vec3 point;
            getPointOnRay(point, ray, i);
            vec3 total_light_rad = {0.0f, 0.0f, 0.0f};
            directIllumInScatter(total_light_rad, point, extinct_coeff, scatter_coeff, view_dir,
                                 so, sl, sample_index);
            // Decrease initial radiance
            float rad_dec = powf((float)K_E, -extinct_coeff * t_seg);
            vec3_scale(init_rad, init_rad, rad_dec);
            vec3_add(init_rad, init_rad, total_light_rad);
        }
        vec3_copy(radiance, init_rad);
    }else
    {
        vec3 init_rad = {0.0f, 0.0f, 0.0f};
        if(in_sr.mat->mat_type == EMISSIVE)
        {
            vec3_scale(init_rad, in_sr.mat->ce, in_sr.mat->ke/1.0f);
        }else
        {
            for(int i = 0; i < sl->num_lights; i++)
            {
                vec3 light_dir;
                getLightDir(light_dir, sl->light_types[i], sl->light_ptrs[i], &in_sr, sample_index);
                float ndotwi = vec3_dot(light_dir, in_sr.normal);
                if(ndotwi <= 0){continue;}
                if(shadowTest(i, sl, so, light_dir, &in_sr)){continue;}
                vec3 light_rad;
                getIncRadiance(light_rad, sl->light_types[i], sl->light_ptrs[i], in_sr.hit_point);
                float light_dist = calcLightDistance(sl->light_types[i], sl->light_ptrs[i], in_sr.hit_point);
                Ray light_ray;
                vec3_copy(light_ray.origin, in_sr.hit_point);
                vec3_copy(light_ray.direction, light_dir);
                ShadeRec light_sr;
                float light_t = intersectTest(&light_sr, so, light_ray);
                float dist_in_medium = 0.0f;
                if(light_t > light_dist)
                {
                    // Case where the light inside the medium or the light isn't physical
                    dist_in_medium = light_dist;
                }else if(light_sr.mat->mat_type == PARTICIPATING)
                {
                    dist_in_medium = light_t;
                }
                float rad_dec = powf((float)K_E, -extinct_coeff * dist_in_medium);
                vec3 dir_illum = {0.0f, 0.0f, 0.0f};
                vec3_scale(dir_illum, light_rad, rad_dec);
                vec3_add(init_rad, init_rad, dir_illum);
            }
        }
        // TODO: sometimes t_exit == T_MAX
        for(float i = t_exit; i > K_EPSILON; i -= t_seg)
        {
            // 4. At each segment, calculate direction illumination
            vec3 point;
            getPointOnRay(point, ray, i);
            ShadeRec tmp_sr;
            vec3_copy(tmp_sr.hit_point, point);
            vec3 total_light_rad = {0.0f, 0.0f, 0.0f};
            directIllumInScatter(total_light_rad, point, extinct_coeff, scatter_coeff, view_dir,
                                 so, sl, sample_index);
            // Decrease initial radiance
            float rad_dec = powf((float)K_E, -extinct_coeff * t_seg);
            vec3_scale(init_rad, init_rad, rad_dec);
            vec3_add(init_rad, init_rad, total_light_rad);
        }
        vec3_copy(radiance, init_rad);
    }
}

/*
  Input: starting point, length, ray
 */
void mediumMarch(vec3 out_rad, const vec3 entry_rad, const Ray ray, const float t_seg,
                 const float medium_dist, const float extinct_coeff,
                 const float scatter_coeff, const SceneObjects *so, const SceneLights *sl,
                 const int sample_index)
{
    vec3 init_rad;
    vec3_copy(init_rad, entry_rad);
    vec3 view_dir;
    vec3_negate(view_dir, ray.direction);
    for(float i = medium_dist; i > K_EPSILON; i -= t_seg)
    {
        // 4. At each segment, calculate direction illumination
        vec3 point;
        getPointOnRay(point, ray, i);
        vec3 total_light_rad = {0.0f, 0.0f, 0.0f};
        directIllumInScatter(total_light_rad, point, extinct_coeff, scatter_coeff, view_dir,
                             so, sl, sample_index);
        // Decrease initial radiance
        float rad_dec = powf((float)K_E, -extinct_coeff * t_seg);
        vec3_scale(init_rad, init_rad, rad_dec);
        vec3_scale(total_light_rad, total_light_rad, (rad_dec * t_seg) + (1.0 - rad_dec)*0.5f*t_seg);
        vec3_add(init_rad, init_rad, total_light_rad);
    }
    vec3_copy(out_rad, init_rad);
}

void fogmarch(vec3 radiance, const Ray ray, const vec3 h_sample, const SceneObjects *so,
              const SceneLights *sl, const int sample_index)
{
    vec3_copy(radiance, BLACK);
    vec3 view_dir;
    vec3_negate(view_dir, ray.direction);
    const float extinct_coeff = 0.002f;
    const float scatter_coeff = extinct_coeff * 0.6f;
    float t_seg = 5.0f;
    ShadeRec min_sr;
    float min_t = intersectTest(&min_sr, so, ray);

    if(min_t != TMAX)
    {
        vec3 init_rad = {0.0f, 0.0f, 0.0f};
        if(min_sr.mat->mat_type == EMISSIVE)
        {
            vec3_scale(init_rad, min_sr.mat->ce, min_sr.mat->ke/1.0f);
        }else
        {
            // Calc surface direct illum
            for(int i = 0; i < sl->num_lights; i++)
            {
                vec3 light_dir;
                getLightDir(light_dir, sl->light_types[i], sl->light_ptrs[i], &min_sr, sample_index);
                float ndotwi = vec3_dot(light_dir, min_sr.normal);
                if(ndotwi <= 0){continue;}
                if(shadowTest(i, sl, so, light_dir, &min_sr)){continue;}
                vec3 light_rad;
                getIncRadiance(light_rad, sl->light_types[i], sl->light_ptrs[i], min_sr.hit_point);
                float light_dist = calcLightDistance(sl->light_types[i], sl->light_ptrs[i], min_sr.hit_point);
                float rad_dec = powf((float)K_E, -extinct_coeff * light_dist);
                vec3 dir_illum = {0.0f, 0.0f, 0.0f};
                vec3_scale(dir_illum, light_rad, rad_dec);
                vec3_add(init_rad, init_rad, dir_illum);
            }
        }
        // ray march back to camera
        mediumMarch(radiance, init_rad, ray, t_seg, min_t, extinct_coeff, scatter_coeff,
                    so, sl, sample_index);
        /*
        for(float i = min_t; i > K_EPSILON; i -= t_seg)
        {
            // 4. At each segment, calculate direction illumination
            vec3 point;
            getPointOnRay(point, ray, i);
            vec3 total_light_rad = {0.0f, 0.0f, 0.0f};
            directIllumInScatter(total_light_rad, point, extinct_coeff, scatter_coeff, view_dir,
                                 so, sl, sample_index);
            // Decrease initial radiance
            float rad_dec = powf((float)K_E, -extinct_coeff * t_seg);
            vec3_scale(init_rad, init_rad, rad_dec);
            vec3_scale(total_light_rad, total_light_rad, (rad_dec * t_seg) + (1.0 - rad_dec)*0.5f*t_seg);
            vec3_add(init_rad, init_rad, total_light_rad);
        }
        vec3_copy(radiance, init_rad);
        */
    }else
    {
        const float max_fog_dist = 1000.0f;
        vec3 init_rad = {0.0f, 0.0f, 0.0f};
        mediumMarch(radiance, init_rad, ray, t_seg, max_fog_dist, extinct_coeff, scatter_coeff,
                    so, sl, sample_index);
        /*
        for(float i = max_fog_dist; i > K_EPSILON; i -= t_seg)
        {
            // 4. At each segment, calculate direction illumination
            vec3 displacement;
            vec3_scale(displacement, ray.direction, i);
            Ray light_ray;
            vec3_add(light_ray.origin, ray.origin, displacement);

            vec3 total_light_rad = {0.0f, 0.0f, 0.0f};
            directIllumInScatter(total_light_rad, light_ray.origin, extinct_coeff, scatter_coeff, view_dir,
                                 so, sl, sample_index);
            // Decrease initial radiance
            float rad_dec = powf((float)K_E, -extinct_coeff * t_seg);
            vec3_scale(init_rad, init_rad, rad_dec);
            vec3_scale(total_light_rad, total_light_rad, (rad_dec * t_seg) + (1.0 - rad_dec)*0.5f*t_seg);
            vec3_add(init_rad, init_rad, total_light_rad);
        }
        vec3_copy(radiance, init_rad);
        */
    }
}

void mediumTrace(vec3 radiance, const Ray ray, const vec3 h_sample, const SceneObjects *so,
                 const SceneLights *sl, const int sample_index)
{
    vec3 view_dir;
    vec3_negate(view_dir, ray.direction);
    const float extinct_coeff = 0.01f;
    const float scatter_coeff = extinct_coeff * 0.5f;
    float t_seg = 5.0f;
    ShadeRec sr;
    float t = intersectTest(&sr, so, ray);
    if(t == TMAX){return;} // Shouldn't happen if the ray is inside a medium
    if(sr.mat->mat_type == PARTICIPATING) // Ray passes through the medium hitting nothing
    {
        Ray new_ray;
        vec3_copy(new_ray.direction, ray.direction);
        getPointOnRay(new_ray.origin, ray, t);

        vec3 entry_rad = {0.0f, 0.0f, 0.0f};
        //trace(entry_rad, new_ray, h_sample, so, sl, sample_index);
    }else // Ray hits something inside the medium
    {

    }
}
