#pragma once

#include "sampling.h"
#include "util/vec.h"
#include "util/ray.h"
#include "util/shaderec.h"
#include "util/constants.h"
#include "scene/scenedata.h"
#include "shading.h"
#include "intersect.h"
#include <assert.h>

extern int MAX_DEPTH;
#define SEPARATE_DIRECT_INDIRECT

enum TraceType
{
    RAYCAST,
    WHITTED,
    PATHTRACE,
    PHOTONMAP
};

// Definition in photonmap.h
struct Photonmap_s;
typedef struct Photonmap_s Photonmap;
struct PhotonQueryVars_s;
typedef struct PhotonQueryVars_s PhotonQueryVars;
void calcPhotonmapComponent(vec3, const vec3, const PhotonQueryVars*,
                                 const Photonmap*, const Photonmap*,
                                 const SceneObjects*, const ShadeRec*, const vec3);

typedef struct TraceArgs_s
{
    Sampler* sampler;
    
    // Data that are constant throughout a single sample
    SceneObjects *objects;
    SceneLights *lights;

    // Data used for photon mapping
    Photonmap *photon_map;
    Photonmap *caustic_map;
    PhotonQueryVars *query_vars;
    vec3 caustic_rad;

    // Currently not really used
    Material* medium_mat;
}TraceArgs;


float raycast(vec3, int, const Ray, TraceArgs *trace_args);
float whittedTrace(vec3, int, const Ray, TraceArgs *trace_args);
float pathTrace(vec3, int, const Ray, TraceArgs *trace_args);
float whittedPhotonTrace(vec3, int, const Ray, TraceArgs *trace_args);
float raycastMedium(vec3, int, const Ray, TraceArgs *trace_args);
float whittedTraceMedium(vec3, int, const Ray, TraceArgs *trace_args);

/*
typedef float (*traceFunc)(vec3, int, const vec3, const Ray,
 const SceneObjects*, const SceneLights*, const int, const unsigned int);
*/
typedef float (*traceFunc)(vec3, int, const Ray, TraceArgs *trace_args);

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
    case PHOTONMAP:
        func = &whittedPhotonTrace;
        break;
    default:
        func = &raycast;
    }
    return func;
}

float raycast(vec3 radiance, int depth, const Ray ray, TraceArgs *trace_args)
{
    const SceneObjects *so = trace_args->objects;
    const SceneLights *sl = trace_args->lights;
    Sampler* sampler = trace_args->sampler;
    vec2 sample;
    Sampler_getSample(sample, sampler);
    vec3 h_sample;
    mapSampleToHemisphere(h_sample, sample);

    vec3_copy(radiance, BLACK);
    ShadeRec min_sr;
    float min_t = intersectTest(&min_sr, so, ray);

    // Shading
    if(min_t < TMAX)
    {
        if(min_sr.mat.tex_flags != NO_TEXTURE)
        {
            updateShadeRecWithTexInfo(&min_sr);
        }
        if(min_sr.mat.mat_type == EMISSIVE)
        {
            vec3_scale(radiance, min_sr.mat.ce, min_sr.mat.ke/1.0f);
            maxToOne(radiance, radiance);
        }else if(min_sr.mat.mat_type == PARTICIPATING)
        {
            /*
            TraceArgs new_trace_args = trace_args;
            new_trace_args.medium_mat = min_sr.mat;
            Ray new_ray;
            getPointOnRay(new_ray.origin, ray, min_t);
            vec3_copy(new_ray.direction, ray.direction);
            raycastMedium(radiance, depth-1, new_ray, new_trace_args);
            */
        }else
        {
            // Add ambient component to radiance
            ambientShading(radiance, sl->amb_light, h_sample, so, &min_sr);

            // Direct illumination
            for(int i = 0; i < sl->num_lights; i++)
            {
                vec3 light_dir;
                vec2 sample;
                Sampler_getSample(sample, sampler);
                getLightDir(light_dir, sl->light_types[i], sl->light_ptrs[i], &min_sr, sample); 
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
                          TraceArgs *trace_args)
{
    Sampler* sampler = trace_args->sampler;
    vec3 reflect_dir, normal;
    vec3_copy(normal, sr->normal);
    calcReflectRayDir(reflect_dir, normal, ray.direction);
    if(vec3_dot(normal, reflect_dir) < 0.0f && sr->mat.mat_type == TRANSPARENT)
    {
        vec3_negate(normal, normal);
    }
    vec3 new_sample;
    Sampler_getCosPowerSample(new_sample, sampler, sr->mat.exp);
    Ray sample_ray;
    vec3_copy(sample_ray.origin, sr->hit_point);
    getVec3InLocalBasis(sample_ray.direction, new_sample, reflect_dir);
    if(vec3_dot(sample_ray.direction, normal) < 0.0f)
    {
        vec3 reflected_sample = {-new_sample[0], -new_sample[1], new_sample[2]};
        getVec3InLocalBasis(sample_ray.direction, reflected_sample, reflect_dir);
    }
    float raydotr = vec3_dot(sample_ray.direction, reflect_dir);
    float phong_lobe = pow(vec3_dot(sample_ray.direction, reflect_dir), sr->mat.exp);
    float pdf = phong_lobe * (vec3_dot(normal, sample_ray.direction));
    vec3 brdf;
    vec3_scale(brdf, sr->mat.cs, phong_lobe);   // Deferring reflectance scaling
    float t = whittedTrace(spec_ref_radiance, depth-1, sample_ray, trace_args);
    vec3_mult(spec_ref_radiance, spec_ref_radiance, brdf);
    vec3_scale(spec_ref_radiance, spec_ref_radiance, vec3_dot(normal, sample_ray.direction) / pdf);
    return t;
}

bool totalInternalReflection(const ShadeRec* sr)
{
    float cos_theta_i = vec3_dot(sr->normal, sr->wo);
    float eta = sr->mat.ior_in / sr->mat.ior_out;

    if(cos_theta_i < 0.0f)
    {
        eta = 1.0f / eta;
    }
    return (1.0f - (1.0f - cos_theta_i * cos_theta_i) / (eta * eta) < 0.0f);
}

float calcTransmitDir(vec3 transmit_dir, const ShadeRec* sr)
{
    /*
    vec3 n;
    vec3_copy(n, sr->normal);
    float cos_theta_i = vec3_dot(n, sr->wo);
    float eta = sr->mat.ior_in / sr->mat.ior_out;
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
    */
    return calcTransmitDir(transmit_dir, sr->normal, sr->wo, sr->mat.ior_in, sr->mat.ior_out);
}

float calcFresnelReflectance(const ShadeRec* sr)
{
    vec3 n;
    vec3_copy(n, sr->normal);
    // NOTE: not sure about the next line
    float cos_theta_i = vec3_dot(n, sr->wo);
    float eta = sr->mat.ior_in / sr->mat.ior_out;    
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

void whittedShade(vec3 radiance, int depth, const Ray ray, TraceArgs *trace_args, const ShadeRec *sr)
{
    const SceneObjects *so = trace_args->objects;
    const SceneLights *sl = trace_args->lights;
    Sampler* sampler = trace_args->sampler;
    if(sr->mat.mat_type == EMISSIVE)
    {
        vec3_scale(radiance, sr->mat.ce, sr->mat.ke/1.0f);
        if(depth == MAX_DEPTH)
        {
            maxToOne(radiance, radiance);
        }
    }else
    {
        vec3 h_sample;
        Sampler_getHemisphereSample(h_sample, sampler);
        ambientShading(radiance, sl->amb_light, h_sample, so, sr);
        // Direct illumination
        // TODO: refactor
        if(sr->mat.mat_type == MATTE || sr->mat.mat_type == PHONG)
        {
            for(int i = 0; i < sl->num_lights; i++)
            {
                vec3 light_dir;
                vec2 sample;
                Sampler_getSample(sample, sampler);
                getLightDir(light_dir, sl->light_types[i], sl->light_ptrs[i], sr, sample);
                float ndotwi = vec3_dot(light_dir, sr->normal);
                if(ndotwi > 0.0f)
                {
                    bool in_shadow = shadowTest(i, sl, so, light_dir, sr);
                    if(!in_shadow)
                    {
                        directIllumShading(radiance, ndotwi, light_dir, sl->light_ptrs[i],
                                           sl->light_types[i], sr);
                    }
                }
            }
        }

        // Indirect illumination
        if(depth > 0 && (sr->mat.mat_type == REFLECTIVE || sr->mat.mat_type == TRANSPARENT))
        {
            float reflect_t;
            vec3 reflected_illum = {0.0f, 0.0f, 0.0f};
            reflect_t = calcSpecRefRadiance(reflected_illum, depth, ray, sr, trace_args);

            if(sr->mat.mat_type == REFLECTIVE)
            {
                vec3_scale(reflected_illum, reflected_illum, sr->mat.ks);
                vec3_add(radiance, radiance, reflected_illum);
                // Problem?
                return;
            }
            float ndotwo = vec3_dot(sr->normal, sr->wo);
            if(!totalInternalReflection(sr))
            {
                // Transmitted radiance
                vec3 transmit_dir = {0, 0, 0};
                float eta = calcTransmitDir(transmit_dir, sr);
                float ndotwt = fabs(vec3_dot(sr->normal, transmit_dir));
                float kr = calcFresnelReflectance(sr);
                float kt = 1.0f - kr;

                vec3 btdf;
                vec3_scale(btdf, WHITE, kt / (eta*eta) / ndotwt);
                Ray transmitted_ray;
                vec3_copy(transmitted_ray.origin, sr->hit_point);
                vec3_copy(transmitted_ray.direction, transmit_dir);

                float transmit_t;
                vec3 transmitted_illum = {0.0f, 0.0f, 0.0f};
                transmit_t = whittedTrace(transmitted_illum, depth-1, transmitted_ray, trace_args);

                vec3_scale(transmitted_illum, transmitted_illum, ndotwt);
                vec3_mult(transmitted_illum, transmitted_illum, btdf);
                // Scaling reflection since there's no total internal reflection
                vec3_scale(reflected_illum, reflected_illum, kr);

                vec3 color_filter_ref, color_filter_trans;
                if(ndotwo > 0.0f)
                {
                    vec3_pow(color_filter_ref, sr->mat.cf_out, reflect_t);
                    vec3_pow(color_filter_trans, sr->mat.cf_in, transmit_t);
                }else
                {
                    vec3_pow(color_filter_ref, sr->mat.cf_in, reflect_t);
                    vec3_pow(color_filter_trans, sr->mat.cf_out, transmit_t);
                }
                vec3_mult(reflected_illum, reflected_illum, color_filter_ref);
                vec3_mult(transmitted_illum, transmitted_illum, color_filter_trans);
                vec3_add(radiance, radiance, transmitted_illum);
            }else
            {
                vec3 color_filter;
                if(ndotwo > 0.0f)
                {
                    vec3_pow(color_filter, sr->mat.cf_out, reflect_t);
                }else
                {
                    vec3_pow(color_filter, sr->mat.cf_in, reflect_t);
                }
                vec3_mult(reflected_illum, reflected_illum, color_filter);
            }
            vec3_add(radiance, radiance, reflected_illum);
        }
    }
}

float whittedTrace(vec3 radiance, int depth, const Ray ray, TraceArgs *trace_args)
{
    const SceneObjects *so = trace_args->objects;
    const SceneLights *sl = trace_args->lights;
    vec3 h_sample;

    vec3_copy(radiance, ORIGIN);
    float min_t = TMAX;
    ShadeRec min_sr;
    min_t = intersectTest(&min_sr, so, ray);
    // Shading
    if(min_t < TMAX)
    {
        if(min_sr.mat.tex_flags != NO_TEXTURE)
        {
            updateShadeRecWithTexInfo(&min_sr);
        }
        whittedShade(radiance, depth, ray, trace_args, &min_sr);
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

float whittedPhotonTrace(vec3 radiance, int depth, const Ray ray, TraceArgs *trace_args)
{
    const SceneObjects *so = trace_args->objects;
    const SceneLights *sl = trace_args->lights;
    Sampler* sampler = trace_args->sampler;
    vec3_copy(radiance, ORIGIN);
    float min_t = TMAX;
    ShadeRec min_sr;
    min_t = intersectTest(&min_sr, so, ray);
    // Shading
    if(min_t < TMAX)
    {
        if(min_sr.mat.tex_flags != NO_TEXTURE)
        {
            updateShadeRecWithTexInfo(&min_sr);
        }
        whittedShade(radiance, depth, ray, trace_args, &min_sr);
        vec3 pm_comp;
        vec3 h_sample;
        Sampler_getHemisphereSample(h_sample, sampler);
        calcPhotonmapComponent(pm_comp, h_sample, trace_args->query_vars, trace_args->photon_map,
                                    trace_args->caustic_map, so, &min_sr, trace_args->caustic_rad);
        vec3_add(radiance, radiance, pm_comp);
    }else
    {
        vec3_copy(radiance, sl->bg_color);
    }
    return min_t;
}

void computeLocalBasis(ShadeRec* sr)
{
    vec3_copy(sr->bsdf.normal, sr->normal);
    vec3_copy(sr->bsdf.tangent, sr->dpdu);
    float dot_product = vec3_dot(sr->bsdf.normal, sr->bsdf.tangent);
    if(dot_product != 0.0f)
    {
        vec3 adjustment;
        vec3_scale(adjustment, sr->bsdf.normal, dot_product);
        vec3_sub(sr->bsdf.tangent, sr->bsdf.tangent, adjustment);
        vec3_normalize(sr->bsdf.tangent, sr->bsdf.tangent);
    }
    vec3_cross(sr->bsdf.binormal, sr->bsdf.normal, sr->bsdf.tangent);
    vec3_normalize(sr->bsdf.binormal, sr->bsdf.binormal);
}

bool isOrthoNormal(const vec3 u, const vec3 v, const vec3 w)
{
    if(vec3_length(u) - 1.0f > K_EPSILON) return false;
    if(vec3_length(v) - 1.0f > K_EPSILON) return false;
    if(vec3_length(w) - 1.0f > K_EPSILON) return false;
    if(vec3_dot(u, v) > K_EPSILON) return false;
    if(vec3_dot(u, w) > K_EPSILON) return false;
    if(vec3_dot(v, w) > K_EPSILON) return false;
    return true;
}

float pathTrace(vec3 radiance, int depth, const Ray ray, TraceArgs *trace_args)
{
    const SceneObjects *so = trace_args->objects;
    const SceneLights *sl = trace_args->lights;
    Sampler* sampler = trace_args->sampler;

    vec3_copy(radiance, ORIGIN);
    float min_t = TMAX;
    ShadeRec min_sr;
    min_t = intersectTest(&min_sr, so, ray);
    // Shading
    if(min_t < TMAX)
    {
        computeLocalBasis(&(min_sr));
        computeScatteringFunc(&(min_sr.bsdf), min_sr.uv, &(min_sr.mat));
        if(min_sr.mat.tex_flags != NO_TEXTURE)
        {
            // TODO move texture fetch into computeScatteringFunc
            updateShadeRecWithTexInfo(&min_sr);
        }
        if(min_sr.mat.mat_type == EMISSIVE)
        {
            //TODO this
#ifdef SEPARATE_DIRECT_INDIRECT
            if(depth == MAX_DEPTH - 1)
            {
                vec3_copy(radiance, BLACK);
            }else
            {
                float ndotwi = clamp(-vec3_dot(min_sr.normal, ray.direction), 0.0f, 1.0f);
                vec3_scale(radiance, min_sr.mat.ce, min_sr.mat.ke * ndotwi);
            }
#else
            vec3_scale(radiance, min_sr.mat.ce, min_sr.mat.ke);
#endif
            //maxToOne(radiance, radiance);
        }else
        {
            if(depth > 0)
            {
#ifdef SEPARATE_DIRECT_INDIRECT
                // Direct illumination
                if(min_sr.mat.mat_type & (MATTE | PHONG) && depth == MAX_DEPTH)
                {
                    for(int i = 0; i < sl->num_lights; i++)
                    {
                        if(sl->light_types[i] == AREALIGHT)
                        {
                            vec3 light_dir;
                            vec2 sample;
                            Sampler_getSample(sample, sampler);
                            getLightDir(light_dir, sl->light_types[i], sl->light_ptrs[i], &min_sr, sample);
                            float ndotwi = vec3_dot(light_dir, min_sr.normal);
                            if(ndotwi > 0)
                            {
                                bool in_shadow = shadowTest(i, sl, so, light_dir, &min_sr);
                                if(!in_shadow)
                                {
                                    directIllumShadingNew(radiance, ndotwi, light_dir, sl->light_ptrs[i],
                                                       sl->light_types[i], &min_sr);
                                }
                            }
                        }else if(sl->light_types[i] == MESHLIGHT)
                        {
                            MeshLight* mesh_light_ptr = (MeshLight*)(sl->light_ptrs[i]);                            
                            // 1. Generate sample
                            vec3 sample, sample_normal;
                            vec3 hit_to_sample;
                            vec3 sample_to_hit;
                            MeshLight_genSample(sample, sample_normal, mesh_light_ptr);
                            // 2. check if hit surface faces sample point. if not go back to 1
                            vec3_sub(hit_to_sample, min_sr.hit_point, sample);
                            // 3. check if sample surface faces hit point. if not go back to 1
                            vec3_negate(sample_to_hit, hit_to_sample);
                            // 4. check if sample point is visible from hit point. if not, terminate
                            if(vec3_dot(sample_normal, hit_to_sample) > 0.0f &&
                                   vec3_dot(min_sr.normal, sample_to_hit) > 0.0f)
                            {
                                Ray ray;
                                vec3_copy(ray.origin, min_sr.hit_point);
                                vec3_normalize(ray.direction, sample_to_hit);
                                float light_dist = vec3_length(sample_to_hit);
                                float t = shadowIntersectTest(so, ray, light_dist);                            
                                // 5. shade
                                if(min_sr.mat.mat_type == MATTE)
                                {
                                    if(t + K_EPSILON >= light_dist)
                                    {
                                        vec3 f;
                                        BSDF_f(f, ray.direction, min_sr.wo, &(min_sr.bsdf));
                                        // Incident radiance
                                        vec3 inc_radiance = {0.0f, 0.0f, 0.0f};
                                        vec3_scale(inc_radiance, mesh_light_ptr->color, mesh_light_ptr->intensity);
                                        vec3 wi, neg_wi;
                                        vec3_normalize(wi, sample_to_hit);
                                        vec3_normalize(neg_wi, hit_to_sample);
                                        float geo_term = vec3_dot(sample_normal, neg_wi) * vec3_dot(min_sr.normal, wi) /
                                            vec3_dot(sample_to_hit, sample_to_hit);

                                        // f * L * G / pdf
                                        vec3 tmp;
                                        vec3_mult(tmp, f, inc_radiance);
                                        float pdf = 1.0f / mesh_light_ptr->surface_area;
                                        vec3_scale(tmp, tmp, geo_term * 1.0f / pdf);
                                        if(tmp[0] < 0.0f || tmp[1] < 0.0f || tmp[2] < 0.0f)
                                        {
                                            printf("negative\n");
                                        }
                                        vec3_add(radiance, radiance, tmp);
                                    }
                                }
                            }
                        }
                    }
                }
#endif
                //if(min_sr.mat.mat_type == REFLECTIVE)
                if(min_sr.mat.mat_type == REFLECTIVE || min_sr.mat.mat_type == TRANSPARENT) 
                {
                    vec3 f;
                    vec3 wi;
                    vec2 sample;
                    Sampler_getSample(sample, sampler);
                    float pdf = BSDF_sample_f(f, wi, min_sr.wo, sample, &(min_sr.bsdf));
                    if(pdf > 0.0f)
                    {
                        Ray sample_ray;
                        vec3_copy(sample_ray.direction, wi);
                        vec3_copy(sample_ray.origin, min_sr.hit_point);
                        
                        vec3 inc_radiance;
                        pathTrace(inc_radiance, depth-1, sample_ray, trace_args);
                        //float ndotwi = vec3_dot(min_sr.normal, wi);

                        vec3 tmp;
                        vec3_mult(tmp, inc_radiance, f);
                        //vec3_scale(tmp, tmp, ndotwi / pdf);
                        vec3_add(radiance, radiance, tmp);
                    }
                }

                if(min_sr.mat.mat_type == MATTE)
                {
                    // sampleF
                    vec3 f;
                    vec3 wi;
                    vec2 sample;
                    Sampler_getSample(sample, sampler);
                    float pdf = BSDF_sample_f(f, wi, min_sr.wo, sample, &(min_sr.bsdf));
                    BSDF* bsdf = &(min_sr.bsdf);
                    if(!isOrthoNormal(bsdf->tangent, bsdf->binormal, bsdf->normal))
                    {
                        printf("not orthonormal\n");
                    }
                    if(pdf > 0.0f)
                    {

                        Ray sample_ray;
                        vec3_copy(sample_ray.direction, wi);
                        vec3_copy(sample_ray.origin, min_sr.hit_point);
                        
                        vec3 inc_radiance;
                        pathTrace(inc_radiance, depth-1, sample_ray, trace_args);
                        float ndotwi = vec3_dot(min_sr.normal, wi);

                        vec3 tmp;
                        vec3_mult(tmp, inc_radiance, f);
                        vec3_scale(tmp, tmp, ndotwi / pdf);
                        vec3_add(radiance, radiance, tmp);
                    }
                }
                /*
                if(min_sr.mat.mat_type == TRANSPARENT)
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
                        // TODO: since I'm using Russian roulette, it probably doesn't make sense
                        // to multiply by kt?
                        //vec3_scale(btdf, WHITE, kt / (eta*eta) / ndotwt);
                        vec3_scale(btdf, WHITE, 1 / (eta*eta) / ndotwt);
                        Ray transmitted_ray;
                        vec3_copy(transmitted_ray.origin, min_sr.hit_point);
                        vec3_copy(transmitted_ray.direction, transmit_dir);

                        transmit_t = pathTrace(transmitted_illum, depth-1, transmitted_ray, trace_args);
                        assert(transmitted_illum[0] >= 0.0f && transmitted_illum[1] >= 0.0f &&
                               transmitted_illum[2] >= 0.0f);
                        vec3_scale(transmitted_illum, transmitted_illum, ndotwt);
                        vec3_mult(transmitted_illum, transmitted_illum, btdf);
                    }

                    vec3 color_filter_ref, color_filter_trans;
                    if(ndotwo > 0.0f)
                    {
                        vec3_pow(color_filter_ref, min_sr.mat.cf_out, reflect_t);
                        vec3_pow(color_filter_trans, min_sr.mat.cf_in, transmit_t);
                    }else
                    {
                        vec3_pow(color_filter_ref, min_sr.mat.cf_in, reflect_t);
                        vec3_pow(color_filter_trans, min_sr.mat.cf_out, transmit_t);
                    }
                    vec3_mult(reflected_illum, reflected_illum, color_filter_ref);
                    vec3_mult(transmitted_illum, transmitted_illum, color_filter_trans);
                    vec3_add(radiance, radiance, transmitted_illum);
                    vec3_add(radiance, radiance, reflected_illum);
                }
                */
            }else
            {
                //vec3_copy(radiance, sl->bg_color);
            }
        }
        // Release BSDF memory
        BSDF_freeBxDFs(&(min_sr.bsdf));
    }else
    {
        //vec3_copy(radiance, sl->bg_color);
        getEnvLightIncRadiance(radiance, ray.direction, sl->env_light);
    }
    return min_t;
}

bool isBlack(const vec3 c)
{
    if(c[0] > 0.0f) return false;
    if(c[1] > 0.0f) return false;
    if(c[2] > 0.0f) return false;
}

inline float powerHeuristic(int nf, float fPdf, int ng, float gPdf) {
    float f = nf * fPdf, g = ng * gPdf;
    return (f * f) / (f * f + g * g);
}

void estimateDirect(vec3 L, Sampler* sampler,
                    const ShadeRec* sr, const int light_index, const SceneLights* sl, const SceneObjects* so)
{
    vec3_copy(L, BLACK);
    vec3 Li, wi;
    // Sample light source
    // Get point on light source
    float light_pdf, scatter_pdf, t;
    vec2 sample;
    Sampler_getSample(sample, sampler);
    vec3 sample_point, sample_normal;
    switch(sl->light_types[light_index])
    {
    case AREALIGHT:
    {
        AreaLight* area_light = (AreaLight*)(sl->light_ptrs[light_index]);
        if(area_light->obj_type == SPHERE)
        {
            Sphere* sphere = (Sphere*)(area_light->obj_ptr);
            vec3 hit_point_to_center;
            vec3_sub(hit_point_to_center, sr->hit_point, sphere->center);
            vec3 z_axis;
            vec3_normalize(z_axis, hit_point_to_center);
            vec3 h_sample;
            mapSampleToHemisphere(h_sample, sample);
            getVec3InLocalBasis(h_sample, h_sample, z_axis); // NOTE: hackish
            vec3_scale(sample_point, h_sample, sphere->radius);
            vec3_add(sample_point, sample_point, sphere->center);
            vec3_copy(sample_normal, h_sample);
            //light_pdf = 1.0 / (calcSphereArea(sphere) * 0.5f);
            light_pdf = fabs(vec3_dot(h_sample, z_axis)) * INV_PI;
        }else if(area_light->obj_type == RECTANGLE)
        {
            Rectangle* rect = (Rectangle*)(area_light->obj_ptr);
            vec3 u, v;
            vec3_scale(u, rect->width, sample[0]);
            vec3_scale(v, rect->height, sample[1]);
            vec3_add(sample_point, rect->point, u);
            vec3_add(sample_point, sample_point, v);
            vec3_copy(sample_normal, rect->normal);
            light_pdf = 1.0f / (vec3_length(rect->width) * vec3_length(rect->height));
        }
        // Calculate wi
        vec3 sample_to_hit_point;
        vec3_sub(sample_to_hit_point, sample_point, sr->hit_point);
        vec3_normalize(wi, sample_to_hit_point);
        // Calculate PDF
        vec3 neg_wi;
        vec3_negate(neg_wi, wi);
        /*
        light_pdf *= 1.0f / vec3_dot(sample_to_hit_point, sample_to_hit_point)
            * fabs(vec3_dot(sample_normal, neg_wi));
            */
        light_pdf *= vec3_dot(sample_to_hit_point, sample_to_hit_point)
            / fabs(vec3_dot(sample_normal, neg_wi));
        // Calculate Li
        vec3_scale(Li, area_light->color, area_light->intensity);
    } break;
    case ENVLIGHT:
    {
        // TODO
    } break;
    default:        
        printf("Invalid light type.\n");
        return;
    }

    vec3 sample_to_hit_point;
    vec3_sub(sample_to_hit_point, sample_point, sr->hit_point);
    if(vec3_dot(sample_to_hit_point, sample_normal) > 0.0f ||
       vec3_dot(sample_to_hit_point, sr->normal) < 0.0f)
    {
        return;
    }

    vec3 f;
    // Calcuate BRDF
    if(light_pdf > 0.0f && !vec3_equal(Li, BLACK))
    {
        BSDF_f(f, wi, sr->wo, &(sr->bsdf));
        vec3_scale(f, f, fabs(vec3_dot(sr->normal, wi)));
        scatter_pdf = BSDF_pdf(wi, sr->wo, &(sr->bsdf));
    }else
        return;

    if(!vec3_equal(f, BLACK))
    {
        bool in_shadow = shadowTest(light_index, sl, so, wi, sr);
        if(!in_shadow)
        {
            vec3 radiance;
            vec3_mult(radiance, f, Li);
            vec3_scale(radiance, radiance, 1.0f / light_pdf);
            vec3_add(L, L, radiance);
        }
    }
}

void uniformSampleOneLight(vec3 L, Sampler* sampler,
                           const ShadeRec* sr, const SceneLights* sl, const SceneObjects* so)
{
    vec3_copy(L, BLACK);
    int num_lights = sl->num_lights;
    // TODO: excluding env light for now
    /*
    if(sl->env_light)
    {
        num_lights++;
    }
    */
    if(num_lights == 0)
    {
        return;
    }
    int light_index;
    float light_pdf = 1.0f / num_lights;
    float rand_float = (float)rand() / (float)RAND_MAX;
    light_index = (int)min(rand_float * num_lights, num_lights - 1);
    vec3 Ld;
    estimateDirect(Ld, sampler, sr, light_index, sl, so);
    vec3_scale(L, Ld, 1.0f / light_pdf);
}

void pathTraceNew(vec3 radiance, int depth, const Ray primary_ray, TraceArgs *trace_args)
{
    const SceneObjects *so = trace_args->objects;
    const SceneLights *sl = trace_args->lights;
    Sampler* sampler = trace_args->sampler;
    
    vec3 L = {0.0f, 0.0f, 0.0f};
    vec3 beta = {1.0f, 1.0f, 1.0f};
    Ray ray = primary_ray;
    bool specular_bounce = false;
    for(int bounces = 0; ; bounces++)
    {
        // Intersect ray with scene
        ShadeRec sr;
        float t = intersectTest(&sr, so, ray);
        
        // Possibly add emitted light at intersection
        if(bounces == 0 || specular_bounce)
        {
            if(t < TMAX)
            {
                if(sr.mat.mat_type == EMISSIVE)
                {
                    // NOTE: questionable
                    float ndotwi = clamp(-vec3_dot(sr.normal, ray.direction), 0.0f, 1.0f);
                    vec3 inc_radiance;
                    vec3_scale(inc_radiance, sr.mat.ce, sr.mat.ke * ndotwi);
                    vec3_mult(inc_radiance, inc_radiance, beta);
                    vec3_add(L, L, inc_radiance);
                }
            }else
            {
                vec3 env_inc_radiance;
                if(sl->env_light)
                    getEnvLightIncRadiance(env_inc_radiance, ray.direction, sl->env_light);
                vec3_mult(env_inc_radiance, env_inc_radiance, beta);
                vec3_add(L, L, env_inc_radiance);
            }
        }

        // Terminate path if ray escaped or maxDepth was reached
        // NOTE: Only do blackbody emission for now
        if(t == TMAX || bounces >= depth || sr.mat.mat_type == EMISSIVE) 
            break;
        
        // Compute scattering functions
        computeLocalBasis(&sr);
        computeScatteringFunc(&(sr.bsdf), sr.uv, &(sr.mat));

        // Sample illumination from lights to find path contribution
        vec3 contrib;
        uniformSampleOneLight(contrib, sampler, &sr, sl, so);
        vec3_mult(contrib, contrib, beta);
        vec3_add(L, L, contrib);

        // Sample BSDF to get new path direction
        vec3 wo, wi, f;
        vec3_negate(wo, ray.direction);
        vec2 sample;
        Sampler_getSample(sample, sampler);
        float pdf = BSDF_sample_f(f, wi, wo, sample, &(sr.bsdf));

        if(vec3_equal(f, BLACK))
        {
            BSDF_freeBxDFs(&(sr.bsdf));                    
            break;
        }
        vec3_scale(f, f, fabs(vec3_dot(wi, sr.normal)) / pdf);
        vec3_mult(beta, beta, f);
        // TODO
        //specular_bounce = isSpecular(&(sr.mat));
        resetRay(&ray, sr.hit_point, wi);
        //vec3_copy(ray.origin, sr.hit_point);
        //vec3_copy(ray.direction, wi);
        
        // Possibly terminate the path with Russian roulette
        if(bounces > 3)
        {
            float q = max(0.05, 1.0f - beta[1]);
            float rand_float = (float)rand() / (float)RAND_MAX;
            if(rand_float < q)
            {
                BSDF_freeBxDFs(&(sr.bsdf));
                break;
            }
            vec3_scale(beta, beta, 1.0f - q);
        }
        BSDF_freeBxDFs(&(sr.bsdf));                
    }
    vec3_copy(radiance, L);
}
