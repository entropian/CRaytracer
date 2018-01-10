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
#include "accelerator/bvh.h"

extern int MAX_DEPTH;
#define SEPARATE_DIRECT_INDIRECT

enum TraceType
{
    RAYCAST,
    WHITTED,
    PATHTRACE
    //PHOTONMAP
};

typedef struct TraceArgs_s
{
    Sampler* sampler;
    
    // Data that are constant throughout a single sample
    SceneObjects *objects;
    SceneLights *lights;
}TraceArgs;


float raycast(vec3, int, const Ray, TraceArgs *trace_args);
float whittedTrace(vec3, int, const Ray, TraceArgs *trace_args);
float pathTrace(vec3, int, const Ray, TraceArgs *trace_args);
//float whittedPhotonTrace(vec3, int, const Ray, TraceArgs *trace_args);
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
        /*
    case RAYCAST:
        func = &raycast;
        break;
    case WHITTED:
        func = &whittedTrace;
        break;
        */
    case PATHTRACE:
        func = &pathTrace;
        break;
        /*
    case PHOTONMAP:
        func = &whittedPhotonTrace;
        break;
        */
    default:
        func = &pathTrace;
    }
    return func;
}

bool totalInternalReflection(const ShadeRec* sr, const float ior_in, const float ior_out)
{
    float cos_theta_i = vec3_dot(sr->normal, sr->wo);
    //float eta = sr->mat.ior_in / sr->mat.ior_out;
    float eta = ior_in / ior_out;

    if(cos_theta_i < 0.0f)
    {
        eta = 1.0f / eta;
    }
    return (1.0f - (1.0f - cos_theta_i * cos_theta_i) / (eta * eta) < 0.0f);
}

float calcTransmitDir(vec3 transmit_dir, const ShadeRec* sr)
{
    if(sr->mat.mat_type != TRANSPARENT)
    {
        fprintf(stderr, "Not refractive material.\n");
        return 0.0f;
    }
    Transparent* trans = (Transparent*)sr->mat.data;
    //return calcTransmitDir(transmit_dir, sr->normal, sr->wo, sr->mat.ior_in, sr->mat.ior_out);
    return calcTransmitDir(transmit_dir, sr->normal, sr->wo, trans->ior_in, trans->ior_out);
}

// Maybe this could be deleted
float calcFresnelReflectance(const ShadeRec* sr)
{
    if(sr->mat.mat_type != TRANSPARENT)
    {
        fprintf(stderr, "Not refractive material.\n");
        return 0.0f;
    }
    Transparent* trans = (Transparent*)sr->mat.data;    
    vec3 n;
    vec3_copy(n, sr->normal);
    // NOTE: not sure about the next line
    float cos_theta_i = vec3_dot(n, sr->wo);
    float eta = trans->ior_in / trans->ior_out;    
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

void transformToLocalBasis(vec3 r, const vec3 a, const BSDF* bsdf)
{
    orthoNormalTransform(r, bsdf->tangent, bsdf->binormal, bsdf->normal, a);
}

void estimateDirect(vec3 L, const vec2 light_sample, const vec2 scatter_sample,
                    const ShadeRec* sr, const int light_index, const SceneLights* sl, const SceneObjects* so)
{
    vec3_copy(L, BLACK);
    vec3 Li, wi;
    // Sample light source
    float light_pdf, scatter_pdf, t;
    vec3 sample_point, sample_normal;
    switch(sl->light_types[light_index])
    {
    case AREALIGHT:
    {
        AreaLight* area_light = (AreaLight*)(sl->light_ptrs[light_index]);
        if(area_light->obj_type == SPHERE)
        {
            // TODO: use new tangent basis for sampling
            Sphere* sphere = (Sphere*)(area_light->obj_ptr);
            vec3 hit_point_to_center;
            vec3_sub(hit_point_to_center, sr->hit_point, sphere->center);
            vec3 z_axis;
            vec3_normalize(z_axis, hit_point_to_center);
            vec3 h_sample;
            mapSampleToHemisphere(h_sample, light_sample);
            getVec3InLocalBasis(h_sample, h_sample, z_axis); // NOTE: hackish
            vec3_scale(sample_point, h_sample, sphere->radius);
            vec3_add(sample_point, sample_point, sphere->center);
            vec3_copy(sample_normal, h_sample);
            //light_pdf = fabs(vec3_dot(h_sample, z_axis)) * INV_PI;
            light_pdf = 1.0f / (2.0f * PI * sphere->radius*sphere->radius)
                * fabs(vec3_dot(h_sample, z_axis)) * INV_PI;
        }else if(area_light->obj_type == RECTANGLE)
        {
            Rectangle* rect = (Rectangle*)(area_light->obj_ptr);
            vec3 u, v;
            vec3_scale(u, rect->width, light_sample[0]);
            vec3_scale(v, rect->height, light_sample[1]);
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
        light_pdf *= vec3_dot(sample_to_hit_point, sample_to_hit_point)
            / fabs(vec3_dot(sample_normal, neg_wi));
        // Calculate Li
        vec3_scale(Li, area_light->color, area_light->intensity);
    } break;
    case ENVLIGHT:
    {
        // TODO: better sampling method
        EnvLight* env_light = (EnvLight*)(sl->light_ptrs[light_index]);
        vec3 h_sample;
        mapSampleToHemisphere(h_sample, light_sample);
        transformToLocalBasis(h_sample, h_sample, &(sr->bsdf));
        vec3 displacement;
        vec3_scale(displacement, h_sample, env_light->world_radius);
        vec3_add(sample_point, sr->hit_point, displacement);
        vec3_negate(sample_normal, h_sample);
        vec3_copy(wi, h_sample);
        light_pdf = fabs(vec3_dot(h_sample, sr->normal)) * INV_PI;
        vec3_scale(Li, env_light->color, env_light->intensity);
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
        bool in_shadow;
        float distance = vec3_length(sample_to_hit_point);
        Ray shadow_ray;
        resetRay(&shadow_ray, sr->hit_point, wi);
        float t = shadowIntersectTest(so, shadow_ray, distance);
        ShadeRec tmp_sr;
        float min_t = intersectTest(&tmp_sr, so, shadow_ray);
        if(t == TMAX && min_t < TMAX)
        {
            t = shadowIntersectTest(so, shadow_ray, distance);
            min_t = intersectTest(&tmp_sr, so, shadow_ray);
        }
        in_shadow = t < distance - K_EPSILON;
        if(!in_shadow)
        {
            vec3 radiance;
            vec3_mult(radiance, f, Li);
            vec3_scale(radiance, radiance, 1.0f / light_pdf);
            vec3_add(L, L, radiance);
        }
    }
}

void uniformSampleOneLight(vec3 L, const vec2 light_sample, const vec2 scatter_sample,
                           const ShadeRec* sr, const SceneLights* sl, const SceneObjects* so)
{
    vec3_copy(L, BLACK);
    int num_lights = sl->num_lights;
    if(num_lights == 0)
    {
        return;
    }
    int light_index;
    float light_pdf = 1.0f / num_lights;
    float rand_float = (float)rand() / (float)RAND_MAX;
    light_index = (int)min(rand_float * num_lights, num_lights - 1);
    vec3 Ld;
    estimateDirect(Ld, light_sample, scatter_sample, sr, light_index, sl, so);
    vec3_scale(L, Ld, 1.0f / light_pdf);
}

float pathTrace(vec3 radiance, int depth, const Ray primary_ray, TraceArgs *trace_args)
{
    const SceneObjects *so = trace_args->objects;
    const SceneLights *sl = trace_args->lights;
    Sampler* sampler = trace_args->sampler;

    // TODO: calculate world radius for env light
    if(sl->env_light)
    {
        if(so->accel == BVH)
        {
            BVHNode* node = (BVHNode*)(so->accel_ptr);
            vec3 dist;
            vec3_sub(dist, node->aabb.max, node->aabb.min);
            //sl->env_light->world_radius = vec3_length(dist) * 0.5f;
            sl->env_light->world_radius = vec3_length(dist) * 2.0f;
        }else if(so->accel == GRID)
        {
            UniformGrid* grid = (UniformGrid*)(so->accel_ptr);
            vec3 dist;
            vec3_sub(dist, grid->aabb.max, grid->aabb.min);
            //sl->env_light->world_radius = vec3_length(dist) * 0.5f;
            sl->env_light->world_radius = vec3_length(dist) * 2.0f;
        }else if(so->accel == BVH4)
        {
            // TODO
        }
    }
    
    
    vec3 L = {0.0f, 0.0f, 0.0f};
    vec3 beta = {1.0f, 1.0f, 1.0f};
    Ray ray = primary_ray;
    bool specular_bounce = false;
    int bounces;
    for(bounces = 0; ; bounces++)
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
                    Emissive* emissive = (Emissive*)sr.mat.data;
                    vec3 inc_radiance;
                    vec3_scale(inc_radiance, emissive->color, emissive->intensity);
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
        vec2 light_sample, scatter_sample;
        Sampler_getSample(light_sample, sampler);
        Sampler_getSample(scatter_sample, sampler);
        if(!(sr.mat.mat_type == REFLECTIVE || sr.mat.mat_type == TRANSPARENT))
        {
            vec3 contrib;        
            uniformSampleOneLight(contrib, light_sample, scatter_sample,  &sr, sl, so);
            vec3_mult(contrib, contrib, beta);
            vec3_add(L, L, contrib);
        }

        // Sample BSDF to get new path direction
        vec3 wo, wi, f;
        vec3_negate(wo, ray.direction);
        vec2 sample;
        Sampler_getSample(sample, sampler);
        float pdf = BSDF_sample_f(f, wi, wo, sample, &(sr.bsdf));
        if(vec3_equal(f, BLACK) || pdf == 0.0f)
        {
            BSDF_freeBxDFs(&(sr.bsdf));                    
            break;
        }
        vec3_scale(f, f, fabs(vec3_dot(wi, sr.normal)) / pdf);
        vec3_mult(beta, beta, f);
        if(sr.mat.mat_type == REFLECTIVE || sr.mat.mat_type == TRANSPARENT)
            specular_bounce = true;
        else
            specular_bounce = false;
        resetRay(&ray, sr.hit_point, wi);
        
        // Possibly terminate the path with Russian roulette
        /*
        if(bounces > 3)
        {
            float q = max(0.05, 1.0f - beta[1]);
            float rand_float = (float)rand() / (float)RAND_MAX;
            if(rand_float < q)
            {
                BSDF_freeBxDFs(&(sr.bsdf));
                break;
            }
            vec3_scale(beta, beta, 1.0f / (1.0f - q));
            // TODO Different from PBRT source
            //vec3_scale(beta, beta, (1.0f - q));
        }
        */
        BSDF_freeBxDFs(&(sr.bsdf));
    }
    if(bounces > 0)
        vec3_scale(L, L, 1.0f / bounces);
    vec3_copy(radiance, L);
    return 0.0f;
}
