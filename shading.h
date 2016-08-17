#pragma once

#include "util/constants.h"
#include "util/vec.h"
#include "util/util.h"
#include "materials.h"
#include "shapes/shapes.h"
#include "lights.h"
#include "scenedata.h"
#include "sampling.h"
#include "intersect.h"

void diffuseBRDF(vec3 f, const ShadeRec* sr)
{
    vec3 reflectance;
    vec3_scale(reflectance, sr->mat->cd, sr->mat->kd);
    vec3_scale(f, reflectance, 1.0f/(float)PI);
}

void diffuseShading(vec3 radiance, const vec3 inc_radiance_cos,
                    const ShadeRec* sr)
{
    vec3 f;
    diffuseBRDF(f, sr);
    vec3 diffuse_comp;
    vec3_mult(diffuse_comp, inc_radiance_cos, f);
    vec3_add(radiance, radiance, diffuse_comp);
}

void specularShading(vec3 radiance, const vec3 light_dir,
                     const vec3 inc_radiance_cos, const ShadeRec* sr)
{
    // ks*cs * cos(theta)^exp * inc_radiance_cos
    // theta = angle between light and reflected view vector
    // reflect vector: 2*dot(n, v)*n - v
    vec3 tmp, wo_neg, reflect_dir;
    vec3_negate(wo_neg, sr->wo);
    float ndotwo = vec3_dot(sr->normal, sr->wo);
    vec3_scale(tmp, sr->normal, ndotwo * 2.0f);
    vec3_add(reflect_dir, wo_neg, tmp);    
    float rdotl = clamp(vec3_dot(reflect_dir, light_dir), 0.0f, 1.0f);
    vec3_scale(tmp, sr->mat->cs, sr->mat->ks * (float)pow(rdotl, sr->mat->exp));
    vec3_mult(tmp, inc_radiance_cos, tmp);
    vec3_add(radiance, radiance, tmp);
}

float AOTest(const vec3 h_sample, const SceneObjects *so, const ShadeRec* sr)
{
    //vec3 h_sample;
    //getNextSample3D(h_sample, samples);
    Ray shadow_ray;
    getVec3InLocalBasis(shadow_ray.direction, h_sample, sr->normal);
    vec3_copy(shadow_ray.origin, sr->hit_point);
    return shadowIntersectTest(so, shadow_ray);
}

void ambientShading(vec3 radiance, const AmbientLight* amb_light, const vec3 sample,
                    const SceneObjects* so, const ShadeRec* sr)
{
    if(sr->mat->mat_type != TRANSPARENT)
    {
        // ka*ca * amb_inc_radiance    
        vec3 amb_inc_radiance;
        vec3_scale(amb_inc_radiance, amb_light->color, amb_light->intensity);
        vec3 reflectance;
        vec3_scale(reflectance, sr->mat->ca, sr->mat->ka);
        if(amb_light->amb_occlusion && AOTest(sample, so, sr) < TMAX)
        {
            // Ambient Occlusion                        
            vec3 min_amb;
            vec3_scale(min_amb, amb_light->color, 0.01f);
            vec3_copy(radiance, min_amb);
        }else
        {
            vec3_mult(radiance, amb_inc_radiance, reflectance);
        }
    }
}

void areaLightShading(vec3 radiance, const float ndotwi, const vec3 light_dir, const AreaLight* area_light_ptr, const ShadeRec* sr)
{
    // diffuse brdf
    vec3 f;
    diffuseBRDF(f, sr);    

    // Incident radiance
    vec3 inc_radiance = {0.0f, 0.0f, 0.0f}, neg_wi, displacement, light_normal;
    vec3_sub(displacement, sr->hit_point, area_light_ptr->sample_point);
    vec3_negate(neg_wi, light_dir);
    getAreaLightNormal(light_normal, area_light_ptr, sr->hit_point);
    if(vec3_dot(neg_wi, light_normal) > 0.0f)
    {
        getIncRadiance(inc_radiance, AREALIGHT, area_light_ptr);
    }else
    {
        return;
    }

    // Geometry term                                    
    float geo_term = vec3_dot(light_normal, neg_wi) * ndotwi /
        vec3_dot(displacement, displacement);

    // f * L * G / PDF
    vec3 tmp;
    vec3_mult(tmp, f, inc_radiance);
    vec3_scale(tmp, tmp, geo_term * 1.0f/area_light_ptr->pdf);
    vec3_add(radiance, radiance, tmp);

    // Specular component
    MatType mat_type = sr->mat->mat_type;
    if(mat_type == PHONG || mat_type == REFLECTIVE || mat_type == TRANSPARENT)
    {
        vec3 reflect_dir, wo_neg;
        vec3_negate(wo_neg, sr->wo);
        float ndotwo = vec3_dot(sr->normal, sr->wo);
        vec3_scale(tmp, sr->normal, 2.0f * ndotwo);
        vec3_add(reflect_dir, wo_neg, tmp);
        float rdotl = vec3_dot(reflect_dir, light_dir);

        vec3_scale(f, sr->mat->cs, sr->mat->ks * (float)pow(rdotl, sr->mat->exp));
        vec3_mult(tmp, f, inc_radiance);
        vec3_scale(tmp, tmp, geo_term * 1.0f/area_light_ptr->pdf);
        vec3_add(radiance, radiance, tmp);
    }
}

void envLightShading(vec3 radiance, const float ndotwi, const vec3 light_dir, const EnvLight* env_light_ptr, const ShadeRec* sr)
{
    vec3 f, inc_radiance_cos, tmp;
    diffuseBRDF(f, sr);
    getIncRadiance(tmp, ENVLIGHT, env_light_ptr);
    vec3_scale(inc_radiance_cos, tmp, ndotwi);
    diffuseShading(radiance, inc_radiance_cos, sr);
    float pdf = ndotwi / (float)PI;
    vec3_scale(radiance, radiance, 1.0f/pdf);

    MatType mat_type = sr->mat->mat_type;
    if(mat_type == PHONG || mat_type == REFLECTIVE || mat_type == TRANSPARENT)        
    {
        vec3 spec_radiance = {0.0f, 0.0f, 0.0f};
        specularShading(spec_radiance, light_dir, inc_radiance_cos, sr);
        vec3_scale(spec_radiance, spec_radiance, 1.0f/pdf);
        vec3_add(radiance, radiance, spec_radiance);
    }
}

// divide vec3 a by its max component if max component > 1
void maxToOne(vec3 r, const vec3 a)
{
    float max;
    max = (a[0] > a[1]) ? a[0] : a[1];
    max = (max > a[2]) ? max : a[2];
    if(max > 1.0f)
    {
        vec3_scale(r, a, 1.0f/max);
    }
}

void directIllumShading(vec3 radiance, const float ndotwi, const vec3 light_dir, const void* light_ptr,
                  const LightType light_type, const ShadeRec* sr)
{
    switch(light_type)
    {
    case AREALIGHT:
    {
        AreaLight* area_light_ptr = (AreaLight*)light_ptr;
        areaLightShading(radiance, ndotwi, light_dir, area_light_ptr, sr);
    } break;
    case ENVLIGHT:
    {
        EnvLight* env_light_ptr = (EnvLight*)light_ptr;
        envLightShading(radiance, ndotwi, light_dir, env_light_ptr, sr);
    } break;
    default:
        // Diffuse component
        // kd*cd/PI * inc_radiance_cos
        vec3 inc_radiance_cos, tmp;
        getIncRadiance(tmp, light_type, light_ptr);
        vec3_scale(inc_radiance_cos, tmp, ndotwi);
        diffuseShading(radiance, inc_radiance_cos,  sr);

        MatType mat_type = sr->mat->mat_type;
        if(mat_type == PHONG || mat_type == REFLECTIVE || mat_type == TRANSPARENT)                    
        {
            // Specular component
            specularShading(radiance, light_dir, inc_radiance_cos, sr);
        }
    }
}

bool shadowTest(const int light_index, const SceneLights* sl, const SceneObjects* so,
                const vec3 light_dir, const ShadeRec* sr)
{
    if(sl->shadow[light_index] && sr->mat->shadow)
    {
        Ray shadow_ray;
        vec3_copy(shadow_ray.origin, sr->hit_point);
        vec3_copy(shadow_ray.direction, light_dir);
        float min_t = shadowIntersectTest(so, shadow_ray);     
        float t = calcLightDistance(sl->light_types[light_index],
                                    sl->light_ptrs[light_index], sr->hit_point);
        if(min_t < t)
        {
            return true;
        }
        return false;
    }
    return false;
}
