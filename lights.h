#pragma once

#include "vec.h"

enum LightType
{
    DIRECTIONAL,
    POINTLIGHT,
    AREALIGHT
};

typedef struct DirLight
{
    float intensity;
    vec3 color;
    vec3 direction;
} DirLight;

typedef struct PointLight
{
    float intensity;
    vec3 color;
    vec3 point;
} PointLight;

typedef struct AreaLight
{
    ObjectType obj_type;
    void* obj_ptr;
    Samples2D* samples2D;
    Samples3D* samples3D;
    float intensity;
    float pdf;
    vec3 sample_point;
    vec3 color;
} AreaLight;

void getAreaLightNormal(vec3 r, const AreaLight* area_light_ptr, const vec3 hit_point)
{
    switch(area_light_ptr->obj_type)
    {
    case RECTANGLE:
    {
        vec3_copy(r, ((Rectangle*)(area_light_ptr->obj_ptr))->normal);
    } break;
    case SPHERE:
    {
        vec3 displacement;
        vec3_sub(displacement, hit_point, ((Sphere*)(area_light_ptr->obj_ptr))->center);
        vec3_normalize(r, displacement);
    } break;
    }
}

void assignDirLight(DirLight *dir_light, const float intensity, const vec3 color, const vec3 direction)
{
    dir_light->intensity = intensity;
    vec3_copy(dir_light->color, color);
    vec3_copy(dir_light->direction, direction);
}

void assignPointLight(PointLight *point_light, const float intensity, const vec3 color, const vec3 point)
{
    point_light->intensity = intensity;
    vec3_copy(point_light->color, color);
    vec3_copy(point_light->point, point);
}

// Calculate and store light direction vector in r
void getLightDir(vec3 r, const LightType light_type, const void* light_ptr, const ShadeRec* sr)
{
    switch(light_type)
    {
    case DIRECTIONAL:
    {
        vec3_copy(r, ((DirLight*)light_ptr)->direction);
    } break;
    case POINTLIGHT:
    {
        vec3 displacement;
        vec3_sub(displacement, ((PointLight*)light_ptr)->point, sr->hit_point);
        vec3_normalize(r, displacement);        
    } break;
    case AREALIGHT:
    {
        // Side effect: calculates and stores the surface sample point 
        AreaLight* area_light_ptr = (AreaLight*)light_ptr;
        switch(area_light_ptr->obj_type)
        {
        case RECTANGLE:
        {
            vec2 sample;
            getNextSample2D(sample, area_light_ptr->samples2D);
            vec3 displacement;
            Rectangle* rect = (Rectangle*)(area_light_ptr->obj_ptr);
            vec3_scale(displacement, rect->width, sample[0]);
            vec3_add(area_light_ptr->sample_point, rect->point, displacement);
            vec3_scale(displacement, rect->height, sample[1]);
            vec3_add(area_light_ptr->sample_point, area_light_ptr->sample_point, displacement);            
            vec3_sub(displacement, area_light_ptr->sample_point, sr->hit_point);
            vec3_normalize(r, displacement);
        } break;
        case SPHERE:
        {
            vec3 h_sample;
            getNextSample3D(h_sample, area_light_ptr->samples3D);
            vec3 u, v, w;
            vec3_copy(w, sr->normal);
            vec3_cross(v, w, JITTERED_UP);
            vec3_normalize(v, v);
            vec3_cross(u, v, w);
            orthoNormalTransform(area_light_ptr->sample_point, u, v, w, h_sample);
            vec3_add(area_light_ptr->sample_point, area_light_ptr->sample_point,
                     ((Sphere*)(area_light_ptr->obj_ptr))->center);
            vec3 displacement;
            vec3_sub(displacement, area_light_ptr->sample_point, sr->hit_point);
            vec3_normalize(r, displacement);
        } break;
        }
    }  break;
    }
}

// Calculate and store incident radiance in r
void getIncRadiance(vec3 r, const LightType light_type, const void* light_ptr)
{
    switch(light_type)
    {
    case DIRECTIONAL:
    {
        vec3_scale(r, ((DirLight*)light_ptr)->color, ((DirLight*)light_ptr)->intensity);
    } break;
    case POINTLIGHT:
    {
        vec3_scale(r, ((PointLight*)light_ptr)->color, ((PointLight*)light_ptr)->intensity);
    } break;
    }
}

float calcLightDistance(const LightType light_type, const void* light_ptr, const vec3 hit_point)
{
    float t;
    switch(light_type)
    {
    case DIRECTIONAL:
    {
        t = TMAX;                                                                    
    } break;
    case POINTLIGHT:
    {
        vec3 light_to_hit_point;
        PointLight* point_light_ptr = (PointLight*)(light_ptr);
        vec3_sub(light_to_hit_point, point_light_ptr->point, hit_point);
        t = vec3_length(light_to_hit_point);
    }break;
    case AREALIGHT:
    {
        vec3 light_to_hit_point;
        AreaLight* area_light_ptr = (AreaLight*)(light_ptr);
        vec3_sub(light_to_hit_point, area_light_ptr->sample_point, hit_point);
        t = vec3_length(light_to_hit_point);
    } break;
    }
    return t;
}

