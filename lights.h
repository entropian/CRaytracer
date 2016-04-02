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
    Samples2D* samples;
    float intensity;
    float inverse_area;
    vec3 sample_point;
    vec3 color;
} AreaLight;

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

void getLightDir(vec3 r, const LightType light_type, const void* light_ptr, const vec3 hit_point)
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
        vec3_sub(displacement, ((PointLight*)light_ptr)->point, hit_point);
        vec3_normalize(r, displacement);        
    } break;
    case AREALIGHT:
    {
        // Side effect: calculates and stores the surface sample point 
        AreaLight* area_light_ptr = (AreaLight*)light_ptr;
        vec2 sample;
        getNextSample2D(sample, area_light_ptr->samples);
        vec3 displacement;
        // Assuming area light is rectangle
        Rectangle* rect = (Rectangle*)(area_light_ptr->obj_ptr);
        vec3_scale(displacement, rect->width, sample[0]);
        vec3_add(area_light_ptr->sample_point, rect->point, displacement);
        vec3_scale(displacement, rect->height, sample[1]);
        vec3_add(area_light_ptr->sample_point, area_light_ptr->sample_point, displacement);
        vec3_sub(displacement, area_light_ptr->sample_point, hit_point);
        vec3_normalize(r, displacement);
    }  break;
    }
}

void getIncRadiance(vec3 r, const LightType light_type, const void* light_ptr)
{
    switch(light_type)
    {
    case DIRECTIONAL:
        vec3_scale(r, ((DirLight*)light_ptr)->color, ((DirLight*)light_ptr)->intensity);
        break;
    case POINTLIGHT:
        vec3_scale(r, ((PointLight*)light_ptr)->color, ((PointLight*)light_ptr)->intensity);
        break;
    }
}
