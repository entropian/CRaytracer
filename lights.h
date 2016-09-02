#pragma once

#include "util/vec.h"
#include "util/math.h"
#include "shapes/shapes.h"
#include "sampling.h"

enum LightType
{
    DIRECTIONAL,
    POINTLIGHT,
    AREALIGHT,
    ENVLIGHT
};

typedef struct DirLight_s
{
    float intensity;
    vec3 color;
    vec3 direction;
}DirLight;

typedef struct PointLight_s
{
    float intensity;
    vec3 color;
    vec3 point;
}PointLight;

typedef struct AreaLight_s
{
    ObjectType obj_type; 
    void* obj_ptr;
    float intensity;
    float pdf;
    Samples2D* samples2D;
    Samples3D* samples3D;
    vec3 sample_point;
    vec3 color;
}AreaLight;

typedef struct EnvLight_s
{
    float intensity;
    Samples3D* samples3D;
    vec3 color;
}EnvLight;

typedef struct AmbientLight_s
{
    bool amb_occlusion;
    float intensity;
    vec3 color;
}AmbientLight;

void getAreaLightNormal(vec3 r, const AreaLight* area_light_ptr, const vec3 hit_point);
void assignDirLight(DirLight *dir_light, const float intensity, const vec3 color, const vec3 direction);
void assignPointLight(PointLight *point_light, const float intensity, const vec3 color, const vec3 point);

// Calculate and store light direction vector in r
void getLightDir(vec3 r, const LightType light_type, const void* light_ptr, const ShadeRec* sr, const int);

// Calculate and store incident radiance in r
void getIncRadiance(vec3 r, const LightType light_type, const void* light_ptr);

float calcLightDistance(const LightType light_type, const void* light_ptr, const vec3 hit_point);

