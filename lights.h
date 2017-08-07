#pragma once

#include "util/vec.h"
#include "util/math.h"
#include "shapes/shapes.h"
#include "sampling.h"
#include "texture.h"

enum LightType
{
    DIRECTIONAL,
    POINTLIGHT,
    AREALIGHT,
    ENVLIGHT,
    MESHLIGHT
};

typedef struct DirLight_s
{
    float intensity;
    vec3 color;
    vec3 direction;
}DirLight;

typedef struct PointLight_s
{
    bool dist_atten;
    float intensity;
    float flux;
    vec3 color;
    vec3 point;
    char *proj_map;
    float proj_coverage;
}PointLight;

typedef struct AreaLight_s
{
    ObjectType obj_type;
    void* obj_ptr;
    float intensity;
    float flux;
    float pdf;
    Samples2D* samples2D;
    Samples3D* samples3D;
    vec3 sample_point;
    vec3 color;
}AreaLight;

enum EnvLightType
{
    CONSTANT,
    CUBEMAP
};

typedef struct EnvLight_s
{
    float intensity;
    Samples3D* samples3D;
    vec3 color;
    /*
      0 -z
      1 z
      2 -x
      3 x
      4 -y
      5 y
     */
    Texture cubemap[6];
    EnvLightType type;
}EnvLight;

typedef struct AmbientLight_s
{
    bool amb_occlusion;
    float intensity;
    vec3 color;
}AmbientLight;

typedef struct MeshLight_s
{
    float power;
    float surface_area;
    float intensity;
    vec3 color;
    int num_triangles;
    float* probability_distribution; // TODO rename
    void** triangles;
    ObjectType obj_type;
}MeshLight;

void genMeshLightSample(vec3 sample, vec3 normal, MeshLight* mesh_light);
void getAreaLightNormal(vec3 r, const AreaLight* area_light_ptr, const vec3 hit_point);
void assignDirLight(DirLight *dir_light, const float intensity, const vec3 color, const vec3 direction);
void assignPointLight(PointLight *point_light, const float intensity, const vec3 color, const vec3 point);

// Calculate and store light direction vector in r
void getLightDir(vec3 r, const LightType light_type, const void* light_ptr, const ShadeRec* sr, const int);

// Calculate and store incident radiance in r
void getIncRadiance(vec3 r, const LightType light_type, const void* light_ptr, const vec3 hit_point);

float calcLightDistance(const LightType light_type, const void* light_ptr, const vec3 hit_point);
void getEnvLightIncRadiance(vec3 r, const vec3 dir, EnvLight* env_light);
void EnvLight_init_cubemap(EnvLight* env_light, char paths[][256]);
void EnvLight_destroy(EnvLight* env_light);
