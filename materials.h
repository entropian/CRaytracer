#pragma once

#include "util/vec.h"
#include "sampling.h"
#include "texture.h"
#include "reflection.h"

enum MatType
{
    INVALID_MAT_TYPE = 0,
    MATTE = 1,
    PHONG = 2,
    REFLECTIVE = 4,
    TRANSPARENT = 8,
    EMISSIVE = 16,
    PARTICIPATING = 32
};

enum TextureType
{
    NO_TEXTURE = 0,
    DIFFUSE = 1,
    NORMAL = 2,
    SPECULAR = 4,
    NOISE = 8
};

typedef struct Material_s
{
    void* data;
    char name[MAX_NAME_LENGTH];
    MatType mat_type;
}Material;

typedef struct Matte_s
{
    vec3 color;
    float sigma;
    Texture* diffuse;   // Diffuse and normal maps
    Texture* normal;
    char diffuse_file_name[MAX_NAME_LENGTH];
    char normal_file_name[MAX_NAME_LENGTH];
}Matte;

typedef struct Reflective_s
{
    vec3 color;
}Reflective;

typedef struct Transparent_s
{
    float ior_in, ior_out;
    vec3 cf_in, cf_out;
}Transparent;

typedef struct Emissive_s
{
    vec3 color;
    float intensity;
}Emissive;

const unsigned int DIFFUSE_MAP_INDEX = 0;
const unsigned int NORMAL_MAP_INDEX = 1;
const unsigned int SPEC_MAP_INDEX = 2;

MatType getMatTypeFromString(const char* str);
void getMaterialDiffuseTexColor(vec3 texel, const Material *mat, const vec2 uv);
void getMaterialNormalTexColor(vec3 texel, const Material *mat, const vec2 uv);
void getMaterialSpecularTexColor(vec3 texel, const Material *mat, const vec2 uv);
void setMaterialDiffuseTexPtr(Material *mat, Texture *tex);
void setMaterialNormalTexPtr(Material *mat, Texture *tex);
void setMaterialSpecularTexPtr(Material *mat, Texture *tex);
void computeScatteringFunc(BSDF* bsdf, const vec2 uv, const Material* mat);
bool Material_hasNormalMap(const Material* mat);
void Material_getNormalMapValue(vec3 value, const Material* mat, const vec2 uv);
