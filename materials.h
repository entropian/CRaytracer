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

typedef struct MaterialNew_s
{
    void* data;
    MatType mat_type;
}MaterialNew;

typedef struct Matte_s
{
    vec3 color;
    float sigma;
    Texture* tex_array[2];   // Diffuse and normal maps
}Matte;

typedef struct Reflective_s
{
    vec3 color;
}Reflective;

typedef struct Transparent_s
{
    float ior_in, ior_out;
}Transparent;

typedef struct Material_s
{
    MatType mat_type;
    unsigned int tex_flags;
    float ka, kd, ks, ke, kr, exp;
    float ior_in, ior_out;
    vec3 ca, cd, cs, ce, cr;
    vec3 cf_in, cf_out;
    float extinct_coeff, scatter_coeff;
    Texture* tex_array[4];
    char name[MAX_NAME_LENGTH];
    Samples3D* h_samples;
}Material;

const unsigned int DIFFUSE_MAP_INDEX = 0;
const unsigned int NORMAL_MAP_INDEX = 1;
const unsigned int SPEC_MAP_INDEX = 2;

MatType getMatTypeFromString(const char* str);
void printMaterial(const Material* mat);
void initMaterial(Material *mat);
void initDefaultMatteMat(Material* mat, const vec3 color);
void initDefaultPhongMat(Material* mat, const vec3 color);
void getMaterialDiffuseTexColor(vec3 texel, const Material *mat, const vec2 uv);
void getMaterialNormalTexColor(vec3 texel, const Material *mat, const vec2 uv);
void getMaterialSpecularTexColor(vec3 texel, const Material *mat, const vec2 uv);
void setMaterialDiffuseTexPtr(Material *mat, Texture *tex);
void setMaterialNormalTexPtr(Material *mat, Texture *tex);
void setMaterialSpecularTexPtr(Material *mat, Texture *tex);
void computeScatteringFunc(BSDF* bsdf, const vec2 uv, const Material* mat);
