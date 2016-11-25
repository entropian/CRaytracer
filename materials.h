#pragma once

#include "util/vec.h"
#include "sampling.h"
#include "texture.h"

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
    bool shadow;
    MatType mat_type;
    unsigned int tex_flags;
    float ka, kd, ks, ke, kr, exp;
    float ior_in, ior_out;    
    vec3 ca, cd, cs, ce, cr;
    vec3 cf_in, cf_out;
    Samples3D* h_samples;
    Texture* tex_array[3];
}Material;

const int DIFFUSE_MAP_INDEX = 0;
const int NORMAL_MAP_INDEX = 1;

MatType getMatTypeFromString(const char* str);
void printMaterial(const Material* mat);
void initMaterial(Material *mat);
void initDefaultMatteMat(Material* mat, const vec3 color);
void initDefaultPhongMat(Material* mat, const vec3 color);
void getMaterialDiffuseTexColor(vec3 texel, const Material *mat, const vec2 uv);
void setMaterialDiffuseTexPtr(Material *mat, Texture *tex);
void setMaterialNormalTexPtr(Material *mat, Texture *tex);
