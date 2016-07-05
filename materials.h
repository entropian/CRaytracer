#pragma once

#include "util/vec.h"
#include "sampling.h"

enum MatType
{
    MATTE,
    PHONG,
    REFLECTIVE,
    TRANSPARENT,
    EMISSIVE
};

// NOTE: this looks like it sucks
typedef struct Material_s
{
    bool shadow;
    MatType mat_type;
    float ka, kd, ks, ke, kr, exp;
    float ior_in, ior_out;    
    vec3 ca, cd, cs, ce, cr;
    vec3 cf_in, cf_out;
    Samples3D* h_samples;
}Material;

void printMaterial(const Material* mat)
{
    printf("shadow %s\n", mat->shadow ? "true" : "false");
    printf("ka %f, kd %f\n", mat->ka, mat->kd);
    printVec3WithText("cd", mat->cd);
    printVec3WithText("cs", mat->cs);
    printf("ior_in %f, ior_out %f\n", mat->ior_in, mat->ior_out);
    printVec3WithText("cf_in", mat->cf_in);
    printVec3WithText("cf_out", mat->cf_out);    
}

void initDefaultMatteMat(Material* mat, const vec3 color)
{
    vec3_copy(mat->cd, color);
    vec3_copy(mat->ca, color);
    mat->kd = 0.6f;
    mat->ka = 1.0f;
    mat->mat_type  = MATTE;
    mat->shadow = true;
    mat->h_samples = genHemisphereSamples(MULTIJITTERED, 1.0f);    
}       

void initDefaultPhongMat(Material* mat, const vec3 color)
{
    vec3_copy(mat->cd, color);
    vec3_copy(mat->cs, color);
    vec3_copy(mat->ca, color);    
    mat->kd = 0.6f;
    mat->ka = 1.0f;
    mat->ks = 0.4f;
    mat->exp = 10.0f;            
    mat->mat_type = PHONG;
    mat->shadow = true;
    mat->h_samples = genHemisphereSamples(MULTIJITTERED, DEFAULT_GLOSSINESS);
}
