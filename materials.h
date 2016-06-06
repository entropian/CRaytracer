#pragma once

#include "util/vec.h"

enum MatType
{
    MATTE,
    PHONG,
    EMISSIVE
};
// NOTE: this looks like it sucks
typedef struct 
{
    bool shadow;
    MatType mat_type;
    float ka, kd, ks, ke, exp;
    vec3 ca, cd, cs, ce;
}Material;

void printMaterial(const Material* mat)
{
    printf("shadow %s\n", mat->shadow ? "true" : "false");
    printf("ka %f, kd %f\n", mat->ka, mat->kd);
}

void initDefaultMatteMat(Material* mat, const vec3 color)
{
    vec3_copy(mat->cd, color);
    vec3_copy(mat->ca, color);
    mat->kd = 0.6f;
    mat->ka = 1.0f;
    mat->mat_type  = MATTE;
    mat->shadow = true;    
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
}
