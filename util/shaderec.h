#pragma once

#include "vec.h"
#include "../materials.h"
#include "math.h"

typedef struct ShadeRec
{
    bool hit_status;
    vec3 hit_point;
    vec3 normal;
    //vec3 shading_normal;
    vec3 dpdu;
    vec3 dpdv;
    vec3 wo;
    vec2 uv;
    Material mat;
    BSDF bsdf;
} ShadeRec;


static void updateShadeRecWithTexInfo(ShadeRec *sr)
{
    if(sr->mat.tex_flags & DIFFUSE)
    {
        getMaterialDiffuseTexColor(sr->mat.cd, &(sr->mat), sr->uv);
    }
    /*
    if(sr->mat.tex_flags & NORMAL)
    {
        vec3 tex_normal, normal, surface_normal;
        vec3_copy(surface_normal, sr->normal);
        getMaterialNormalTexColor(tex_normal, &(sr->mat), sr->uv);
        orthoNormalTransform(normal, sr->dpdu, sr->dpdv, sr->normal, tex_normal);
        vec3_normalize(normal, normal);
        vec3 displacement;
        vec3_sub(displacement, normal, sr->normal);
        if(vec3_equal(sr->normal, BLACK))
        {
            vec3_copy(sr->normal, surface_normal);
        }
    }
    */
    if(sr->mat.tex_flags & SPECULAR)
    {
        // TEMP:
        vec3 spec_color;
        getMaterialSpecularTexColor(spec_color, &(sr->mat), sr->uv);
        sr->mat.ks = spec_color[0];
        sr->mat.kd = 1.0f - sr->mat.ks;
    }
}
