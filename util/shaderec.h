#pragma once

#include "vec.h"
#include "../materials.h"
#include "math.h"

typedef struct ShadeRec
{
    bool hit_status;
    vec3 hit_point;
    vec3 normal;
    vec3 shading_normal;
    vec3 dpdu;
    vec3 dpdv;
    vec3 wo;
    vec2 uv;
    Material *mat;                // Points to the material member of an object    
} ShadeRec;

static void updateShadeRecWithTexInfo(ShadeRec *sr)
{
    if(sr->mat->tex_flags & DIFFUSE)
    {
        getMaterialDiffuseTexColor(sr->mat->cd, sr->mat, sr->uv);
    }
    if(sr->mat->tex_flags & SPECULAR)
    {
        // TEMP:
        vec3 spec_color;
        getMaterialSpecularTexColor(spec_color, sr->mat, sr->uv);
        sr->mat->ks = spec_color[0];
        sr->mat->kd = 1.0f - sr->mat->ks;
    }
}
