#pragma once

#include "vec.h"
#include "../materials.h"
#include "math.h"

typedef struct ShadeRec
{
    bool hit_status;
    vec3 hit_point;
    vec3 normal;
    vec3 wo;
    vec2 uv;
    Material *mat;                // Points to the material member of an object    
} ShadeRec;

static void updateShadeRecWithTexInfo(ShadeRec *sr)
{
    // Update diffuse color with diffuse texture value
    if(sr->mat->tex_flags & DIFFUSE)
    {
        getMaterialDiffuseTexColor(sr->mat->cd, sr->mat, sr->uv);
    }
    /*
    // Update normal with normal texture value
    if(sr->mat->tex_flags & NORMAL)
    {
        vec3 texel;
        getMaterialNormalTexColor(texel, sr->mat, sr->uv);
        getVec3InLocalBasis(sr->normal, texel, sr->normal);
        vec3_normalize(sr->normal, sr->normal);
    }
    */
}
