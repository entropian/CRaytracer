#pragma once

#include "vec.h"
#include "../materials.h"

typedef struct ShadeRec
{
    bool hit_status;
    vec3 hit_point;
    vec3 normal;
    vec3 wo;
    vec2 uv;
    Material *mat;                // Points to the material member of an object    
} ShadeRec;
