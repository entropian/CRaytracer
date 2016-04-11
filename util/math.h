#pragma once

#include "vec.h"

void orthoNormalTransform(vec3 r, const vec3 u, const vec3 v, const vec3 w, const vec3 a)
{
    vec3 u_comp, v_comp, w_comp;
    vec3_scale(u_comp, u, a[0]);
    vec3_scale(v_comp, v, a[1]);
    vec3_scale(w_comp, w, a[2]);
    vec3_add(r, u_comp, v_comp);
    vec3_add(r, r, w_comp);    
}

void getVec3InLocalBasis(vec3 r, const vec3 a, const vec3 normal)
{
    vec3 u, v, w;
    vec3_copy(w, normal);
    vec3_cross(v, w, JITTERED_UP);
    vec3_normalize(v, v);
    vec3_cross(u, v, w);
    orthoNormalTransform(r, u, v, w, a);
}
