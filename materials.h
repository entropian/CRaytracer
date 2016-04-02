#pragma once

enum MatType
{
    MATTE,
    PHONG,
    EMISSIVE
};
// NOTE: this looks like it sucks
typedef struct Material
{
    bool shadow;
    MatType mat_type;
    float ka, kd, ks, ke, exp;
    vec3 ca, cd, cs, ce;
}Material;
