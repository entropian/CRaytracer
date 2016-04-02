#pragma once

enum MatType
{
    MATTE,
    PHONG,
    EMISSIVE
};

typedef struct Material
{
    bool shadow;
    MatType mat_type;
    float ka, kd, ks, ke, exp;
    vec3 ca, cd, cs, ce;
}Material;
