#pragma once

enum MatType
{
    MATTE,
    PHONG
};

typedef struct Material
{
    bool shadow;
    MatType mat_type;
    float ka, kd, ks, exp;
    vec3 ca, cd, cs;
}Material;
