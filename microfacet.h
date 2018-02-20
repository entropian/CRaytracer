#pragma once
#include "util/vec.h"

enum FacetDistribType
{
    BECKMANN,
    TROWBRIDGEREITZ
};

typedef struct MicrofacetDistribution_s
{
    // Distribution could be either Beckmann or TrowbridgeReitz
    // Both distributions use the same data    
    float alphax, alphay;
    FacetDistribType type;
}MicrofacetDistribution;

float MicrofacetDistribution_D(const vec3 wh, const MicrofacetDistribution* distrib);
void MicrofacetDistribution_sample_wh(vec3 wh, const vec3 wo, const vec2 sample, const MicrofacetDistribution* mf);
float MicrofacetDistribution_G1(const vec3 w, const MicrofacetDistribution* distrib);
float MicrofacetDistribution_G(const vec3 wo, const vec3 wi, const MicrofacetDistribution* distrib);
float MicrofacetDistribution_pdf(const vec3 wo, const vec3 wh, const MicrofacetDistribution* distrib);


inline float BeckmannRoughnessToAlpha(float roughness)
{
    roughness = max(roughness, 1e-3);
    float x = logf(roughness);
    return 1.62142f + 0.819955f * x + 0.1734f * x * x +
        0.0171201f * x * x * x + 0.000640711f * x * x * x * x;
}
