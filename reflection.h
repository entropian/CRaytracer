#pragma once
#include "util/vec.h"
#include "mempool.h"
#include <stdio.h>
#include "util/util.h"
#include "microfacet.h"

enum BxDFFlags
{
    BSDF_NONE = 0,
    BSDF_REFLECTION = 1,
    BSDF_TRANSMISSION = 2,
    BSDF_DIFFUSE = 4,
    BSDF_GLOSSY = 8,
    BSDF_SPECULAR = 16,
    BSDF_ALL = BSDF_REFLECTION | BSDF_TRANSMISSION | BSDF_DIFFUSE | BSDF_GLOSSY | BSDF_SPECULAR,
};

enum BxDFType
{
    LAMBERTIAN,
    ORENNAYAR,
    SPECULAR_REFLECTION,
    SPECULAR_TRANSMISSION,
    MICROFACET_REFLECTION,
    MICROFACET_FRESNEL,
    FRESNEL_BLEND,
    FRESNEL_BLEND_DIFFUSE,
    FRESNEL_BLEND_SPECULAR
};

inline BxDFFlags getBxDFFlagsFromType(const BxDFType type)
{
    switch(type)
    {
    case LAMBERTIAN:
    case ORENNAYAR:
        return (BxDFFlags)(BSDF_REFLECTION | BSDF_DIFFUSE);
    case SPECULAR_REFLECTION:
        return (BxDFFlags)(BSDF_REFLECTION | BSDF_SPECULAR);
    case SPECULAR_TRANSMISSION:
        return (BxDFFlags)(BSDF_REFLECTION | BSDF_TRANSMISSION | BSDF_SPECULAR);
    case MICROFACET_REFLECTION:
        return (BxDFFlags)(BSDF_REFLECTION | BSDF_GLOSSY);
    case MICROFACET_FRESNEL:
        return (BxDFFlags)(BSDF_REFLECTION | BSDF_TRANSMISSION | BSDF_GLOSSY);
    case FRESNEL_BLEND:
        //return (BxDFFlags)(BSDF_REFLECTION | BSDF_DIFFUSE | BSDF_GLOSSY);
        return (BxDFFlags)(BSDF_REFLECTION | BSDF_DIFFUSE);
    case FRESNEL_BLEND_DIFFUSE:
        return (BxDFFlags)(BSDF_REFLECTION | BSDF_DIFFUSE);
    case FRESNEL_BLEND_SPECULAR:
        return (BxDFFlags)(BSDF_REFLECTION | BSDF_GLOSSY);        
    default:
        return BSDF_NONE;
    }
}

typedef struct Lambertian_s
{
    vec3 cd;
}Lambertian;
void Lambertian_f(vec3 f, const vec3 wi, const vec3 wo, const Lambertian* l);
float Lambertian_pdf(const vec3 wi, const vec3 wo);
float Lambertian_sample_f(vec3 f, vec3 wi,
                          const vec3 wo, const vec2 sample, const Lambertian* l);

typedef struct OrenNayar_s
{
    vec3 r;
    float a, b;
}OrenNayar;
void OrenNayar_f(vec3 f, const vec3 wi, const vec3 wo, const OrenNayar* on);
float OrenNayar_pdf(const vec3 wi, const vec3 wo);
float OrenNayar_sample_f(vec3 f, vec3 wi, const vec3 wo, const vec2 sample, const OrenNayar* on);


typedef struct SpecularReflection_s
{
    vec3 cr;
    // TODO: fresnel
}SpecularReflection;
float SpecularReflection_sample_f(vec3 f, vec3 wi,
                                  const vec3 wo, const vec2 sample, const SpecularReflection* spec_ref);

typedef struct SpecularTransmission_s
{
    float ior_in, ior_out;
    vec3 cf_in, cf_out;
}SpecularTransmission;
float SpecularTransmission_sample_f(vec3 f, vec3 wi,
                                    const vec3 wo, const vec2 sample, const SpecularTransmission* spec_trans);
void BxDF_f(vec3 f, const vec3 wi, const vec3 wo, const void* bxdf, const BxDFType type);
float BxDF_sample_f(vec3 f, vec3 wi, const vec3 wo, const vec2 sample, const void* bxdf, const BxDFType type);
float BxDF_pdf(const vec3 wi, const vec3 wo, const void* bxdf, const BxDFType type);


typedef struct MicrofacetReflection_s
{
    vec3 color;
    float ior_in, ior_out;

    bool is_metal;
    vec3 k, etaI, etaT; // Used for fresnel if material is metallic
    
    MicrofacetDistribution distrib;
    // TODO: add enum to signify whether to use fresnel dielectric or conductor
}MicrofacetReflection;
void MicrofacetReflection_f(vec3 f, const vec3 wi, const vec3 wo, const MicrofacetReflection* mr);
float MicrofacetReflection_sample_f(vec3 f, vec3 wi, const vec3 wo, const vec2 sample,
                                    const MicrofacetReflection* mr);
float MicrofacetReflection_pdf(const vec3 wi, const vec3 wo, const MicrofacetReflection* mr);

typedef struct MicrofacetFresnel_s
{
    vec3 color;
    float ior_in, ior_out;
    MicrofacetDistribution distrib;
}MicrofacetFresnel;
void MicrofacetFresnel_f(vec3 f, const vec3 wi, const vec3 wo, const MicrofacetFresnel * mt);
float MicrofacetFresnel_sample_f(vec3 f, vec3 wi, const vec3 wo, const vec2 sample,
                                      const MicrofacetFresnel mt);
float MicrofacetFresnel_pdf(const vec3 wi, const vec3 wo, const MicrofacetFresnel* mt);

typedef struct FresnelBlend_s
{
    vec3 kd, ks;
    float ior_in, ior_out;
    MicrofacetDistribution distrib;
}FresnelBlend;
void FresnelBlend_f(vec3 f, const vec3 wi, const vec3 wo, const FresnelBlend * fb);
float FresnelBlend_sample_f(vec3 f, vec3 wi, const vec3 wo, const vec2 sample,
                                      const FresnelBlend mt);
float FresnelBlend_pdf(const vec3 wi, const vec3 wo, const FresnelBlend* mt);

typedef struct FresnelBlendDiffuse_s
{
    vec3 kd, ks;
}FresnelBlendDiffuse;
void FresnelBlendDiffuse_f(vec3 f, const vec3 wi, const vec3 wo, const FresnelBlendDiffuse * fbd);
float FresnelBlendDiffuse_sample_f(vec3 f, vec3 wi, const vec3 wo, const vec2 sample,
                                   const FresnelBlendDiffuse* fbd);
float FresnelBlendDiffuse_pdf(const vec3 wi, const vec3 wo, const FresnelBlendDiffuse* fbd);

typedef struct FresnelBlendSpecular_s
{
    vec3 ks;
    float ior_in, ior_out;
    MicrofacetDistribution distrib;
}FresnelBlendSpecular;
void FresnelBlendSpecular_f(vec3 f, const vec3 wi, const vec3 wo, const FresnelBlendSpecular * fbs);
float FresnelBlendSpecular_sample_f(vec3 f, vec3 wi, const vec3 wo, const vec2 sample,
                                    const FresnelBlendSpecular* fbs);
float FresnelBlendSpecular_pdf(const vec3 wi, const vec3 wo, const FresnelBlendSpecular* fbs);

const int MAX_BXDF = 8;
typedef struct BSDF_s
{
    void* bxdfs[MAX_BXDF];
    BxDFType types[MAX_BXDF];
    int num_bxdf;
    vec3 normal;
    vec3 tangent;
    vec3 binormal;
}BSDF;

bool initBSDFMem(const int num_threads, const int num_depth);
void freeBSDFMem();
void* allocateBxDF();

void freeBxDF(void** bxdf);

void BSDF_f(vec3 f, const vec3 wi, const vec3 wo, const BSDF* bsdf, const BxDFFlags excluded_bxdf);
float BSDF_sample_f(vec3 f, vec3 wi, BxDFFlags* sampled_flags,
                    const vec3 wo, const vec2 sample, const BSDF* bsdf);
float BSDF_pdf(const vec3 wi, const vec3 wo, const BSDF* bsdf);
void BSDF_addLambertian(BSDF* bsdf, const vec3 cd);
void BSDF_addOrenNayar(BSDF* bsdf, const vec3 r, const float sigma);
void BSDF_addSpecularReflection(BSDF* bsdf, const vec3 cr);
void BSDF_addSpecularTransmission(BSDF* bsdf, const float ior_in, const float ior_out, const vec3 cf_in,
                                  const vec3 cf_out);
void BSDF_addMicrofacetReflection(BSDF* bsdf, const vec3 color, const float ior_in, const float ior_out,
                                  const float alphax, const float alphay, const FacetDistribType type);
void BSDF_addMicrofacetReflectionMetal(BSDF* bsdf, const vec3 color, const vec3 etaT, const vec3 etaI, const vec3 k,
                                       const float alphax, const float alphay, const FacetDistribType type);
void BSDF_addMicrofacetFresnel(BSDF* bsdf, const vec3 color, const float ior_in, const float ior_out,
                                    const float alphax, const float alphay, const FacetDistribType type);
void BSDF_addFresnelBlend(BSDF* bsdf, const vec3 kd, const vec3 ks, const float ior_in, const float ior_out,
                          const float alphax, const float alphay, const FacetDistribType type);
void BSDF_addFresnelBlendDiffuse(BSDF* bsdf, const vec3 kd, const vec3 ks);
void BSDF_addFresnelBlendSpecular(BSDF* bsdf, const vec3 ks, const float ior_in, const float ior_out,
                                  const float alphax, const float alphay, const FacetDistribType type);
void BSDF_freeBxDFs(BSDF* bsdf);

