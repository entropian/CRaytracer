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
    MICROFACET_TRANSMISSION
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
        return (BxDFFlags)(BSDF_TRANSMISSION | BSDF_SPECULAR);
    case MICROFACET_REFLECTION:
        return (BxDFFlags)(BSDF_REFLECTION | BSDF_GLOSSY);
    case MICROFACET_TRANSMISSION:
        return (BxDFFlags)(BSDF_TRANSMISSION | BSDF_GLOSSY);
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
    MicrofacetDistribution distrib;
    // TODO: add enum to signify whether to use fresnel dielectric or conductor
}MicrofacetReflection;
void MicrofacetReflection_f(vec3 f, const vec3 wi, const vec3 wo, const MicrofacetReflection* mr);
float MicrofacetReflection_sample_f(vec3 f, vec3 wi, const vec3 wo, const vec2 sample,
                                    const MicrofacetReflection* mr);
float MicrofacetReflection_pdf(const vec3 wi, const vec3 wo, const MicrofacetReflection* mr);

typedef struct MicrofacetTransmission_s
{
    vec3 color;
    float ior_in, ior_out;
    MicrofacetDistribution distrib;
}MicrofacetTransmission;
void MicrofacetTransmission_f(vec3 f, const vec3 wi, const vec3 wo, const MicrofacetTransmission * mt);
float MicrofacetTransmission_sample_f(vec3 f, vec3 wi, const vec3 wo, const vec2 sample,
                                      const MicrofacetTransmission mt);
float MicrofacetTransmission_pdf(const vec3 wi, const vec3 wo, const MicrofacetTransmission* mt);

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
void BSDF_addMicrofacetTransmission(BSDF* bsdf, const vec3 color, const float ior_in, const float ior_out,
                                    const float alphax, const float alphay, const FacetDistribType type);
void BSDF_freeBxDFs(BSDF* bsdf);

