#pragma once
#include "util/vec.h"
#include "mempool.h"
#include <stdio.h>
#include "util/util.h"

inline float cosTheta(const vec3 w) { return w[2]; }
inline float cos2Theta(const vec3 w) { return w[2] * w[2]; }
inline float absCosTheta(const vec3 w) { return fabs(w[2]); }
inline float sin2Theta(const vec3 w) {
    return max((float)0, (float)1 - cos2Theta(w));
}

inline float sinTheta(const vec3 w) { return sqrt(sin2Theta(w)); }

inline float tanTheta(const vec3 w) { return sinTheta(w) / cosTheta(w); }

inline float tan2Theta(const vec3 w) {
    return sin2Theta(w) / cos2Theta(w);
}

inline float cosPhi(const vec3 w) {
    float sin_theta = sinTheta(w);
    return (sin_theta == 0) ? 1 : clamp(w[0] / sin_theta, -1.0f, 1.0f);
}

inline float sinPhi(const vec3 w) {
    float sin_theta = sinTheta(w);
    return (sin_theta == 0) ? 0 : clamp(w[1] / sin_theta, -1.0f, 1.0f);
}

inline float cos2Phi(const vec3 w) { return cosPhi(w) * cosPhi(w); }

inline float sin2Phi(const vec3 w) { return sinPhi(w) * sinPhi(w); }

inline float cosDPhi(const vec3 wa, const vec3 wb) {
    return clamp(
        (wa[0] * wb[0] + wa[1] * wb[1]) / sqrt((wa[0] * wa[0] + wa[1] * wa[1]) *
                                                (wb[0] * wb[0] + wb[1] * wb[1])),
        -1, 1);
}

enum BxDFType
{
    LAMBERTIAN,
    ORENNAYAR,
    SPECULAR_REFLECTION,
    SPECULAR_TRANSMISSION,
};

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
// TODO
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

void BSDF_f(vec3 f, const vec3 wi, const vec3 wo, const BSDF* bsdf);
float BSDF_sample_f(vec3 f, vec3 wi,
                    const vec3 wo, const vec2 sample, const BSDF* bsdf);
float BSDF_pdf(const vec3 wi, const vec3 wo, const BSDF* bsdf);
void BSDF_addLambertian(BSDF* bsdf, const vec3 cd);
void BSDF_addOrenNayar(BSDF* bsdf, const vec3 r, const float sigma);
void BSDF_addSpecularReflection(BSDF* bsdf, const vec3 cr);
void BSDF_addSpecularTransmission(BSDF* bsdf, const float ior_in, const float ior_out, const vec3 cf_in,
                                  const vec3 cf_out);
void BSDF_freeBxDFs(BSDF* bsdf);

