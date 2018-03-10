#include "reflection.h"
#include "sampling.h"
#include "util/ray.h"
#include "util/math.h"

float cosHemispherePdf(const vec3 wi, const vec3 wo)
{
    float dot_product = fabs(cosTheta(wi)) * INV_PI;
    if(dot_product == 0.0f)
    {
        //printf("dot product 0\n");
    }    
    //return sameHemisphere(wi, wo) ? fabs(vec3_dot(wi, normal)) * INV_PI: 0.0f;
    return sameHemisphere(wi, wo) ? absCosTheta(wi) * INV_PI: 0.0f;    
}

float cosHemisphereSample(vec3 wi, const vec3 wo, const vec2 sample)
{
    mapSampleToHemisphere(wi, sample);
    return cosHemispherePdf(wi, wo);
}
#if 0
bool refract(vec3 wt, const vec3 n, const vec3 wi, const float eta)
{
/*
    // Compute $\cos \theta_\roman{t}$ using Snell's law
    Float cosThetaI = Dot(n, wi);
    Float sin2ThetaI = std::max(Float(0), Float(1 - cosThetaI * cosThetaI));
    Float sin2ThetaT = eta * eta * sin2ThetaI;

    // Handle total internal reflection for transmission
    if (sin2ThetaT >= 1) return false;
    Float cosThetaT = std::sqrt(1 - sin2ThetaT);
    *wt = eta * -wi + (eta * cosThetaI - cosThetaT) * Vector3f(n);
    return true;
 */
    float cos_theta_i = vec3_dot(n, wi);
    float sin_2_theta_i = max(0.0f, 1.0f - cos_theta_i * cos_theta_i);
    float sin_2_theta_t = eta * eta * sin_2_theta_i;

    if(sin_2_theta_i >= 1.0f) return false;
    float cos_theta_t = sqrtf(1.0f - sin_2_theta_t);
    vec3 scaled_wi, scaled_n;
    vec3_scale(scaled_wi, wi, eta);
    vec3_scale(scaled_n, n, eta * cos_theta_i - cos_theta_t);
    vec3_add(wt, scaled_wi, scaled_n);
    return true;
}
#endif

// Dielectric
float calcFresnelDielectric(const vec3 normal, const vec3 wo, float etaT, float etaI)
{
    vec3 n;
    vec3_copy(n, normal);
    float cos_theta_i = vec3_dot(n, wo);
    if(cos_theta_i < 0.0f)
    {
        float tmp;
        tmp = etaI;
        etaI = etaT;
        etaT = tmp;
        vec3_negate(n, n);
        cos_theta_i = fabs(cos_theta_i);
    }
    float sin_theta_i = sqrtf(max(0.0f, 1.0f - cos_theta_i * cos_theta_i));
    float sin_theta_t = etaI / etaT * sin_theta_i;

    if(sin_theta_t >= 1.0f) return 1.0f;
    float cos_theta_t = sqrtf(max(0.0f, 1- sin_theta_t * sin_theta_t));
    float r_parl = ((etaT * cos_theta_i) - (etaI * cos_theta_t)) /
                   ((etaT * cos_theta_i) + (etaI * cos_theta_t));
    float r_perp = ((etaI * cos_theta_i) - (etaT * cos_theta_t)) /
                   ((etaI * cos_theta_i) + (etaT * cos_theta_t));
    return (r_parl * r_parl + r_perp * r_perp) / 2;
}

void calcFresnelConductor(vec3 fresnel, const vec3 normal, const vec3 wo,
                          const vec3 etaT, const vec3 etaI, const vec3 k)
{
/*
    cosThetaI = Clamp(cosThetaI, -1, 1);
    Spectrum eta = etat / etai;
    Spectrum etak = k / etai;

    Float cosThetaI2 = cosThetaI * cosThetaI;
    Float sinThetaI2 = 1. - cosThetaI2;
    Spectrum eta2 = eta * eta;
    Spectrum etak2 = etak * etak;

    Spectrum t0 = eta2 - etak2 - sinThetaI2;
    Spectrum a2plusb2 = Sqrt(t0 * t0 + 4 * eta2 * etak2);
    Spectrum t1 = a2plusb2 + cosThetaI2;
    Spectrum a = Sqrt(0.5f * (a2plusb2 + t0));
    Spectrum t2 = (Float)2 * cosThetaI * a;
    Spectrum Rs = (t1 - t2) / (t1 + t2);

    Spectrum t3 = cosThetaI2 * a2plusb2 + sinThetaI2 * sinThetaI2;
    Spectrum t4 = t2 * sinThetaI2;
    Spectrum Rp = Rs * (t3 - t4) / (t3 + t4);

    return 0.5 * (Rp + Rs);
 */
    float cos_theta_i = vec3_dot(wo, normal);
    // TEMP
    /*
    cos_theta_i = 0.845039606f;
    vec3 etaI = {1.0f, 1.0f, 1.0f};
    vec3 etaT = {0.14282006, 0.37414363, 1.43944442};
    vec3 k = {3.97471833, 2.38065982, 1.59981036};
    */
    vec3 eta;
    vec3_div(eta, etaT, etaI);
    vec3 etaK;
    vec3_div(etaK, k, etaI);

    float cos_theta_i_2 = cos_theta_i * cos_theta_i;
    float sin_theta_i_2 = 1.0f - cos_theta_i_2;
    vec3 eta2;
    vec3_mult(eta2, eta, eta);
    vec3 etaK2;
    vec3_mult(etaK2, etaK, etaK);

    vec3 t0;
    vec3_sub(t0, eta2, etaK2);
    vec3_sub_c(t0, t0, sin_theta_i_2);
    vec3 a2plusb2, t0_2, eta2_2;
    vec3_mult(t0_2, t0, t0);
    vec3_mult(eta2_2, eta2, etaK2);
    vec3_scale(eta2_2, eta2_2, 4.0f);
    vec3_add(a2plusb2, t0_2, eta2_2);
    vec3_sqrt(a2plusb2, a2plusb2);
    vec3 t1;
    vec3_add_c(t1, a2plusb2, cos_theta_i_2);
    vec3 a;
    vec3_add(a, a2plusb2, t0);
    vec3_scale(a, a, 0.5f);
    vec3_sqrt(a, a);
    vec3 t2;
    vec3_scale(t2, a, 2.0f * cos_theta_i);
    vec3 rs, t1_minus_t2, t1_plus_t2;
    vec3_sub(t1_minus_t2, t1, t2);
    vec3_add(t1_plus_t2, t1, t2);
    vec3_div(rs, t1_minus_t2, t1_plus_t2);
    vec3 t3;
    vec3_scale(t3, a2plusb2, cos_theta_i_2);
    vec3_add_c(t3, t3, sin_theta_i_2 * sin_theta_i_2);
    vec3 t4;
    vec3_scale(t4, t2, sin_theta_i_2);
    vec3 rp, t3_minus_t4, t3_plus_t4;
    vec3_sub(t3_minus_t4, t3, t4);
    vec3_add(t3_plus_t4, t3, t4);
    vec3_div(rp, t3_minus_t4, t3_plus_t4);
    vec3_mult(rp, rp, rs);
    vec3_add(fresnel, rp, rs);
    vec3_scale(fresnel, fresnel, 0.5f);
}

void Lambertian_f(vec3 f, const vec3 wi, const vec3 wo, const Lambertian* l)
{
    vec3_copy(f, l->cd);
    vec3_scale(f, f, INV_PI);
}

float Lambertian_pdf(const vec3 wi, const vec3 wo)
{
    /*
    if(!sameHemisphere(wi, wo))
    {
        printf("not same hemisphere\n");
    }
    */
    //float dot_product = fabs(vec3_dot(wi, normal)) * INV_PI;
    return cosHemispherePdf(wi, wo);
}

float Lambertian_sample_f(vec3 f, vec3 wi,
                         const vec3 wo, const vec2 sample, const Lambertian* l)
{
    vec3 tmp_wo;
    vec3_copy(tmp_wo, wo);
    if(tmp_wo[2] < 0.0f)
    {
        tmp_wo[2] *= -1.0f;
    }
    vec3_copy(f, l->cd);
    vec3_scale(f, f, INV_PI);
    return cosHemisphereSample(wi, tmp_wo, sample);
}

void OrenNayar_f(vec3 f, const vec3 wi, const vec3 wo, const OrenNayar* on)
{
    float sin_theta_i = sinTheta(wi);
    float sin_theta_o = sinTheta(wo);
    // Compute cosine term of Oren-Nayar model
    float max_cos = 0.0f;
    if (sin_theta_i > 1e-4 && sin_theta_o > 1e-4) {
        float sin_phi_i = sinPhi(wi), cos_phi_i = cosPhi(wi);
        float sin_phi_o = sinPhi(wo), cos_phi_o = cosPhi(wo);
        float d_cos = cos_phi_i * cos_phi_o + sin_phi_i * sin_phi_o;
        max_cos = max((float)0, d_cos);
    }

    // Compute sine and tangent terms of Oren-Nayar model
    float sin_alpha, tan_beta;
    if (absCosTheta(wi) > absCosTheta(wo)) {
        sin_alpha = sin_theta_o;
        tan_beta = sin_theta_i / absCosTheta(wi);
    } else {
        sin_alpha = sin_theta_i;
        tan_beta = sin_theta_o / absCosTheta(wo);
    }
    vec3_scale(f, on->r, (on->a + on->b * max_cos * sin_alpha * tan_beta) * INV_PI);
}

float OrenNayar_pdf(const vec3 wi, const vec3 wo)
{
    return cosHemispherePdf(wi, wo);
}

float OrenNayar_sample_f(vec3 f, vec3 wi, const vec3 wo, const vec2 sample, const OrenNayar* on)
{
    vec3 tmp_wo;
    vec3_copy(tmp_wo, wo);
    if(tmp_wo[2] < 0.0f)
    {
        tmp_wo[2] *= -1.0f;
    }
    float pdf = cosHemisphereSample(wi, tmp_wo, sample);
    OrenNayar_f(f, wi, wo, on);
    return pdf;
}

float SpecularReflection_sample_f(vec3 f, vec3 wi,
                                  const vec3 wo, const vec2 sample, const SpecularReflection* spec_ref)
{
    // Assume in tangent space
    vec3_assign(wi, -wo[0], -wo[1], wo[2]);
    vec3_scale(f, spec_ref->cr, 1.0f / fabs(cosTheta(wi)));
    return 1.0f;
}

float SpecularTransmission_sample_f(vec3 f, vec3 wi,
                                    const vec3 wo, const vec2 sample, const SpecularTransmission* spec_trans)
{
    vec3 normal = {0.0f, 0.0f, 1.0f};
    float kr = calcFresnelDielectric(normal, wo, spec_trans->ior_in, spec_trans->ior_out);
    float rand_float = (float)rand() / (float)RAND_MAX;
    if(rand_float <= kr)
    {
        vec3_assign(wi, -wo[0], -wo[1], wo[2]);
        vec3_scale(f, WHITE, kr / fabs(cosTheta(wi)));
        //vec3_scale(f, WHITE, 1.0f / fabs(cosTheta(wi)));
        return kr;
    }else
    {
        float eta = spec_trans->ior_out / spec_trans->ior_in;
        if(vec3_dot(normal, wo) < 0.0f)
        {
            vec3_negate(normal, normal);
            eta = 1.0f / eta;
        }
        if(!refract(wi, wo, normal, eta))
            return 0.0f;
        vec3_scale(f, WHITE, 1.0f - kr);
        //vec3_scale(f, WHITE, 1.0f);
        vec3_scale(f, f, (eta*eta) / fabs(cosTheta(wi)));
        return 1.0f - kr;
    }
}

void MicrofacetReflection_f(vec3 f, const vec3 wi, const vec3 wo, const MicrofacetReflection* mf)
{
    float cos_theta_o = absCosTheta(wo), cos_theta_i = absCosTheta(wi);
    vec3 wh;
    vec3_add(wh, wi, wo);
    if(cos_theta_i == 0.0f || cos_theta_o == 0.0f) return;
    if(wh[0] == 0.0f && wh[1] == 0.0f && wh[2] == 0.0f) return;
    vec3_normalize(wh, wh);
    vec3 fresnel_reflectance;
    if(!mf->is_metal)
    {
        // NOTE not sure why it's 1-kr and not just kr
        float kr = 1.0f - calcFresnelDielectric(wh, wi, mf->ior_in, mf->ior_out);
        vec3_assign(fresnel_reflectance, kr, kr, kr);
    }else
    {
        calcFresnelConductor(fresnel_reflectance, wh, wi, mf->etaT, mf->etaI, mf->k);
    }
    /*
      // Old
    float scale_factor = MicrofacetDistribution_D(wh, &(mf->distrib)) *
        MicrofacetDistribution_G(wo, wi, &(mf->distrib)) * (1.0f - fresnel_reflectance) /
        (4.0f * cos_theta_i * cos_theta_o);
        vec3_scale(f, mf->color, scale_factor);
    */
    float scale_factor = MicrofacetDistribution_D(wh, &(mf->distrib)) *
        MicrofacetDistribution_G(wo, wi, &(mf->distrib)) /
        (4.0f * cos_theta_i * cos_theta_o);    
    vec3_scale(f, mf->color, scale_factor);
    vec3_mult(f, f, fresnel_reflectance);
}

float MicrofacetReflection_sample_f(vec3 f, vec3 wi, const vec3 wo, const vec2 sample,
                                    const MicrofacetReflection* mr)
{
    if(wo[2] == 0.0f) return 0.0f;
    vec3 wh;
    MicrofacetDistribution_sample_wh(wh, wo, sample, &(mr->distrib));
    vec3 neg_wo;
    vec3_negate(neg_wo, wo);
    calcReflectRayDir(wi, wh, neg_wo);
    if(!sameHemisphere(wo, wi))
    {
        return 0.0f;
    }    
    MicrofacetReflection_f(f, wi, wo, mr);
    return MicrofacetDistribution_pdf(wo, wh, &(mr->distrib)) / (4.0f * vec3_dot(wo, wh));
}

float MicrofacetReflection_pdf(const vec3 wi, const vec3 wo, const MicrofacetReflection* mr)
{
    if(!sameHemisphere(wi, wo)) return 0.0f;
    vec3 wh;
    vec3_add(wh, wo, wi);
    vec3_normalize(wh, wh);
    return MicrofacetDistribution_pdf(wo, wh, &(mr->distrib)) / (4.0f * vec3_dot(wo, wh));    
}

// TODO name doesn't quite make sense
void MicrofacetFresnel_f(vec3 f, const vec3 wi, const vec3 wo, const MicrofacetFresnel * mt)
{
    if(sameHemisphere(wo, wi)) return;

    float cos_theta_o = cosTheta(wo);
    float cos_theta_i = cosTheta(wi);
    if(cos_theta_i == 0.0f || cos_theta_o == 0.0f)
    {
        vec3_copy(f, BLACK);
        return;
    }
    float eta = cos_theta_o > 0.0f ? (mt->ior_in / mt->ior_out) : (mt->ior_out / mt->ior_in);    
    vec3 wh, scaled_wi;
    vec3_scale(scaled_wi, wi, eta);
    vec3_add(wh, wo, scaled_wi);
    vec3_normalize(wh, wh);
    if(wh[2] < 0.0f) vec3_negate(wh, wh);

    float fresnel_reflectance = calcFresnelDielectric(wh, wo, mt->ior_in, mt->ior_out);
    float sqrt_denom = vec3_dot(wo, wh) + eta * vec3_dot(wi, wh);
    // TODO transport mode?
    float denom = cos_theta_i * cos_theta_o * sqrt_denom * sqrt_denom;
    const MicrofacetDistribution* distrib = &(mt->distrib);
    float factor = 1.0f / eta;
    float numerator = MicrofacetDistribution_D(wh, distrib) *
        MicrofacetDistribution_G(wo, wi, distrib) * eta * eta * fabs(vec3_dot(wi, wh)) * fabs(vec3_dot(wo, wh))
        * factor * factor;
    vec3 tmp;
    vec3_scale(tmp, mt->color, 1.0f - fresnel_reflectance);
    vec3_scale(f, tmp, fabs(numerator / denom));

}

float MicrofacetFresnel_sample_f(vec3 f, vec3 wi, const vec3 wo, const vec2 sample,
                                      const MicrofacetFresnel* mt)
{
    /*
    if (wo.z == 0) return 0.;
    Vector3f wh = distribution->Sample_wh(wo, u);
    Float eta = CosTheta(wo) > 0 ? (etaA / etaB) : (etaB / etaA);
    if (!Refract(wo, (Normal3f)wh, eta, wi)) return 0;
    *pdf = Pdf(wo, *wi);
    return f(wo, *wi);
     */
    /*
    if(wo[2] == 0) return 0.0f;
    vec3 wh;
    MicrofacetDistribution_sample_wh(wh, wo, sample, &(mt->distrib));
    float eta = cosTheta(wo) > 0.0f ? (mt->ior_out / mt->ior_in) :
        (mt->ior_in / mt->ior_out);
    if(!refract(wi, wo, wh, eta))
        return 0.0f;
    MicrofacetFresnel_f(f, wi, wo, mt);
    if(f[0] < 0.0f || f[1] < 0.0f || f[2] < 0.0f)
    {
            MicrofacetFresnel_f(f, wi, wo, mt);
    }
    return MicrofacetFresnel_pdf(wi, wo, mt);
    */
    if(wo[2] == 0.0f) return 0.0f;
    vec3 wh;
    MicrofacetDistribution_sample_wh(wh, wo, sample, &(mt->distrib));
    printVec3WithText("wh", wh);
    float kr = calcFresnelDielectric(wh, wo, mt->ior_in, mt->ior_out);
    float rand_float = (float)rand() / (float)RAND_MAX;
    if(rand_float <= kr)
    {
        vec3 neg_wo;
        vec3_negate(neg_wo, wo);
        calcReflectRayDir(wi, wh, neg_wo);
        if(!sameHemisphere(wo, wi))
            return 0.0f;
        MicrofacetReflection mr;
        vec3_copy(mr.color, mt->color);
        mr.ior_in = mt->ior_in;
        mr.ior_out = mt->ior_out;
        mr.distrib = mt->distrib;
        MicrofacetReflection_f(f, wi, wo, &mr);
        //vec3_scale(f, f, 1.0f / (1.0f - kr));
        return MicrofacetDistribution_pdf(wo, wh, &(mr.distrib)) / (4.0f * vec3_dot(wo, wh));
    }else
    {
        float eta = cosTheta(wo) > 0.0f ? (mt->ior_out / mt->ior_in) :
            (mt->ior_in / mt->ior_out);
        if(!refract(wi, wo, wh, eta))
            return 0.0f;
        MicrofacetFresnel_f(f, wi, wo, mt);
        //vec3_scale(f, f, 1.0f / (1.0f - kr));
        return MicrofacetFresnel_pdf(wi, wo, mt);
    }
}

float MicrofacetFresnel_pdf(const vec3 wi, const vec3 wo, const MicrofacetFresnel* mt)
{
    /*
    if (SameHemisphere(wo, wi)) return 0;
    // Compute $\wh$ from $\wo$ and $\wi$ for microfacet transmission
    Float eta = CosTheta(wo) > 0 ? (etaB / etaA) : (etaA / etaB);
    Vector3f wh = Normalize(wo + wi * eta);

    // Compute change of variables _dwh\_dwi_ for microfacet transmission
    Float sqrtDenom = Dot(wo, wh) + eta * Dot(wi, wh);
    Float dwh_dwi =
        std::abs((eta * eta * Dot(wi, wh)) / (sqrtDenom * sqrtDenom));
    return distribution->Pdf(wo, wh) * dwh_dwi;
     */
    if(sameHemisphere(wo, wi)) return 0.0f;
    float eta = cosTheta(wo) > 0.0f ? (mt->ior_in / mt->ior_out) : (mt->ior_out / mt->ior_in);    
    vec3 wh, scaled_wi;
    vec3_scale(scaled_wi, wi, eta);
    vec3_add(wh, wo, scaled_wi);
    vec3_normalize(wh, wh);
    float sqrt_denom = vec3_dot(wo, wh) + eta * vec3_dot(wi, wh);
    float dwh_dwi = fabs((eta * eta * vec3_dot(wi, wh)) / (sqrt_denom * sqrt_denom));
    return MicrofacetDistribution_pdf(wo, wh, &(mt->distrib)) * dwh_dwi;
}

void BxDF_f(vec3 f, const vec3 wi, const vec3 wo, const void* bxdf, const BxDFType type)
{
    switch(type)
    {
    case LAMBERTIAN:
    {
        return Lambertian_f(f, wi, wo, (Lambertian*)bxdf);
    } break;
    case ORENNAYAR:
    {
        return OrenNayar_f(f, wi, wo, (OrenNayar*)bxdf);
    } break;
    case SPECULAR_REFLECTION:
    case SPECULAR_TRANSMISSION:
        break;
    case MICROFACET_REFLECTION:
    {
        return MicrofacetReflection_f(f, wi, wo, (MicrofacetReflection*)bxdf);
    } break;
    case MICROFACET_FRESNEL:
    {
        return MicrofacetFresnel_f(f, wi, wo, (MicrofacetFresnel*)bxdf);
    } break;
    }    
}

float BxDF_sample_f(vec3 f, vec3 wi, const vec3 wo, const vec2 sample, const void* bxdf, const BxDFType type)
{
    switch(type)
    {
    case LAMBERTIAN:
    {
        return Lambertian_sample_f(f, wi, wo, sample, (Lambertian*)bxdf);
    } break;
    case ORENNAYAR:
    {
        return OrenNayar_sample_f(f, wi, wo, sample, (OrenNayar*)bxdf);
    } break;
    case SPECULAR_REFLECTION:
    {
        return SpecularReflection_sample_f(f, wi, wo, sample, (SpecularReflection*)bxdf);
    } break;
    case SPECULAR_TRANSMISSION:
    {
        return SpecularTransmission_sample_f(f, wi, wo, sample, (SpecularTransmission*)bxdf);
    } break;
    case MICROFACET_REFLECTION:
    {
        return MicrofacetReflection_sample_f(f, wi, wo, sample, (MicrofacetReflection*)bxdf);
    } break;
    case MICROFACET_FRESNEL:
    {
        return MicrofacetFresnel_sample_f(f, wi, wo, sample, (MicrofacetFresnel*)bxdf);
    } break;
    }
    return 0.0f;
}

float BxDF_pdf(const vec3 wi, const vec3 wo, const void* bxdf, const BxDFType type)
{
    switch(type)
    {
    case LAMBERTIAN:
    {
        return Lambertian_pdf(wi, wo);
    } break;
    case ORENNAYAR:
    {
        return OrenNayar_pdf(wi, wo);
    }
    case SPECULAR_REFLECTION:
    {
        return 0.0f;
    } break;
    case SPECULAR_TRANSMISSION:
    {
        return 0.0f;
    } break;
    case MICROFACET_REFLECTION:
    {
        MicrofacetReflection* mr = (MicrofacetReflection*)bxdf;
        return MicrofacetReflection_pdf(wi, wo, mr);
    }
    case MICROFACET_FRESNEL:
    {
        MicrofacetFresnel* mt = (MicrofacetFresnel*)bxdf;
        return MicrofacetFresnel_pdf(wi, wo, mt);
    }
    }
    return 0.0f;
}

void BSDF_f(vec3 f, const vec3 wi, const vec3 wo, const BSDF* bsdf, const BxDFFlags excluded)
{
    vec3 wi_local, wo_local;
    orthoNormalTransform(wi_local, bsdf->tangent, bsdf->binormal, bsdf->normal, wi);
    orthoNormalTransform(wo_local, bsdf->tangent, bsdf->binormal, bsdf->normal, wo);
    vec3_copy(f, BLACK);
    for(int i = 0; i < bsdf->num_bxdf; i++)
    {
        BxDFType type = bsdf->types[i];
        if(!(getBxDFFlagsFromType(type) & excluded))
        {
            vec3 cur_f = {0.0f, 0.0f, 0.0f};        
            BxDF_f(cur_f, wi, wo, bsdf->bxdfs[i], bsdf->types[i]);
            vec3_add(f, f, cur_f);
        }
    }
}

float BSDF_pdf(const vec3 wi, const vec3 wo, const BSDF* bsdf)
{
    vec3 wi_local, wo_local;
    transposeTransform(wo_local, bsdf->tangent, bsdf->binormal, bsdf->normal, wo);
    transposeTransform(wo_local, bsdf->tangent, bsdf->binormal, bsdf->normal, wo);
    float pdf = 0;
    for(int i = 0; i < bsdf->num_bxdf; i++)
    {
        pdf += BxDF_pdf(wi, wo, bsdf->bxdfs[i], bsdf->types[i]);
    }
    return pdf;
}

float BSDF_sample_f(vec3 f, vec3 wi, BxDFFlags* sampled_flags,
                    const vec3 wo, const vec2 sample, const BSDF* bsdf)
{
    if(bsdf->num_bxdf == 0)
    {
        vec3_copy(f, BLACK);
        *sampled_flags = BSDF_NONE;
        printf("0 bxdf\n");
        return 0.0f;
    }
    // Choose BxDF
    int bxdf_index = (int)(sample[0] * bsdf->num_bxdf);
    BxDFType bxdf_type = bsdf->types[bxdf_index];
    *sampled_flags = getBxDFFlagsFromType(bxdf_type);
    vec2 remapped_sample = {sample[0] * bsdf->num_bxdf - bxdf_index, sample[1]};

    // Transform wo to tangent space
    vec3 wi_local, wo_local;
    transposeTransform(wo_local, bsdf->tangent, bsdf->binormal, bsdf->normal, wo);
    if(wo_local[2] == 0.0f)
    {
        printf("wo parallel\n");
        return 0.0f;
    }
    
    // Sample BxDF
    vec3_copy(f, BLACK);
    float pdf = BxDF_sample_f(f, wi_local, wo_local, remapped_sample,
                              bsdf->bxdfs[bxdf_index], bsdf->types[bxdf_index]);
    if(pdf == 0.0f)
    {
        vec3_copy(f, BLACK);
        return pdf;
    }
    orthoNormalTransform(wi, bsdf->tangent, bsdf->binormal, bsdf->normal, wi_local);

    // Add pdfs from other BxDFs
    //if(!(bsdf->types[bxdf_index] == SPECULAR_REFLECTION || bsdf->types[bxdf_index] == SPECULAR_TRANSMISSION))
    if(!(*sampled_flags & BSDF_SPECULAR))
    {
        for(int i = 0; i < bsdf->num_bxdf; i++)
        {
            if(i != bxdf_index)
            {
                pdf += BxDF_pdf(wi_local, wo_local, bsdf->bxdfs[i], bsdf->types[i]);
            }
        }
    }

    // Compute f for rest of the BxDFs
    for(int i = 0; i < bsdf->num_bxdf; i++)
    {
        if(i != bxdf_index)
        {
            vec3 cur_f = {0.0f, 0.0f, 0.0f};
            BxDF_f(cur_f, wi_local, wo_local, bsdf->bxdfs[i], bsdf->types[i]);
            vec3_add(f, f, cur_f);
        }
    }
    return pdf;
}

void BSDF_addLambertian(BSDF* bsdf, const vec3 cd)
{
    Lambertian* l = (Lambertian*)allocateBxDF();
    vec3_copy(l->cd, cd);
    bsdf->bxdfs[bsdf->num_bxdf] = l;
    bsdf->types[bsdf->num_bxdf] = LAMBERTIAN;
    bsdf->num_bxdf++;
}


void BSDF_addOrenNayar(BSDF* bsdf, const vec3 r, const float sigma)
{
    float sigma_rad = degToRad(sigma);
    float sigma2 = sigma_rad * sigma_rad;
    OrenNayar* on = (OrenNayar*)allocateBxDF();
    on->a = 1.0f - (sigma2 / (2.0f * (sigma2 + 0.33f)));
    on->b = 0.45f * sigma2 / (sigma2 + 0.09f);
    vec3_copy(on->r, r);

    bsdf->bxdfs[bsdf->num_bxdf] = on;
    bsdf->types[bsdf->num_bxdf] = ORENNAYAR;
    bsdf->num_bxdf++;
}

void BSDF_addSpecularReflection(BSDF* bsdf, const vec3 cr)
{

    SpecularReflection* spec_ref = (SpecularReflection*)allocateBxDF();
    vec3_copy(spec_ref->cr, cr);
    bsdf->bxdfs[bsdf->num_bxdf] = spec_ref;
    bsdf->types[bsdf->num_bxdf] = SPECULAR_REFLECTION;
    bsdf->num_bxdf++;
}

void BSDF_addSpecularTransmission(BSDF* bsdf, const float ior_in, const float ior_out, const vec3 cf_in,
                                  const vec3 cf_out)
{

    SpecularTransmission* spec_trans = (SpecularTransmission*)allocateBxDF();
    spec_trans->ior_in = ior_in;
    spec_trans->ior_out = ior_out;
    vec3_copy(spec_trans->cf_in, cf_in);
    vec3_copy(spec_trans->cf_out, cf_out);
    bsdf->bxdfs[bsdf->num_bxdf] = spec_trans;
    bsdf->types[bsdf->num_bxdf] = SPECULAR_TRANSMISSION;
    bsdf->num_bxdf++;
}

void BSDF_addMicrofacetReflection(BSDF* bsdf, const vec3 color, const float ior_in, const float ior_out,
                                  const float alphax, const float alphay, const FacetDistribType type)
{
    MicrofacetReflection* mr = (MicrofacetReflection*)allocateBxDF();
    vec3_copy(mr->color, color);
    //mr->distrib.alphax = alphax;
    //mr->distrib.alphay = alphay;
    mr->distrib.alphax = BeckmannRoughnessToAlpha(alphax);
    mr->distrib.alphay = BeckmannRoughnessToAlpha(alphay);
    mr->ior_in = ior_in;
    mr->ior_out = ior_out;
    mr->distrib.type = type;
    mr->is_metal = false;
    bsdf->bxdfs[bsdf->num_bxdf] = mr;
    bsdf->types[bsdf->num_bxdf] = MICROFACET_REFLECTION;
    bsdf->num_bxdf++;
}

void BSDF_addMicrofacetReflectionMetal(BSDF* bsdf, const vec3 color, const vec3 etaT, const vec3 etaI,
                                       const vec3 k, const float alphax, const float alphay,
                                       const FacetDistribType type)
{
    MicrofacetReflection* mr = (MicrofacetReflection*)allocateBxDF();
    vec3_copy(mr->color, color);
    mr->distrib.alphax = BeckmannRoughnessToAlpha(alphax);
    mr->distrib.alphay = BeckmannRoughnessToAlpha(alphay);
    vec3_copy(mr->etaT, etaT);
    vec3_copy(mr->etaI, etaI);
    vec3_copy(mr->k, k);
    mr->distrib.type = type;
    mr->is_metal = true;
    bsdf->bxdfs[bsdf->num_bxdf] = mr;
    bsdf->types[bsdf->num_bxdf] = MICROFACET_REFLECTION;
    bsdf->num_bxdf++;
}

void BSDF_addMicrofacetFresnel(BSDF* bsdf, const vec3 color, const float ior_in, const float ior_out,
                                    const float alphax, const float alphay, const FacetDistribType type)
{
    MicrofacetFresnel* mt = (MicrofacetFresnel*)allocateBxDF();
    vec3_copy(mt->color, color);
    mt->distrib.alphax = alphax;
    mt->distrib.alphay = alphay;
    mt->ior_in = ior_in;
    mt->ior_out = ior_out;
    mt->distrib.type = type;
    bsdf->bxdfs[bsdf->num_bxdf] = mt;
    bsdf->types[bsdf->num_bxdf] = MICROFACET_FRESNEL;
    bsdf->num_bxdf++;    
}

void BSDF_freeBxDFs(BSDF* bsdf)
{
    for(int i = 0; i < bsdf->num_bxdf; i++)
    {
        freeBxDF(&(bsdf->bxdfs[i]));
    }
    bsdf->num_bxdf = 0;
}

static MemPool bsdf_mem_pool;

bool initBSDFMem(const int num_threads, const int num_depth)
{
    float bxdf_size = max(sizeof(SpecularTransmission), sizeof(MicrofacetReflection));
    MemPool_init(&bsdf_mem_pool, bxdf_size, num_threads * num_depth * MAX_BXDF);    
}

void freeBSDFMem()
{
    MemPool_destroy(&bsdf_mem_pool);
}

void* allocateBxDF()
{
    return MemPool_requestElement(&bsdf_mem_pool);
}

void freeBxDF(void** bxdf)
{
    MemPool_releaseElement(&bsdf_mem_pool, bxdf);
}
