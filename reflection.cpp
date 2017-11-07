#include "reflection.h"
#include "sampling.h"
#include "util/ray.h"
#include "util/math.h"

bool sameHemisphere(const vec3 a, const vec3 b)
{
    return a[2] * b[2] > 0.0f;
}

float calcFresnelReflectance(const vec3 normal, const vec3 wo, const float ior_in, const float ior_out)
{
    vec3 n;
    vec3_copy(n, normal);
    // NOTE: not sure about the next line
    float cos_theta_i = vec3_dot(n, wo);
    float eta = ior_in / ior_out;    
    if(cos_theta_i < 0.0f)
    {
        cos_theta_i = -cos_theta_i;
        vec3_negate(n, n);
        eta = 1.0f / eta;
    }
    float tmp = 1.0f - (1.0f - cos_theta_i * cos_theta_i) / (eta * eta);
    if(tmp < 0.0f)
    {
        return 1.0f;
    }
    float cos_theta_t = (float)sqrt(tmp);
    float r_parallel = (eta * cos_theta_i - cos_theta_t) / (eta * cos_theta_i + cos_theta_t);
    float r_perpendicular = (cos_theta_i - eta * cos_theta_t) / (cos_theta_i + eta * cos_theta_t);
    float fresnel_reflectance = 0.5f * (r_parallel * r_parallel + r_perpendicular * r_perpendicular);
    return fresnel_reflectance;
}

void Lambertian_f(vec3 f, const vec3 wi, const vec3 wo, const Lambertian* l)
{
    vec3_copy(f, l->cd);
    vec3_scale(f, f, INV_PI);
}

float Lambertian_pdf(const vec3 wi, const vec3 wo)
{
    vec3 normal = {0.0f, 0.0f, 1.0f};
    //return sameHemisphere(wi, wo) ? fabs(vec3_dot(wi, wo)) * INV_PI: 0.0f;
    if(!sameHemisphere(wi, wo))
    {
        printf("not same hemisphere\n");
    }
    float dot_product = fabs(vec3_dot(wi, normal)) * INV_PI;
    if(dot_product == 0.0f)
    {
        printf("dot product 0\n");
    }
    
    return sameHemisphere(wi, wo) ? fabs(vec3_dot(wi, normal)) * INV_PI: 0.0f;
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
    mapSampleToHemisphere(wi, sample);
    return Lambertian_pdf(wi, tmp_wo);
}

float SpecularReflection_sample_f(vec3 f, vec3 wi,
                                  const vec3 wo, const vec2 sample, const SpecularReflection* spec_ref)
{
    // Assume in tangent space
    vec3_assign(wi, -wo[0], -wo[1], wo[2]);
    vec3_copy(f, spec_ref->cr);
    return 1.0f;
}

float SpecularTransmission_sample_f(vec3 f, vec3 wi,
                                    const vec3 wo, const vec2 sample, const SpecularTransmission* spec_trans)
{
    vec3 normal = {0.0f, 0.0f, 1.0f};
    float kr = calcFresnelReflectance(normal, wo, spec_trans->ior_in, spec_trans->ior_out);
    float rand_float = (float)rand() / (float)RAND_MAX;
    if(rand_float <= kr)
    {
        vec3_assign(wi, -wo[0], -wo[1], wo[2]);
        vec3_copy(f, WHITE);
    }else
    {
        float eta = calcTransmitDir(wi, normal, wo, spec_trans->ior_in, spec_trans->ior_out);
        float ndotwi = fabs(vec3_dot(normal, wi));
        //float kt = 1.0f - kr;
        //vec3_scale(f, WHITE, kt / (eta*eta) / ndotwi);
        //vec3_scale(f, WHITE, kt / (eta*eta));
        vec3_scale(f, WHITE, (eta*eta));
    }
    return 1.0f;
}

void BxDF_f(vec3 f, const vec3 wi, const vec3 wo, const void* bxdf, const BxDFType type)
{
    switch(type)
    {
    case LAMBERTIAN:
    {
        return Lambertian_f(f, wi, wo, (Lambertian*)bxdf);
    } break;
    case SPECULAR_REFLECTION:
    case SPECULAR_TRANSMISSION:
        break;
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
    case SPECULAR_REFLECTION:
    {
        return SpecularReflection_sample_f(f, wi, wo, sample, (SpecularReflection*)bxdf);
    } break;
    case SPECULAR_TRANSMISSION:
    {
        return SpecularTransmission_sample_f(f, wi, wo, sample, (SpecularTransmission*)bxdf);
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
    case SPECULAR_REFLECTION:
    {
        return 0.0f;
    } break;
    case SPECULAR_TRANSMISSION:
    {
        return 0.0f;
    } break;
    }
    return 0.0f;
}

void BSDF_f(vec3 f, const vec3 wi, const vec3 wo, const BSDF* bsdf)
{
    vec3 wi_local, wo_local;
    orthoNormalTransform(wi_local, bsdf->tangent, bsdf->binormal, bsdf->normal, wi);
    orthoNormalTransform(wo_local, bsdf->tangent, bsdf->binormal, bsdf->normal, wo);
    vec3_copy(f, BLACK);
    for(int i = 0; i < bsdf->num_bxdf; i++)
    {
        vec3 cur_f = {0.0f, 0.0f, 0.0f};
        BxDF_f(f, wi, wo, bsdf->bxdfs[i], bsdf->types[i]);
        vec3_add(f, f, cur_f);
    }
}

float BSDF_pdf(const vec3 wi, const vec3 wo, const BSDF* bsdf)
{
    vec3 wi_local, wo_local;
    orthoNormalTransform(wi_local, bsdf->tangent, bsdf->binormal, bsdf->normal, wi);
    orthoNormalTransform(wo_local, bsdf->tangent, bsdf->binormal, bsdf->normal, wo);
    float pdf = 0;
    for(int i = 0; i < bsdf->num_bxdf; i++)
    {
        pdf += BxDF_pdf(wi, wo, bsdf->bxdfs[i], bsdf->types[i]);
    }
    return pdf;
}

float BSDF_sample_f(vec3 f, vec3 wi,
                    const vec3 wo, const vec2 sample, const BSDF* bsdf)
{
    if(bsdf->num_bxdf == 0)
    {
        vec3_copy(f, BLACK);
        printf("0 bxdf\n");
        return 0.0f;
    }
    // Choose BxDF
    int bxdf_index = (int)(sample[0] * bsdf->num_bxdf);    
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
        //printf("pdf 0\n");
        return pdf;
    }
    //transposeTransform(wi, bsdf->tangent, bsdf->binormal, bsdf->normal, wi_local);
    orthoNormalTransform(wi, bsdf->tangent, bsdf->binormal, bsdf->normal, wi_local);

    // Add pdfs from other BxDFs
    if(bsdf->types[bxdf_index] == SPECULAR_REFLECTION ||
       bsdf->types[bxdf_index] == SPECULAR_TRANSMISSION)
        return pdf;
    for(int i = 0; i < bsdf->num_bxdf; i++)
    {
        if(i != bxdf_index)
        {
            pdf += BxDF_pdf(wi_local, wo_local, bsdf->bxdfs[i], bsdf->types[i]);
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

void BSDF_freeBxDFs(BSDF* bsdf)
{
    for(int i = 0; i < bsdf->num_bxdf; i++)
    {
        freeBxDF(&(bsdf->bxdfs[i]));
        bsdf->num_bxdf = 0;
    }
}

static MemPool bsdf_mem_pool;

bool initBSDFMem(const int num_threads, const int num_depth)
{
    MemPool_init(&bsdf_mem_pool, sizeof(SpecularTransmission), num_threads * num_depth * MAX_BXDF);    
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
