#include "materials.h"
#include "util/constants.h"
#include <string.h>

// Metal properties
const vec3 COPPER_N = {0.19999069, 0.92208463, 1.09987593};
const vec3 GOLD_ETA = {0.14282006, 0.37414363, 1.43944442};
const vec3 GOLD_K = {3.90463543, 2.44763327, 2.13765264};
const vec3 SILVER_ETA = {0.154935181, 0.116475478, 0.138087392};
const vec3 SILVER_K = {4.81810093, 3.11561656, 2.1424017};
const vec3 BERYLLIUM_ETA = {4.17617416, 3.1783011, 2.77819276};
const vec3 BERYLLIUM_K = {3.82729554, 3.00373626, 2.86292768};
const vec3 CHROMIUM_ETA = {4.36040831, 2.9105196, 1.65118635};
const vec3 CHROMIUM_K = {5.19538164, 4.22238398, 3.74699736};
const vec3 CESIUM_ETA = {2.14034843, 1.69870293, 1.65889668};
const vec3 CESIUM_K = {0.0f, 0.0f, 0.0f};
const vec3 COPPER_ETA = {0.19999069, 0.92208463, 1.09987593};
const vec3 COPPER_K = {3.90463543, 2.44763327, 2.13765264};
const vec3 MERCURY_ETA = {2.39383841, 1.43696785, 0.907622635};
const vec3 MERCURY_K = {6.31419611, 4.36266136, 3.41453838};

MatType getMatTypeFromString(const char* str)
{
    if(strcmp(str, "MATTE") == 0)
    {
        return MATTE;
    }else if(strcmp(str, "MIRROR") == 0)
    {
        return MIRROR;
    }else if(strcmp(str, "TRANSPARENT") == 0)
    {
        return TRANSPARENT;
    }else if(strcmp(str, "EMISSIVE") == 0)
    {
        return EMISSIVE;
    }else if(strcmp(str, "PLASTIC") == 0)
    {
        return PLASTIC;
    }else if(strcmp(str, "GLASS") == 0)
    {
        return GLASS;
    }else if(strcmp(str, "METAL") == 0)
    {
        return METAL;
    }else
    {
        return INVALID_MAT_TYPE;
    }
}

void Matte_getDiffuseColor(vec3 color, const Matte* matte, const vec2 uv)
{
    if(matte->diffuse)
    {
        getTexColor(color, matte->diffuse, uv);
    }else
    {
        vec3_copy(color, matte->color);
    }
}
void Matte_getNormalMapColor(vec3 color, const Matte* matte, const vec2 uv)
{
    if(matte->normal)
    {
        getTexColor(color, matte->normal, uv);
    }
}

void Mirror_getColor(vec3 color, const Mirror* mirror)
{
    vec3_copy(color, mirror->color);
}

void procMetalType(Metal* metal, const char* type)
{
    if(strcmp(type, "GOLD") == 0)
    {
        vec3_copy(metal->eta, GOLD_ETA);
        vec3_copy(metal->k, GOLD_K);
    }else if(strcmp(type, "SILVER") == 0)
    {
        vec3_copy(metal->eta, SILVER_ETA);
        vec3_copy(metal->k, SILVER_K);
    }else if(strcmp(type, "BERYLLIUM") == 0)
    {
        vec3_copy(metal->eta, BERYLLIUM_ETA);
        vec3_copy(metal->k, BERYLLIUM_K);
    }else if(strcmp(type, "CHROMIUM") == 0)
    {
        vec3_copy(metal->eta, CHROMIUM_ETA);
        vec3_copy(metal->k, CHROMIUM_K);
    }else if(strcmp(type, "CESIUM") == 0)
    {
        vec3_copy(metal->eta, CESIUM_ETA);
        vec3_copy(metal->k, CESIUM_K);
    }else if(strcmp(type, "COPPER") == 0)
    {
        vec3_copy(metal->eta, COPPER_ETA);
        vec3_copy(metal->k, COPPER_K);

    }else if(strcmp(type, "MERCURY") == 0)
    {
        vec3_copy(metal->eta, MERCURY_ETA);
        vec3_copy(metal->k, MERCURY_K);                
    }else
    {
        fprintf(stderr, "Invalid metal type.\n");
    }
}

void computeScatteringFunc(BSDF* bsdf, const vec2 uv, const Material* mat)
{
    bsdf->num_bxdf = 0;
    switch(mat->mat_type)
    {
    case MATTE:
    {
        vec3 color;
        Matte* matte = (Matte*)mat->data;
        Matte_getDiffuseColor(color, matte, uv);
        if(matte->sigma == 0.0f)
        {
            BSDF_addLambertian(bsdf, color);
        }else
        {
            BSDF_addOrenNayar(bsdf, color, matte->sigma);
        }        
    } break;
    case MIRROR:
    {
        vec3 color;
        Mirror* ref = (Mirror*)mat->data;        
        Mirror_getColor(color, ref);
        //printVec3WithText("ref color", color);
        BSDF_addSpecularReflection(bsdf, color);
    } break;
    case TRANSPARENT:
    {
        Transparent* trans = (Transparent*)mat->data;
        BSDF_addSpecularTransmission(bsdf, trans->ior_in, trans->ior_out, trans->cf_in, trans->cf_out);
    } break;
    case PLASTIC:
    {
        Plastic* plastic = (Plastic*)mat->data;
        if(!vec3_equal(plastic->kd, BLACK))
            BSDF_addLambertian(bsdf, plastic->kd);
        if(!vec3_equal(plastic->ks, BLACK))
        {
            if(plastic->roughness == 0.0f)
            {
                BSDF_addSpecularReflection(bsdf, plastic->ks);
            }else
            {
                //printf("roughness %f\n", plastic->roughness);
                BSDF_addMicrofacetReflection(bsdf, plastic->ks, 1.5f, 1.0f,
                                             plastic->roughness, plastic->roughness, BECKMANN);
            }
        }
    } break;
    case GLASS:
    {
        Glass* glass = (Glass*)mat->data;
        vec3 color = {1.0f, 1.0f, 1.0f};
        BSDF_addMicrofacetFresnel(bsdf, color, glass->ior_in, glass->ior_out, glass->uroughness,
                                  glass->vroughness, BECKMANN);
    } break;
    case METAL:
    {
        Metal* metal = (Metal*)mat->data;
        vec3 color = {1.0f, 1.0f, 1.0f};
        vec3 etaI = {1.0f, 1.0f, 1.0f};
        BSDF_addMicrofacetReflectionMetal(bsdf, color, metal->eta, etaI, metal->k, metal->uroughness,
                                          metal->vroughness, BECKMANN);
    } break;
    }
}

bool Material_hasNormalMap(const Material* mat)
{
    switch(mat->mat_type)
    {
    case MATTE:
        return true;
    case MIRROR:
        return false;
    case TRANSPARENT:
        return false;
    }
    return false;
}

void Material_getNormalMapValue(vec3 value, const Material* mat, const vec2 uv)
{
    switch(mat->mat_type)
    {
    case MATTE:
    {
        Matte* matte = (Matte*)mat->data;
        Matte_getNormalMapColor(value, matte, uv);
    } break;
    }
}
