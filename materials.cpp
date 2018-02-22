#include "materials.h"
#include "util/constants.h"
#include <string.h>

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
        // TODO: add texture for matte
        if(!vec3_equal(plastic->kd, BLACK))
            BSDF_addLambertian(bsdf, plastic->kd);
        if(!vec3_equal(plastic->ks, BLACK))
        {
            if(plastic->roughness == 0.0f)
            {
                BSDF_addSpecularReflection(bsdf, plastic->ks);
            }else
            {
                BSDF_addMicrofacetReflection(bsdf, plastic->ks, 1.5f, 1.0f,
                                             plastic->roughness, plastic->roughness, BECKMANN);
            }
        }
    } break;
    case GLASS:
    {
        /*
          typedef struct Glass_s
          {
          vec3 kr, kt;
          float ior_in, ior_out;
          float uroughness, vroughness;
          }Glass;
         */
        Glass* glass = (Glass*)mat->data;
        bool is_perfect_specular = glass->uroughness == 0.0f && glass->vroughness == 0.0f;
        float remapped_uroughness = BeckmannRoughnessToAlpha(glass->uroughness);
        float remapped_vroughness = BeckmannRoughnessToAlpha(glass->vroughness);
        if(!vec3_equal(glass->kr, BLACK))
        {
            if(is_perfect_specular)
            {
                BSDF_addSpecularReflection(bsdf, glass->kr);
            }else
            {
                BSDF_addMicrofacetReflection(bsdf, glass->kr, glass->ior_in, glass->ior_out,
                                             remapped_uroughness, remapped_vroughness, BECKMANN);
            }
        }
        if(!vec3_equal(glass->kt, BLACK))
        {
            if(is_perfect_specular)
            {
                // TODO resolve transmission color issue
                BSDF_addSpecularTransmission(bsdf, glass->ior_in, glass->ior_out, glass->kr, glass->kr);
            }else
            {
                BSDF_addMicrofacetTransmission(bsdf, glass->kt, glass->ior_in, glass->ior_out,
                                               remapped_uroughness, remapped_vroughness, BECKMANN);
            }
        }
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
