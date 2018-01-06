#include "materials.h"
#include "util/constants.h"
#include <string.h>

MatType getMatTypeFromString(const char* str)
{
    if(strcmp(str, "MATTE") == 0)
    {
        return MATTE;
    }else if(strcmp(str, "PHONG") == 0)
    {
        return PHONG;
    }else if(strcmp(str, "REFLECTIVE") == 0)
    {
        return REFLECTIVE;
    }else if(strcmp(str, "TRANSPARENT") == 0)
    {
        return TRANSPARENT;
    }else if(strcmp(str, "EMISSIVE") == 0)
    {
        return EMISSIVE;
    }else if(strcmp(str, "PARTICIPATING") == 0)
    {
        return PARTICIPATING;
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

void Reflective_getColor(vec3 color, const Reflective* ref)
{
    vec3_copy(color, ref->color);
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
    case REFLECTIVE:
    {
        vec3 color;
        Reflective* ref = (Reflective*)mat->data;        
        Reflective_getColor(color, ref);
        //printVec3WithText("ref color", color);
        BSDF_addSpecularReflection(bsdf, color);
    } break;
    case TRANSPARENT:
    {
        Transparent* trans = (Transparent*)mat->data;
        BSDF_addSpecularTransmission(bsdf, trans->ior_in, trans->ior_out, trans->cf_in, trans->cf_out);
    } break;
    }
}

bool Material_hasNormalMap(const Material* mat)
{
    switch(mat->mat_type)
    {
    case MATTE:
        return true;
    case REFLECTIVE:
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
