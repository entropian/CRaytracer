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

void printMaterial(const Material* mat)
{
    printf("shadow %s\n", mat->shadow ? "true" : "false");
    printf("ka %f, kd %f\n", mat->ka, mat->kd);
    printVec3WithText("cd", mat->cd);
    printVec3WithText("cs", mat->cs);
    printf("ior_in %f, ior_out %f\n", mat->ior_in, mat->ior_out);
    printVec3WithText("cf_in", mat->cf_in);
    printVec3WithText("cf_out", mat->cf_out);
}

void initDefaultMatteMat(Material* mat, const vec3 color)
{
    vec3_copy(mat->cd, color);
    vec3_copy(mat->ca, color);
    mat->kd = 0.6f;
    mat->ka = 1.0f;
    mat->mat_type  = MATTE;
    mat->shadow = true;
    mat->h_samples = genHemisphereSamples(MULTIJITTERED, 1.0f);
    mat->tex_flags = NO_TEXTURE;
}

void initDefaultPhongMat(Material* mat, const vec3 color)
{
    vec3_copy(mat->cd, color);
    vec3_copy(mat->cs, color);
    vec3_copy(mat->ca, color);    
    mat->kd = 0.6f;
    mat->ka = 1.0f;
    mat->ks = 0.4f;
    mat->exp = 10.0f;            
    mat->mat_type = PHONG;
    mat->shadow = true;
    mat->h_samples = genHemisphereSamples(MULTIJITTERED, DEFAULT_GLOSSINESS);
    mat->tex_flags = NO_TEXTURE;    
}
