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
    printf("ks %f, ke %f\n", mat->ks, mat->ke);
    printVec3WithText("ca", mat->cd);
    printVec3WithText("cd", mat->cd);
    printVec3WithText("cs", mat->cs);
    printVec3WithText("ce", mat->cs);
    printf("ior_in %f, ior_out %f\n", mat->ior_in, mat->ior_out);
    printVec3WithText("cf_in", mat->cf_in);
    printVec3WithText("cf_out", mat->cf_out);
}

void initMaterial(Material *mat)
{
    mat->shadow = false;
    mat->mat_type = INVALID_MAT_TYPE;
    mat->tex_flags = 0;
    mat->ka = mat->kd = mat->ks = mat->ke = mat->kr = mat->exp = 0.0f;
    mat->ior_in = mat->ior_out = 0.0f;
    vec3_copy(mat->ca, BLACK);
    vec3_copy(mat->cd, BLACK);
    vec3_copy(mat->cs, BLACK);
    vec3_copy(mat->cr, BLACK);
    vec3_copy(mat->ce, BLACK);
    vec3_copy(mat->cf_in, BLACK);
    vec3_copy(mat->cf_out, BLACK);
    mat->h_samples = NULL;
    for(int i = 0; i < 3; i++)
    {
        mat->tex_array[i] = NULL;
    }
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

void getMaterialDiffuseTexColor(vec3 texel, const Material *mat, const vec2 uv)
{
    getTexColor(texel, mat->tex_array[DIFFUSE_MAP_INDEX], uv);
}

void getMaterialNormalTexColor(vec3 texel, const Material *mat, const vec2 uv)
{
    getTexColor(texel, mat->tex_array[NORMAL_MAP_INDEX], uv);
}

void setMaterialDiffuseTexPtr(Material *mat, Texture *tex)
{
    mat->tex_array[DIFFUSE_MAP_INDEX] = tex;
    mat->tex_flags |= DIFFUSE;
}

void setMaterialNormalTexPtr(Material *mat, Texture *tex)
{
    mat->tex_array[NORMAL_MAP_INDEX] = tex;
    mat->tex_flags |= NORMAL;
}
