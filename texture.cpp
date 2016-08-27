#include "texture.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

bool loadTexture(Texture* tex, const char* file_name)
{
    tex->data = stbi_load(file_name, &(tex->width), &(tex->height), &(tex->comp), STBI_rgb);
    if(!tex->data)
    {
        fprintf(stderr, "Failed to load %s.\n", file_name);
        return false;
    }
    return true;
}

void freeTexture(Texture* tex)
{
    if(tex->data)
    {
        free(tex->data);
    }
    tex->width = tex->height = tex->comp = 0;
}

bool getTexColor(vec3 color_out, const Texture* texture, const vec2 texcoord)
{    
    float u_float = texcoord[0];
    if(u_float < 0.0f)
    {
        u_float = -u_float;
    }
    if(u_float > 1.0f)
    {
        u_float = u_float - (int)u_float;
        u_float = 1.0f - u_float;
    }
    u_float = u_float * texture->width;
    int u;
    if(u_float - (int)u_float > 0.5f)
    {
        u = (int)u_float + 1;
    }else
    {
        u = (int)u_float;
    }
    u = u % texture->width;

    float v_float = texcoord[1];
    if(v_float < 0.0f)
    {
        v_float = -v_float;
    }
    if(v_float > 1.0f)
    {
        v_float = v_float - (int)v_float;
        v_float = 1.0f - v_float;
    }
    v_float = 1.0f - v_float;
    v_float = v_float * texture->height;

    int v;
    if(v_float - (int)v_float > 0.5f)
    {
        v = (int)v_float + 1;
    }else
    {
        v = (int)v_float;
    }
    v = v % texture->height;

    int row_len = texture->width * 3; // NOTE: is 3 gonna be okay?
    unsigned char* texel_ptr = texture->data + v * row_len + u * 3;
    vec3_assign(color_out, texel_ptr[0]/255.0f, texel_ptr[1]/255.0f, texel_ptr[2]/255.0f);
    return true;
}
