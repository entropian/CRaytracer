#pragma once

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

typedef struct Texture_s
{
    int width;
    int height;
    int comp;
    unsigned char* data;    
}Texture;

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

