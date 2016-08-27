#pragma once
#include "util/vec.h"

typedef struct Texture_s
{
    int width;
    int height;
    int comp;
    unsigned char* data;    
}Texture;

bool loadTexture(Texture* tex, const char* file_name);
void freeTexture(Texture* tex);
bool getTexColor(vec3 color_out, const Texture* texture, const vec2 texcoord);
