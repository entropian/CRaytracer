#pragma once
#include "util/vec.h"
#include "util/constants.h"

typedef struct Texture_s
{
    int width;
    int height;
    int comp;
    unsigned char* data;
    bool is_float;
    //char name[MAX_NAME_LENGTH];
}Texture;

bool loadTexture(Texture* tex, const char* file_name);
void freeTexture(Texture* tex);
bool getTexColor(vec3 color_out, const Texture* texture, const vec2 texcoord);
