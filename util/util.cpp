#include <stdio.h>
#include "util.h"

bool getNextTokenInFile(char buffer[], FILE* fp)
{
    char c = fgetc(fp);
    while((c ==  ' ' || c == '\n') && c != EOF)
    {
        c = fgetc(fp);
    }
    if(c != EOF)
    {
        int i;
        for(i = 0 ; c != ' ' && c != '\n' && c != EOF; i++)
        {
            buffer[i] = c;
            c = fgetc(fp);
        }
        buffer[i] = '\0';
        return true;
    }else
    {
        return false;
    }
}

extern void toneMap(vec3, const vec3);
void genImageFromColorBuffer(unsigned char* image,
                             const float* color_buffer, const int num_pixels, const int num_samples)
{
    for(int i = 0; i < num_pixels; i++)
    {
        int index = i * 3;
        vec3 color;
        color[0] = color_buffer[index];
        color[1] = color_buffer[index+1];
        color[2] = color_buffer[index+2];
        vec3_scale(color, color, 1.0f / (float)num_samples);
        toneMap(color, color);

        image[index] = (unsigned char)(color[0] * 255.0f);
        image[index+1] = (unsigned char)(color[1] * 255.0f);
        image[index+2] = (unsigned char)(color[2] * 255.0f);
    }
}

