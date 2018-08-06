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

void fixBufferInfAndNaN(float* color_buffer, const int num_pixels, const int width)
{
    for(int i = 0; i < num_pixels; i++)
    {
        int index = i * 3;
        if(isinf(color_buffer[0]) || isnan(color_buffer[0]))
        {
            int count = 0;
            vec3 neighbor_sum = {0.0f, 0.0f, 0.0f};
            int current_index = index - width * 3;
            if(current_index > -1)
            {
                neighbor_sum[0] += color_buffer[current_index];
                neighbor_sum[1] += color_buffer[current_index];
                neighbor_sum[2] += color_buffer[current_index];
                count++;
            }

            current_index = index - 3;
            if(current_index > -1)
            {
                neighbor_sum[0] += color_buffer[current_index];
                neighbor_sum[1] += color_buffer[current_index];
                neighbor_sum[2] += color_buffer[current_index];
                count++;                
            }

            current_index = index + 3;
            if(current_index < num_pixels * 3)
            {
                neighbor_sum[0] += color_buffer[current_index];
                neighbor_sum[1] += color_buffer[current_index];
                neighbor_sum[2] += color_buffer[current_index];
                count++;                
            }

            current_index = index + width * 3;
            if(current_index < num_pixels * 3)
            {
                neighbor_sum[0] += color_buffer[current_index];
                neighbor_sum[1] += color_buffer[current_index];
                neighbor_sum[2] += color_buffer[current_index];
                count++;                
            }

            vec3 new_value;
            vec3_scale(new_value, neighbor_sum, 1.0f / count);
            color_buffer[index] = new_value[0];
            color_buffer[index+1] = new_value[1];
            color_buffer[index+2] = new_value[2];
        }
    }
}

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

