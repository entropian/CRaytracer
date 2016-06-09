#pragma once

#include "vec.h"

inline float min(float a, float b)
{
    return (a < b) ? a : b;
}

inline float max(float a, float b)
{
    return (a > b) ? a : b;
}

inline float alt_abs(const float a)
{
    return (a < 0.0f) ? -a : a;
}

inline float clamp(const float x, const float min, const float max)
{
    return (x < min ? min : (x > max ? max : x));
}

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

void printVec3WithText(const char* text, const vec3 v)
{
    printf("%s %f, %f, %f\n", text, v[0], v[1], v[2]);
}
