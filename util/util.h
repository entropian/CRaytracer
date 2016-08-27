#pragma once

#include <stdio.h>
#include "vec.h"

bool getNextTokenInFile(char buffer[], FILE* fp);

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


inline void printVec3WithText(const char* text, const vec3 v)
{
    printf("%s %f, %f, %f\n", text, v[0], v[1], v[2]);
}

