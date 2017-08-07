#pragma once
#include <string.h>
#include <stdio.h>
#include "vec.h"
#include "constants.h"

bool getNextTokenInFile(char buffer[], FILE* fp);

inline float min(float a, float b)
{
    return (a < b) ? a : b;
}

inline float max(float a, float b)
{
    return (a > b) ? a : b;
}

inline float max3(float a, float b, float c)
{
    return max(a, max(b, c));
}
 
inline float min3(float a, float b, float c)
{
    return min(a, min(b, c));
}

inline float alt_abs(float a)
{
    return (a < 0.0f) ? -a : a;
}

inline float clamp(float x, float min, float max)
{
    return (x < min ? min : (x > max ? max : x));
}


inline void printVec3WithText(const char* text, const vec3 v)
{
    printf("%s %f, %f, %f\n", text, v[0], v[1], v[2]);
}

inline int isZeroVector(vec3 v)
{
    //if (fabs(v[0]) < K_EPSILON && fabs(v[1]) < K_EPSILON && fabs(v[2]) < K_EPSILON) return true;
    if (fabs(v[0]) < K_EPSILON && fabs(v[1]) < K_EPSILON && fabs(v[2]) < K_EPSILON) return true;   
    return false;
}

static void stringCopy(char* dest, const int max_len, const char* src)
{
#ifdef _MSC_VER
    strcpy_s(dest, max_len, src);
#else
    strcpy(dest, src);
#endif
}

static void stringNCopy(char* dest, const int max_len, const char*src, const int len)
{
#ifdef _MSC_VER
    strncpy_s(dest, max_len, src, len);
#else
    strncpy(dest, src, len);
#endif
}
