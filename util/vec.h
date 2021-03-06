#pragma once

#include <math.h>

// TODO: figure out this stuff
static const double PI = 3.14159265358979323846264338327950288;
static const double EPS = 1e-8;
static const double EPS2 = EPS*EPS;
static const double EPS3 = EPS*EPS*EPS;
static const double INV_PI = 1.0 / PI;

typedef float vec2[2];
typedef float vec3[3];
typedef float vec4[4];

inline void vec3FromVec2(vec3 r, const vec2 a, const float b)
{
    r[0] = a[0];
    r[1] = a[1];
    r[2] = b;
}

inline void vec4FromVec3(vec4 r, const vec3 a, const float b)
{
    r[0] = a[0];
    r[1] = a[1];
    r[2] = a[2];
    r[3] = b;
}

inline void vec2FromVec3(vec2 r, const vec3 a)
{
    r[0] = a[0];
    r[1] = a[1];
}

inline void vec3FromVec4(vec3 r, const vec4 a)
{
    r[0] = a[0];
    r[1] = a[1];
    r[2] = a[2];
}

inline void vec3_add_c(vec3 r, const vec3 a, const float c)
{
    r[0] = a[0] + c;
    r[1] = a[1] + c;
    r[2] = a[2] + c;
}

inline void vec2_add(vec2 r, vec2 a, vec2 b)
{
    r[0] = a[0] + b[0];    
    r[1] = a[1] + b[1];    
}

inline void vec3_add(vec3 r, const vec3 a, const vec3 b)
{
    r[0] = a[0] + b[0];    
    r[1] = a[1] + b[1];    
    r[2] = a[2] + b[2];    
}

inline void vec4_add(vec4 r, const vec4 a, const vec4 b)
{
    r[0] = a[0] + b[0];    
    r[1] = a[1] + b[1];    
    r[2] = a[2] + b[2];
    r[3] = a[3] + b[3];        
}

inline void vec3_sub_c(vec3 r, const vec3 a, const float c)
{
    r[0] = a[0] - c;
    r[1] = a[1] - c;
    r[2] = a[2] - c;
}

inline void vec2_sub(vec2 r, const vec2 a, const vec2 b)
{
    r[0] = a[0] - b[0];    
    r[1] = a[1] - b[1];    
}

inline void vec3_sub(vec3 r, const vec3 a, const vec3 b)
{
    r[0] = a[0] - b[0];    
    r[1] = a[1] - b[1];    
    r[2] = a[2] - b[2];
}

inline void vec4_sub(vec4 r, const vec4 a, const vec4 b)
{
    r[0] = a[0] - b[0];    
    r[1] = a[1] - b[1];    
    r[2] = a[2] - b[2];
    r[3] = a[3] - b[3];        
}

inline void vec2_scale(vec2 r, const vec2 a, const float scale)
{
    r[0] = a[0] * scale;
    r[1] = a[1] * scale;
}

inline void vec3_scale(vec3 r, const vec3 a, const float scale)
{
    r[0] = a[0] * scale;
    r[1] = a[1] * scale;
    r[2] = a[2] * scale;
}

inline void vec4_scale(vec4 r, const vec4 a, const float scale)
{
    r[0] = a[0] * scale;
    r[1] = a[1] * scale;
    r[2] = a[2] * scale;
    r[3] = a[3] * scale;    
}

inline void vec2_mult(vec2 r, const vec2 a, const vec2 b)
{
    r[0] = a[0] * b[0];
    r[1] = a[1] * b[1];
}

inline void vec3_mult(vec3 r, const vec3 a, const vec3 b)
{
    r[0] = a[0] * b[0];
    r[1] = a[1] * b[1];
    r[2] = a[2] * b[2];
}

inline void vec4_mult(vec4 r, const vec4 a, const vec4 b)
{
    r[0] = a[0] * b[0];
    r[1] = a[1] * b[1];
    r[2] = a[2] * b[2];
    r[3] = a[3] * b[3];    
}

inline void vec3_div(vec3 r, const vec3 a, const vec3 b)
{
    r[0] = a[0] / b[0];
    r[1] = a[1] / b[1];
    r[2] = a[2] / b[2];
}

inline float vec2_dot(const vec2 a, const vec2 b)
{
    return a[0]*b[0] + a[1]*b[1];
}

inline float vec3_dot(const vec3 a, const vec3 b)
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

inline float vec4_dot(const vec4 a, const vec4 b)
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + a[3]*b[3];
}


inline float vec2_dot(const vec4 a, const float x, const float y)
{
    return a[0]*x + a[1]*y;
}

inline float vec3_dot(const vec4 a, const float x, const float y, const float z)
{
    return a[0]*x + a[1]*y + a[2]*z;
}

inline float vec4_dot(const vec4 a, const float x, const float y, const float z, const float w)
{
    return a[0]*x + a[1]*y + a[2]*z + a[3]*w;
}

inline void vec2_clear(vec2 a)
{
    a[0] = a[1] = 0.0f;
}

inline void vec3_clear(vec3 a)
{
    a[0] = a[1] = a[2] = 0.0f;
}

inline void vec4_clear(vec4 a)
{
    a[0] = a[1] = a[2] = a[3] = 0.0f;
}

inline void vec2_normalize(vec2 r, const vec2 a)
{
    float length = (float)sqrt(vec2_dot(a, a));
    if(length == 0)
    {
        vec2_clear(r);
        return;
    }
    vec2_scale(r, a, 1.0f/length);
}

inline void vec3_normalize(vec3 r, const vec3 a)
{
    float length = (float)sqrt(vec3_dot(a, a));
    if(length == 0)
    {
        vec3_clear(r);
        return;
    }
    vec3_scale(r, a, 1.0f/length);
}

inline void vec4_normalize(vec4 r, vec4 a)
{
    float length = (float)sqrt(vec4_dot(a, a));
    if(length == 0)
    {
        vec4_clear(r);
        return;
    }
    vec4_scale(r, a, 1.0f/length);
}

inline void vec3_cross(vec3 r, const vec3 a, const vec3 b)
{
    r[0] = a[1]*b[2] - a[2]*b[1];
    r[1] = a[2]*b[0] - a[0]*b[2];
    r[2] = a[0]*b[1] - a[1]*b[0];
}

inline float vec2_length(const vec2 a)
{
    return (float)sqrt(vec2_dot(a, a));
}

inline float vec3_length(const vec3 a)
{
    return (float)sqrt(vec3_dot(a, a));
}

inline float vec4_length(const vec4 a)
{
    return (float)sqrt(vec4_dot(a, a));
}

inline void vec2_copy(vec2 r, const vec2 a)
{
    r[0] = a[0];
    r[1] = a[1];
}

inline void vec3_copy(vec3 r, const vec3 a)
{
    r[0] = a[0];
    r[1] = a[1];
    r[2] = a[2];
}

inline void vec4_copy(vec4 r, const vec4 a)
{
    r[0] = a[0];
    r[1] = a[1];
    r[2] = a[2];
    r[3] = a[3];    
}

inline void vec2_assign(vec2 r, const float x, const float y)
{
    r[0] = x;
    r[1] = y;
}

inline void vec3_assign(vec3 r, const float x, const float y, const float z)
{
    r[0] = x;
    r[1] = y;
    r[2] = z;
}

inline void vec4_assign(vec4 r, const float x, const float y, const float z, const float a)
{
    r[0] = x;
    r[1] = y;
    r[2] = z;
    r[3] = a;
}

inline void vec2_negate(vec2 r, const vec2 a)
{
    r[0] = -a[0];
    r[1] = -a[1];
}

inline void vec3_negate(vec3 r, const vec3 a)
{
    r[0] = -a[0];
    r[1] = -a[1];
    r[2] = -a[2];
}

inline void vec4_negate(vec4 r, const vec4 a)
{
    r[0] = -a[0];
    r[1] = -a[1];
    r[2] = -a[2];
    r[3] = -a[3];    
}

inline void vec2_pow(vec2 r, const vec2 a, const float b)
{
    r[0] = powf(a[0], b);
    r[1] = powf(a[1], b);
}

inline void vec3_pow(vec3 r, const vec3 a, const float b)
{
    r[0] = powf(a[0], b);
    r[1] = powf(a[1], b);
    r[2] = powf(a[2], b);
}

inline void vec4_pow(vec4 r, const vec4 a, const float b)
{
    r[0] = powf(a[0], b);
    r[1] = powf(a[1], b);
    r[2] = powf(a[2], b);
    r[3] = powf(a[3], b);    
}

inline bool vec2_equal(const vec2 a, const vec2 b)
{
    if(a[0] != b[0]){return false;}
    if(a[1] != b[1]){return false;}
    return true;
}

inline bool vec3_equal(const vec3 a, const vec3 b)
{
    if(a[0] != b[0]){return false;}
    if(a[1] != b[1]){return false;}
    if(a[2] != b[2]){return false;}
    return true;
}

inline bool vec4_equal(const vec4 a, const vec4 b)
{
    if(a[0] != b[0]){return false;}
    if(a[1] != b[1]){return false;}
    if(a[2] != b[2]){return false;}
    if(a[3] != b[3]){return false;}    
    return true;
}

inline bool vec3_less(const vec3 a, const vec3 b)
{
    if(a[0] < b[0] && a[1] < b[1] && a[2] < b[2])
    {
        return true;
    }
    return false;
}

inline void vec3_sqrt(vec3 r, const vec3 a)
{
    r[0] = sqrtf(a[0]);
    r[1] = sqrtf(a[1]);
    r[2] = sqrtf(a[2]);
}

