#pragma once
#include <xmmintrin.h>
#include <emmintrin.h>


typedef struct __declspec(align(16)) vec3_4_s
{
    __m128 x, y, z;
}vec3_4;


inline void vec3_4_assign(vec3_4* vec,const float x[4], const float y[4], const float z[4])
{
    vec->x = _mm_load_ps(x);
    vec->y = _mm_load_ps(y);
    vec->z = _mm_load_ps(z);    
}

inline void vec3_4_copy(vec3_4* dst, const vec3_4* src)
{
    dst->x = src->x;
    dst->y = src->y;
    dst->z = src->z;    
}

inline void vec3_4_add(vec3_4* dst, const vec3_4* a, const vec3_4* b)
{
    dst->x = _mm_add_ps(a->x, b->x);
    dst->y = _mm_add_ps(a->y, b->y);
    dst->z = _mm_add_ps(a->z, b->z);    
}

inline void vec3_4_sub(vec3_4* dst, const vec3_4* a, const vec3_4* b)
{
    dst->x = _mm_sub_ps(a->x, b->x);
    dst->y = _mm_sub_ps(a->y, b->y);
    dst->z = _mm_sub_ps(a->z, b->z);    
}

inline void vec3_4_mult(vec3_4* dst, const vec3_4* a, const vec3_4* b)
{
    dst->x = _mm_mul_ps(a->x, b->x);
    dst->y = _mm_mul_ps(a->y, b->y);
    dst->z = _mm_mul_ps(a->z, b->z);    
}

inline void vec3_4_div(vec3_4* dst, const vec3_4* a, const vec3_4* b)
{
    dst->x = _mm_div_ps(a->x, b->x);
    dst->y = _mm_div_ps(a->y, b->y);
    dst->z = _mm_div_ps(a->z, b->z);        
}
