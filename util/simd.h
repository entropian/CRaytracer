#pragma once
#include <xmmintrin.h>
#include <emmintrin.h>

#ifdef _MSC_VER
typedef struct __declspec(align(16)) vec3_4_s
#else
typedef struct __attribute__((aligned(16))) vec3_4_s
#endif
{
    __m128 x, y, z;
}vec3_4;

#ifdef _MSC_VER
#define CACHE_ALIGN __declspec(align(16))
#else
#define CACHE_ALIGN __attribute__((aligned(16)))
#endif

inline void vec3_4_assign(vec3_4* vec,const float x[4], const float y[4], const float z[4])
{
    vec->x = _mm_load_ps(x);
    vec->y = _mm_load_ps(y);
    vec->z = _mm_load_ps(z);    
}

inline void vec3_4_assignv(vec3_4* vec, const vec3 v0, const vec3 v1, const vec3 v2, const vec3 v3)
{
    CACHE_ALIGN float buffer[4];
    buffer[0] = v0[0];
    buffer[1] = v1[0];
    buffer[2] = v2[0];
    buffer[3] = v3[0];
    vec->x = _mm_load_ps(buffer);

    buffer[0] = v0[1];
    buffer[1] = v1[1];
    buffer[2] = v2[1];
    buffer[3] = v3[1];
    vec->y = _mm_load_ps(buffer);

    buffer[0] = v0[2];
    buffer[1] = v1[2];
    buffer[2] = v2[2];
    buffer[3] = v3[2];
    vec->z = _mm_load_ps(buffer);
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
