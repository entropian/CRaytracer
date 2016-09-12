#include <cmath>
#include <cassert>
#include "../util/constants.h"
#include "../util/vec.h"
#include "../util/intvector.h"
#include "../aabb.h"
#include "../shapes/objecttype.h"
#include "../shapes/shapes.h"
#include "../shapes/instanced.h"
#include <xmmintrin.h>

// NOTES:
// Instead of if statement, trying doing both cases?
__m128 rayIntersectAABB4(const float bbox[24], const Ray ray)
{
    // TODO: change it so that this is done once per ray    
    __m128 a = _mm_set1_ps(1.0f / ray.direction[0]);
    __m128 b = _mm_set1_ps(1.0f / ray.direction[1]);
    __m128 c = _mm_set1_ps(1.0f / ray.direction[2]);

    __m128 ray_ox = _mm_set1_ps(ray.origin[0]);
    __m128 ray_oy = _mm_set1_ps(ray.origin[1]);
    __m128 ray_oz = _mm_set1_ps(ray.origin[2]);
    
    __m128 min_x = _mm_load_ps(&(bbox[0]));
    __m128 min_y = _mm_load_ps(&(bbox[4]));
    __m128 min_z = _mm_load_ps(&(bbox[8]));

    __m128 max_x = _mm_load_ps(&(bbox[12]));
    __m128 max_y = _mm_load_ps(&(bbox[16]));
    __m128 max_z = _mm_load_ps(&(bbox[20]));

    __m128 tx_min, tx_max;
    if(ray.direction[0] >= 0.0f)
    {
        tx_min = _mm_mul_ps(_mm_sub_ps(min_x, ray_ox), a);
        tx_max = _mm_mul_ps(_mm_sub_ps(max_x, ray_ox), a);
    }else
    {
        tx_min = _mm_mul_ps(_mm_sub_ps(max_x, ray_ox), a);
        tx_max = _mm_mul_ps(_mm_sub_ps(min_x, ray_ox), a);
    }

    __m128 ty_min, ty_max;
    if(ray.direction[1] >= 0.0f)
    {
        ty_min = _mm_mul_ps(_mm_sub_ps(min_y, ray_oy), b);
        ty_max = _mm_mul_ps(_mm_sub_ps(max_y, ray_oy), b);
    }else
    {
        ty_min = _mm_mul_ps(_mm_sub_ps(max_y, ray_oy), b);
        ty_max = _mm_mul_ps(_mm_sub_ps(min_y, ray_oy), b);
    }

    __m128 tz_min, tz_max;
    if(ray.direction[2] >= 0.0f)
    {
        tz_min = _mm_mul_ps(_mm_sub_ps(min_z, ray_oz), c);
        tz_max = _mm_mul_ps(_mm_sub_ps(max_z, ray_oz), c);            
    }else
    {
        tz_min = _mm_mul_ps(_mm_sub_ps(max_z, ray_oz), c);
        tz_max = _mm_mul_ps(_mm_sub_ps(min_z, ray_oz), c);            
    }

    __m128 t0, t1;
    t0 = _mm_max_ps(_mm_max_ps(tx_min, ty_min), tz_min);
    t1 = _mm_min_ps(_mm_min_ps(tx_max, ty_max), tz_max);

    __m128 epsilon = _mm_set1_ps(K_EPSILON);
    __m128 comp_mask = _mm_cmplt_ps(t0, t1);
    __m128 epsilon_mask = _mm_cmpgt_ps(t1, epsilon);
    __m128 hit_mask = _mm_and_ps(comp_mask, epsilon_mask);

    // Mix t0 and t1
    epsilon_mask = _mm_cmpgt_ps(t0, epsilon);
    t0 = _mm_and_ps(t0, epsilon_mask);
    t1 = _mm_andnot_ps(epsilon_mask, t1);
    __m128 t = _mm_or_ps(t0, t1);

    t = _mm_and_ps(t, hit_mask);    
    __m128 tmax = _mm_set1_ps(TMAX);
    tmax = _mm_andnot_ps(hit_mask, tmax);
    t = _mm_or_ps(t, tmax);
    
    return t;
}

typedef struct BVHNode4_s
{
    float bbox[2*3*4];

}BVHNode4;
