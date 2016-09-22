#include <cmath>
#include <cassert>
#include "../util/constants.h"
#include "../util/vec.h"
#include "../shapes/objecttype.h"
#include "../shapes/shapes.h"
#include "../shapes/instanced.h"
#include "../aabb.h"
#include <xmmintrin.h>
#include "bvh.h"

/*
  This QBVH implementation is based on the paper
  Shallow Bounding Volume Hierarchies for Fast SIMD Ray Tracing of Incoherent Rays
  by Dammertz, Hanika, and Keller
 */

#define CACHE_ALIGN __declspec(align(16))

__m128 rayIntersectAABB4(__m128* less_than, const float bbox[24], const Ray ray)
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
    *less_than = _mm_cmplt_ps(t0, epsilon);
    t0 = _mm_and_ps(t0, epsilon_mask);
    t1 = _mm_andnot_ps(epsilon_mask, t1);
    __m128 t = _mm_or_ps(t0, t1);

    t = _mm_and_ps(t, hit_mask);    
    __m128 tmax = _mm_set1_ps(TMAX);
    tmax = _mm_andnot_ps(hit_mask, tmax);
    t = _mm_or_ps(t, tmax);
    
    return t;
}

typedef _declspec(align(16)) struct BVHNode4_s
{
    float bbox[2*3*4]; // 96 bytes
    void *child[4]; // 16 bytes
    char type[4]; // 4 bytes
    int axis0, axis1, axis2; // 12 bytes
    // Total 128 bytes
}BVHNode4;

void BVH4_build(BVHNode4 **tree, Object_t objects[], int num_obj)
{
    assert(num_obj > 0);

    const int MIN_OBJECTS_PER_LEAF = 4;
    BVHNode4* pNode = (BVHNode4*)_aligned_malloc(sizeof(BVHNode4) + 1, 16);
    *tree = pNode;
    // Partion objects into two sets, then partition each of those two sets to produce four sets
    // Calculate the bounding box for each of the four sets
    AABB tmp_aabb = calcBoundingVolume(objects, num_obj);
    float x_extent = tmp_aabb.max[0] - tmp_aabb.min[0];
    float y_extent = tmp_aabb.max[1] - tmp_aabb.min[1];
    float z_extent = tmp_aabb.max[2] - tmp_aabb.min[2];    
    int axis_index;
    if(x_extent > y_extent && x_extent > z_extent)
    {
        axis_index = 0;
    }else if(y_extent > z_extent)
    {
        axis_index = 1;
    }else
    {
        axis_index = 2;
    }
    pNode->axis0 = axis_index;
    int middle_split = partitionObjects(objects, num_obj, axis_index);
    
    tmp_aabb = calcBoundingVolume(objects, middle_split);
    x_extent = tmp_aabb.max[0] - tmp_aabb.min[0];
    y_extent = tmp_aabb.max[1] - tmp_aabb.min[1];
    z_extent = tmp_aabb.max[2] - tmp_aabb.min[2];    
    if(x_extent > y_extent && x_extent > z_extent)
    {
        axis_index = 0;
    }else if(y_extent > z_extent)
    {
        axis_index = 1;
    }else
    {
        axis_index = 2;
    }
    pNode->axis1 = axis_index;
    int left_split = partitionObjects(objects, middle_split, axis_index);

    tmp_aabb = calcBoundingVolume(&(objects[left_split]), num_obj - middle_split);
    x_extent = tmp_aabb.max[0] - tmp_aabb.min[0];
    y_extent = tmp_aabb.max[1] - tmp_aabb.min[1];
    z_extent = tmp_aabb.max[2] - tmp_aabb.min[2];    
    if(x_extent > y_extent && x_extent > z_extent)
    {
        axis_index = 0;
    }else if(y_extent > z_extent)
    {
        axis_index = 1;
    }else
    {
        axis_index = 2;
    }
    pNode->axis2 = axis_index;
    int right_split = partitionObjects(&(objects[middle_split]), num_obj - middle_split, axis_index);
    right_split += middle_split;

    int index, length;

    // First AABB
    index = 0;
    length = left_split;
    tmp_aabb = calcBoundingVolume(objects, left_split);
    
    int offset = 0;
    pNode->bbox[offset + 0] = tmp_aabb.min[0];
    pNode->bbox[offset + 4] = tmp_aabb.min[1];
    pNode->bbox[offset + 8] = tmp_aabb.min[2];
    pNode->bbox[offset + 12] = tmp_aabb.max[0];
    pNode->bbox[offset + 16] = tmp_aabb.max[1];
    pNode->bbox[offset + 20] = tmp_aabb.max[2];

    // Second AABB
    index = left_split;
    length = middle_split - left_split;
    tmp_aabb = calcBoundingVolume(&(objects[left_split]), middle_split - left_split);

    offset = 1;
    pNode->bbox[offset + 0] = tmp_aabb.min[0];
    pNode->bbox[offset + 4] = tmp_aabb.min[1];
    pNode->bbox[offset + 8] = tmp_aabb.min[2];
    pNode->bbox[offset + 12] = tmp_aabb.max[0];
    pNode->bbox[offset + 16] = tmp_aabb.max[1];
    pNode->bbox[offset + 20] = tmp_aabb.max[2];

    // Third AABB
    index = middle_split;
    length = right_split - middle_split;
    tmp_aabb = calcBoundingVolume(&(objects[middle_split]), right_split - middle_split);    

    offset = 2;
    pNode->bbox[offset + 0] = tmp_aabb.min[0];
    pNode->bbox[offset + 4] = tmp_aabb.min[1];
    pNode->bbox[offset + 8] = tmp_aabb.min[2];
    pNode->bbox[offset + 12] = tmp_aabb.max[0];
    pNode->bbox[offset + 16] = tmp_aabb.max[1];
    pNode->bbox[offset + 20] = tmp_aabb.max[2];

    // Fourth AABB
    index = right_split;
    length = num_obj - right_split;
    tmp_aabb = calcBoundingVolume(&(objects[right_split]), num_obj - right_split);

    offset = 3;
    pNode->bbox[offset + 0] = tmp_aabb.min[0];
    pNode->bbox[offset + 4] = tmp_aabb.min[1];
    pNode->bbox[offset + 8] = tmp_aabb.min[2];
    pNode->bbox[offset + 12] = tmp_aabb.max[0];
    pNode->bbox[offset + 16] = tmp_aabb.max[1];
    pNode->bbox[offset + 20] = tmp_aabb.max[2];        
    

    // Check the number of objects in each set
    // if <= 4, make child a leaf
    // For leafs, since the relevant portion of the object array is already sorted,
    // child[i] can point to the first object, and type[i] be the number of objects in this leaf
    if(left_split > 4)
    {
        pNode->type[0] = 0;
        BVH4_build((BVHNode4**)&(pNode->child[0]), objects, left_split);
    }else
    {
        pNode->child[0] = objects;
        pNode->type[0] = left_split;
    }

    if((middle_split - left_split) > 4)
    {
        pNode->type[1] = 0;        
        BVH4_build((BVHNode4**)&(pNode->child[1]), &(objects[left_split]), middle_split - left_split);
    }else
    {
        pNode->child[1] = &(objects[left_split]);
        pNode->type[1] = middle_split - left_split;
    }

    if((right_split - middle_split) > 4)
    {
        pNode->type[2] = 0;        
        BVH4_build((BVHNode4**)&(pNode->child[2]), &(objects[middle_split]), right_split - middle_split);
    }else
    {
        pNode->child[2] = &(objects[middle_split]);
        pNode->type[2] = right_split - middle_split;
    }

    if((num_obj - right_split) > 4)
    {
        pNode->type[3] = 0;        
        BVH4_build((BVHNode4**)&(pNode->child[3]), &(objects[right_split]), num_obj - right_split);
    }else
    {
        pNode->child[3] = &(objects[right_split]);
        pNode->type[3] = num_obj - right_split;        
    }            
}

float BVH4IntersectTest(ShadeRec* sr, const BVHNode4* tree, const Ray ray)
{
    CACHE_ALIGN float bbox_result[5] = {TMAX, TMAX, TMAX, TMAX, TMAX};
    CACHE_ALIGN float less_than[4];
    __m128* tmp_ptr = (__m128*)less_than;
    __m128 tmp = rayIntersectAABB4(tmp_ptr, tree->bbox, ray);
    _mm_store_ps(bbox_result, tmp);
    int indices[4];
    
    if(ray.direction[tree->axis0] > 0.0f)
    {
        if(ray.direction[tree->axis1] > 0.0f)
        {
            indices[0] = 0;
            indices[1] = 1;
        }else
        {
            indices[0] = 1;
            indices[1] = 0;
        }
        if(ray.direction[tree->axis2] > 0.0f)
        {
            indices[2] = 2;
            indices[3] = 3;
        }else
        {
            indices[2] = 3;
            indices[3] = 2;
        }
    }else
    {
        if(ray.direction[tree->axis2] > 0.0f)
        {
            indices[0] = 2;
            indices[1] = 3;
        }else
        {
            indices[0] = 3;
            indices[1] = 2;
        }
        if(ray.direction[tree->axis1] > 0.0f)
        {
            indices[2] = 0;
            indices[3] = 1;
        }else
        {
            indices[2] = 1;
            indices[3] = 0;
        }
    }
   
    float min_t = TMAX;
    for(int i = 0; i < 4; i++)
    {        
        char index = indices[i];
<<<<<<< HEAD
        if(bbox_result[index] < min_t || (less_than[index] && bbox_result[index] < TMAX))
=======
        if(bbox_result[index] < min_t || (less_than[index] && bbox_result[index] < TMAX)) // key line?
>>>>>>> 4e9b9c8062f1b2e557427376eda04b5ca86400c9
        {
            float t;
            ShadeRec tmp_sr;
            if(tree->type[index] == 0)
            {
                t = BVH4IntersectTest(&tmp_sr, (BVHNode4*)(tree->child[index]), ray);
                if(t < min_t)
                {                    
                    min_t = t;
                    *sr = tmp_sr;
                }                
            }else // child[i] is a leaf
            {
                Object_t* obj_ptr = (Object_t*)(tree->child[index]);
                for(int j = 0; j < tree->type[index]; j++)
                {
                    t = rayIntersectObject(&tmp_sr, *obj_ptr, ray);
                    obj_ptr++;
                    if(t < min_t)
                    {
                        min_t = t;
                        *sr = tmp_sr;
                    }                    
                }
            }

        }
    }
    return min_t;
}


float BVH4ShadowIntersectTest(const BVHNode4* tree, const Ray ray, const float light_dist)
{
    CACHE_ALIGN float bbox_result[5] = {TMAX, TMAX, TMAX, TMAX, TMAX};
    CACHE_ALIGN float less_than[4];    
    __m128* tmp_ptr = (__m128*)less_than;
    __m128 tmp = rayIntersectAABB4(tmp_ptr, tree->bbox, ray);    
    _mm_store_ps(bbox_result, tmp);
    int indices[4];
    
    if(ray.direction[tree->axis0] > 0.0f)
    {
        if(ray.direction[tree->axis1] > 0.0f)
        {
            indices[0] = 0;
            indices[1] = 1;
        }else
        {
            indices[0] = 1;
            indices[1] = 0;
        }
        if(ray.direction[tree->axis2] > 0.0f)
        {
            indices[2] = 2;
            indices[3] = 3;
        }else
        {
            indices[2] = 3;
            indices[3] = 2;
        }
    }else
    {
        if(ray.direction[tree->axis2] > 0.0f)
        {
            indices[0] = 2;
            indices[1] = 3;
        }else
        {
            indices[0] = 3;
            indices[1] = 2;
        }
        if(ray.direction[tree->axis1] > 0.0f)
        {
            indices[2] = 0;
            indices[3] = 1;
        }else
        {
            indices[2] = 1;
            indices[3] = 0;
        }
    }
    
    for(int i = 0; i < 4; i++)
    {
        char index = indices[i];
<<<<<<< HEAD
        if(bbox_result[index] < light_dist || (less_than[index] && bbox_result[index] < TMAX))
=======
        if(bbox_result[index] < light_dist || (less_than[index] && bbox_result[index] < TMAX)) // Key line?
>>>>>>> 4e9b9c8062f1b2e557427376eda04b5ca86400c9
        {
            float t = TMAX;
            if(tree->type[index] == 0)
            {
                t = BVH4ShadowIntersectTest((BVHNode4*)(tree->child[index]), ray, light_dist);
                if(t < light_dist)
                {
                    return t;
                }                
            }else // child[i] is a leaf
            {
                Object_t* obj_ptr = (Object_t*)(tree->child[index]);
                for(int j = 0; j < tree->type[index]; j++)
                {
                    t = shadowRayIntersectObject(*obj_ptr, ray);
                    obj_ptr++;
                    if(t < light_dist)
                    {
                        return t;
                    }                    
                }
            }

        }
    }
    return TMAX;
}

