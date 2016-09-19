#include <cmath>
#include <cassert>
#include "../util/constants.h"
#include "../util/vec.h"
#include "../shapes/objecttype.h"
#include "../shapes/shapes.h"
#include "../shapes/instanced.h"
#include <xmmintrin.h>
#include "bvh.h"

#define CACHE_ALIGN __declspec(align(16))

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

typedef _declspec(align(16)) struct BVHNode4_s
{
    float bbox[2*3*4]; // 96 bytes
    //void *child1, *child2, *child3, *child4;
    void *child[4];
    //char type1, type2, type3, type4; // 0 = inner node; 1 = leaf
    char type[4];
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
    int right_split = partitionObjects(&(objects[middle_split]), num_obj - middle_split, axis_index);
    right_split += middle_split;

    /*
    char tabs[56];
    int i; 
    for(i = 0; i < depth; i++)
    {
        tabs[i] = '\t';
    }
    tabs[i] = '\0';
    */
    int index, length;

    // First AABB
    index = 0;
    length = left_split;
    tmp_aabb = calcBoundingVolume(objects, left_split);
    //fprintf(fp, "%s%f %f %f\n", tabs, tmp_aabb.min[0], tmp_aabb.min[1], tmp_aabb.min[2]);
    //fprintf(fp, "%s%f %f %f\n", tabs, tmp_aabb.max[0], tmp_aabb.max[1], tmp_aabb.max[2]);
    
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
    //fprintf(fp, "%s%f %f %f\n", tabs, tmp_aabb.min[0], tmp_aabb.min[1], tmp_aabb.min[2]);
    //fprintf(fp, "%s%f %f %f\n", tabs, tmp_aabb.max[0], tmp_aabb.max[1], tmp_aabb.max[2]);    
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
    //fprintf(fp, "%s%f %f %f\n", tabs, tmp_aabb.min[0], tmp_aabb.min[1], tmp_aabb.min[2]);
    //fprintf(fp, "%s%f %f %f\n", tabs, tmp_aabb.max[0], tmp_aabb.max[1], tmp_aabb.max[2]);
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
    //fprintf(fp, "%s%f %f %f\n", tabs, tmp_aabb.min[0], tmp_aabb.min[1], tmp_aabb.min[2]);
    //fprintf(fp, "%s%f %f %f\n", tabs, tmp_aabb.max[0], tmp_aabb.max[1], tmp_aabb.max[2]);
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
    __m128 tmp = rayIntersectAABB4(tree->bbox, ray);
    _mm_store_ps(bbox_result, tmp);
    float min_t = TMAX;
    for(int i = 0; i < 4; i++)
    {        
        if(bbox_result[i] < TMAX)
        {
            float t;
            ShadeRec tmp_sr;
            if(tree->type[i] == 0)
            {
                t = BVH4IntersectTest(&tmp_sr, (BVHNode4*)(tree->child[i]), ray);
                if(t < min_t)
                {
                    min_t = t;
                    *sr = tmp_sr;
                }                
            }else // child[i] is a leaf
            {
                Object_t* obj_ptr = (Object_t*)(tree->child[i]);
                for(int j = 0; j < tree->type[i]; j++)
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

float BVH4ShadowIntersectTest(const BVHNode4* tree, const Ray ray)
{
    CACHE_ALIGN float bbox_result[5] = {TMAX, TMAX, TMAX, TMAX, TMAX};
    __m128 tmp = rayIntersectAABB4(tree->bbox, ray);
    _mm_store_ps(bbox_result, tmp);
    for(int i = 0; i < 4; i++)
    {        
        if(bbox_result[i] < TMAX)
        {
            float t;
            if(tree->type[i] == 0)
            {
                t = BVH4ShadowIntersectTest((BVHNode4*)(tree->child[i]), ray);
                if(t < TMAX)
                {
                    return t;
                }                
            }else // child[i] is a leaf
            {
                Object_t* obj_ptr = (Object_t*)(tree->child[i]);
                for(int j = 0; j < tree->type[i]; j++)
                {
                    t = shadowRayIntersectObject(*obj_ptr, ray);
                    obj_ptr++;
                    if(t < TMAX)
                    {
                        return t;
                    }                    
                }
            }

        }
    }
    return TMAX;
}

/*
bool BVH4IntersectTest(const BVHNode4* tree, const Ray ray, const int depth)
{
    CACHE_ALIGN float t[5] = {TMAX, TMAX, TMAX, TMAX, TMAX};
    
    __m128 bbox_result = rayIntersectAABB4(tree->bbox, ray);
    _mm_store_ps(t, bbox_result);
    bool hit_status = false;
    for(int i = 0; i < 4 && !hit_status; i++)
    {        
        if(t[i] < TMAX)
        {
            if(tree->type[i] == 0)
            {
                if(BVH4IntersectTest((BVHNode4*)(tree->child[i]), ray, depth + 1))
                {
                    hit_status = true;
                }
            }else if(tree->type[i] == 1)
            {
                hit_status = true;
            }
        }
        if(hit_status)
        {
            return true;
        }
    }
    return false;
}
*/

void BVH4_print(const BVHNode4 *tree, const int depth, FILE *fp)
{
    char tabs[56];
    int i;
    for(i = 0; i < depth; i++)
    {
        tabs[i] = '\t';
    }
    tabs[i] = '\0';

    for(int i = 0; i < 4; i++)
    {
        fprintf(fp, "%s%f %f %f\n", tabs, tree->bbox[i], tree->bbox[i+4], tree->bbox[i+8]);
        fprintf(fp, "%s%f %f %f\n", tabs, tree->bbox[i+12], tree->bbox[i+16], tree->bbox[i+20]);        
    }

    for(int i = 0; i < 4; i++)
    {
        if(tree->type[i] == 0)
        {
            BVH4_print((BVHNode4*)(tree->child[i]), depth+1, fp);
        }
    }
}
