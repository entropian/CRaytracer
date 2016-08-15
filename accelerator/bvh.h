#pragma once
#include <cmath>
#include <cassert>
#include "../util/constants.h"
#include "../util/vec.h"
#include "../util/intvector.h"
#include "../aabb.h"
#include "../shapes/objecttype.h"
#include "../shapes/shapes.h"
#include "../shapes/instanced.h"

static const int QUICKSORT_THRESHOLD = 30;

enum BVHNodeType
{
    NODE,
    LEAF
};

typedef struct BVHNode
{
    BVHNodeType type;
    int num_obj;
    BVHNode *left, *right;    
    Object_t obj;
    AABB aabb;
} BVHNode;

AABB calcBoundingVolume(const Object_t objects[], const int num_obj)
{
    AABB r, tmp_aabb;
    vec3_assign(r.min, FLT_MAX, FLT_MAX, FLT_MAX);
    vec3_assign(r.max, -HUGEVALUE, -HUGEVALUE, -HUGEVALUE);    
    for(int i = 0; i < num_obj; i++)
    {
        getObjectAABB(&tmp_aabb, objects[i]);
        addToAABB(&r, &tmp_aabb);
    }
    return r;
}

bool bvhExchangeElements(float centroid[], Object_t objects[], const int a, const int b, const int num_obj)
{
    if(a >= num_obj || b >= num_obj){
        fprintf(stderr, "Invalid index.\n");
        return false;
    }
    float tmp_centroid = centroid[a];
    Object_t tmp_obj = objects[a];
    centroid[a] = centroid[b];
    objects[a] = objects[b];
    centroid[b] = tmp_centroid;
    objects[b] = tmp_obj;
    return true;
}

int bvhQuicksortPartition(float centroid[], Object_t objects[], const int num_obj)
{
    int rand_index = rand() % num_obj;
    bvhExchangeElements(centroid, objects, rand_index, num_obj-1, num_obj);
    
    int i = -1;
    float pivot = centroid[num_obj-1];
    for(int j = 0; j < num_obj - 1; j++)
    {
        if(centroid[j] <= pivot)
        {
            i++;
            bvhExchangeElements(centroid, objects, j, i, num_obj);
        }
    }
    bvhExchangeElements(centroid, objects, i+1, num_obj-1, num_obj);
    return i+1;
}

void bvhQuicksort(float centroid[], Object_t objects[], const int num_obj)
{
    if(num_obj > 1)
    {
        int pivot_index = bvhQuicksortPartition(centroid, objects, num_obj);
        bvhQuicksort(centroid, objects, pivot_index);
        bvhQuicksort(&(centroid[pivot_index]), &(objects[pivot_index]), num_obj - pivot_index);
    }
}

int partitionObjects(Object_t objects[], int num_obj, const int axis_index)
{
    float* centroid = (float*)malloc(sizeof(float) * num_obj);    
    for(int i = 0; i < num_obj; ++i)
    {
        AABB tmp_aabb;
        getObjectAABB(&tmp_aabb, objects[i]);
        centroid[i] = (tmp_aabb.min[axis_index] + tmp_aabb.max[axis_index]) * 0.5f;
    }
    if(num_obj > QUICKSORT_THRESHOLD)
    {
        bvhQuicksort(centroid, objects, num_obj);
    }else
    {
        for(int i = 1; i < num_obj; ++i)
        {
            for(int j = i; j > 0 && centroid[j-1] > centroid[j]; --j)
            {
                float tmp = centroid[j];
                centroid[j] = centroid[j-1];
                centroid[j-1] = tmp;

                Object_t tmp_obj = objects[j];
                objects[j] = objects[j-1];
                objects[j-1] = tmp_obj;
            }
        }
    }
    free(centroid);
    return num_obj / 2;
}

void BVH_build(BVHNode **tree, Object_t objects[], int num_obj)
{
    assert(num_obj > 0);

    const int MIN_OBJECTS_PER_LEAF = 1;
    BVHNode* pNode = (BVHNode*)malloc(sizeof(BVHNode));
    *tree = pNode;
    pNode->aabb = calcBoundingVolume(objects, num_obj);
    float x_extent = pNode->aabb.max[0] - pNode->aabb.min[0];
    float y_extent = pNode->aabb.max[1] - pNode->aabb.min[1];
    float z_extent = pNode->aabb.max[2] - pNode->aabb.min[2];
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

    if(num_obj <= MIN_OBJECTS_PER_LEAF)
    {
        pNode->type = LEAF;
        pNode->num_obj = num_obj;
        pNode->obj = objects[0];
        pNode->left = NULL;
        pNode->right = NULL;
    }else
    {
        pNode->type = NODE;
        int k = partitionObjects(objects, num_obj, axis_index);
        BVH_build(&(pNode->left), objects, k);
        BVH_build(&(pNode->right), &(objects[k]), num_obj - k);        
    }
}

void BVH_destroy(BVHNode *tree)
{
    if(tree->type == NODE)
    {
        BVH_destroy(tree->left);
        BVH_destroy(tree->right);
        free(tree->left);
        free(tree->right);
    }
}

void printBVH(BVHNode* tree, int* leaf_count, int depth)
{
    for(int i = 0; i < depth; ++i)
    {
        printf("\t");
    }
    if(tree->type == NODE)
    {
        ++(*leaf_count);
        printf("NODE %d depth %d\n", *leaf_count, depth);
    }else if(tree->type == LEAF)
    {
        ++(*leaf_count);
        printf("LEAF %d depth %d ", *leaf_count, depth);
        printObjType(tree->obj.type);
    }

    if(tree->left){printBVH(tree->left, leaf_count, depth+1);}
    if(tree->right){printBVH(tree->right, leaf_count, depth+1);}
}

