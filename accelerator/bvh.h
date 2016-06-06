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

int partitionObjects(Object_t objects[], int num_obj, const int axis_index)
{
    float* centroid = (float*)malloc(sizeof(float) * num_obj);
    for(int i = 0; i < num_obj; ++i)
    {
        AABB tmp_aabb;
        getObjectAABB(&tmp_aabb, objects[i]);
        centroid[i] = (tmp_aabb.min[axis_index] + tmp_aabb.max[axis_index]) * 0.5f;
    }

    int sorted = 0;
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
    free(centroid);
    return num_obj / 2;
}

void BVH_build(BVHNode **tree, Object_t objects[], int num_obj, const int axis_index)
{
    assert(num_obj > 0);

    const int MIN_OBJECTS_PER_LEAF = 1;
    BVHNode* pNode = (BVHNode*)malloc(sizeof(BVHNode));
    *tree = pNode;
    pNode->aabb = calcBoundingVolume(objects, num_obj);

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
        BVH_build(&(pNode->left), objects, k, (axis_index + 1) % 3);
        BVH_build(&(pNode->right), &(objects[k]), num_obj - k, (axis_index + 1) % 3);        
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

