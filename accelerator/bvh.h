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
    void* obj_ptr;
    ObjectType obj_type;
    BVHNode *left, *right;
    AABB aabb;
} BVHNode;

AABB calcBoundingVolume(const void* obj_ptrs[], const ObjectType obj_types[], const int num_obj)
{
    AABB r, tmp_aabb;
    vec3_assign(r.min, FLT_MAX, FLT_MAX, FLT_MAX);
    vec3_assign(r.max, -HUGEVALUE, -HUGEVALUE, -HUGEVALUE);    
    for(int i = 0; i < num_obj; i++)
    {
        getObjectAABB(&tmp_aabb, obj_ptrs[i], obj_types[i]);
        addToAABB(&r, &tmp_aabb);
    }
    return r;
}

int partitionObjects(void* obj_ptrs[], ObjectType obj_types[], int num_obj, const int axis_index)
{
    vec3 axis;
    switch(axis_index)
    {
    case 0:
        vec3_copy(axis, RIGHT);        
        break;
    case 1:
        vec3_copy(axis, UP);        
        break;
    case 2:
        vec3_copy(axis, BACKWARD);
        break;
    }

    float* centroid = (float*)malloc(sizeof(float) * num_obj);
    for(int i = 0; i < num_obj; ++i)
    {
        AABB tmp_aabb;
        getObjectAABB(&tmp_aabb, obj_ptrs[i], obj_types[i]);
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

            void* tmp_void = obj_ptrs[j];
            obj_ptrs[j] = obj_ptrs[j-1];
            obj_ptrs[j-1] = tmp_void;

            ObjectType tmp_type = obj_types[j];
            obj_types[j] = obj_types[j-1];
            obj_types[j-1] = tmp_type;
        }
    }
    free(centroid);
    return num_obj / 2;
}

void BVH_create(BVHNode **tree, void* obj_ptrs[], ObjectType obj_types[], int num_obj, const int axis_index)
{
    assert(num_obj > 0);

    const int MIN_OBJECTS_PER_LEAF = 1;
    BVHNode* pNode = (BVHNode*)malloc(sizeof(BVHNode));
    *tree = pNode;
    pNode->aabb = calcBoundingVolume((const void**)obj_ptrs, obj_types, num_obj);

    if(num_obj <= MIN_OBJECTS_PER_LEAF)
    {
        pNode->type = LEAF;
        pNode->num_obj = num_obj;
        pNode->obj_ptr = obj_ptrs[0];
        pNode->left = NULL;
        pNode->right = NULL;
    }else
    {
        pNode->type = NODE;
        int k = partitionObjects(obj_ptrs, obj_types, num_obj, axis_index);
        BVH_create(&(pNode->left), obj_ptrs, obj_types, k, (axis_index + 1) % 3);
        BVH_create(&(pNode->right), &(obj_ptrs[k]), &(obj_types[k]), num_obj - k, (axis_index + 1) % 3);
    }
}

void printBVH(BVHNode* tree, int* leaf_count, int depth)
{
    if(tree->type == NODE)
    {
        ++(*leaf_count);        
        printf("NODE %d depth %d\n", *leaf_count, depth);
    }else if(tree->type == LEAF)
    {
        ++(*leaf_count);
        printf("LEAF %d depth %d\n", *leaf_count, depth);
        printVec3WithText("aabb.min", tree->aabb.min);
        printVec3WithText("aabb.max", tree->aabb.max);
        if(tree->obj_type != SPHERE)
        {
            printf("WHAT THE HELL\n");
        }
    }else
    {
        printf("WTF NEITHER?\n");
    }
    if(tree->left){printBVH(tree->left, leaf_count, depth+1);}
    if(tree->right){printBVH(tree->right, leaf_count, depth+1);}
}
