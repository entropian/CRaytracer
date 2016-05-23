#pragma once
#include <cmath>
#include "../util/constants.h"
#include "../util/vec.h"
#include "../util/intvector.h"
#include "../aabb.h"
#include "../shapes/shapes.h"
#include "../shapes/instanced.h"

typedef struct GridCell
{
    int length;
    IntVector indicies;
} GridCell;

typedef struct UniformGrid
{
    int num_cells;
    int nx, ny, nz;    // number of cells along each dimension
    float wx, wy, wz;    // length of the bounding box along each dimension
    IntVector* cells;
    AABB aabb;
} UniformGrid;

void UniformGrid_destroy(UniformGrid* ug)
{
    for(int i = 0; i < ug->num_cells; i++)
    {
        IntVector_destroy(&(ug->cells[i]));
    }
}

void getSceneObjAABBs(AABB* aabb_array, int* num_aabb,
                      const void* obj_ptrs[], const ObjectType obj_types[], const int num_obj)
{   
    int aabb_index = 0;
    for(int i = 0; i < num_obj; ++i)
    {
        if(isGridObjType(obj_types[i]))
        {
            getObjectAABB(&(aabb_array[aabb_index++]), obj_ptrs[i], obj_types[i]);
        }
    }
    *num_aabb = aabb_index;
}

AABB calcUniformGridAABB(const AABB aabb_array[], const int num_aabb)
{
    AABB r;
    vec3_assign(r.min, FLT_MAX, FLT_MAX, FLT_MAX);
    vec3_assign(r.max, -HUGEVALUE, -HUGEVALUE, -HUGEVALUE);
    for(int i = 0; i < num_aabb; ++i)
    {
        addToAABB(&r, &(aabb_array[i]));
    }
    return r;
}

// Assume non-grid objects are in the beginning of the array
//void UniformGrid_create(UniformGrid* rg, const void* obj_ptrs[], const ObjectType obj_types[],
//                        const int num_obj, const int num_non_grid_obj, const int multiplier)
UniformGrid* UniformGrid_create(void* obj_ptrs[], ObjectType obj_types[],
                        int* num_obj, const int num_non_grid_obj, const int multiplier)    
{
    // Calc aabb for every obj and store in a array
    // use array to calc grid aabb
    // later use array to place obj in cells
    UniformGrid* ug = (UniformGrid*)malloc(sizeof(UniformGrid));
    AABB* aabb_array = (AABB*)malloc(sizeof(AABB) * *num_obj);
    int num_aabb;
    getSceneObjAABBs(aabb_array, &num_aabb, (const void**)obj_ptrs, obj_types, *num_obj);
    ug->aabb = calcUniformGridAABB(aabb_array, num_aabb);
    
    ug->wx = ug->aabb.max[0] - ug->aabb.min[0];
    ug->wy = ug->aabb.max[1] - ug->aabb.min[1];
    ug->wz = ug->aabb.max[2] - ug->aabb.min[2];
    printVec3WithText("aabb.max", ug->aabb.max);
    printVec3WithText("aabb.min", ug->aabb.min);    

    //float s = pow(ug->wx * ug->wy * ug->wz / *num_obj, 0.33333);
    float s = pow(ug->wx * ug->wy * ug->wz / *num_obj, 0.33333);
    ug->nx = ug->wx * multiplier / s + 1;
    ug->ny = ug->wy * multiplier / s + 1;
    ug->nz = ug->wz * multiplier / s + 1;
    ug->num_cells = ug->nx * ug->ny * ug->nz;

    ug->cells = (IntVector*)malloc(ug->num_cells * sizeof(IntVector));
    for(int i = 0; i < ug->num_cells; ++i)
    {
        ug->cells[i] = IntVector_create();
    }

    for(int i = 0; i < num_aabb; ++i)
    {
        int min_ix = clamp(((aabb_array[i].min[0] - ug->aabb.min[0]) / ug->wx) * ug->nx, 0, ug->nx - 1);
        int min_iy = clamp(((aabb_array[i].min[1] - ug->aabb.min[1]) / ug->wy) * ug->ny, 0, ug->ny - 1);
        int min_iz = clamp(((aabb_array[i].min[2] - ug->aabb.min[2]) / ug->wz) * ug->nz, 0, ug->nz - 1);

        int max_ix = clamp(((aabb_array[i].max[0] - ug->aabb.min[0]) / ug->wx) * ug->nx, 0, ug->nx - 1);
        int max_iy = clamp(((aabb_array[i].max[1] - ug->aabb.min[1]) / ug->wy) * ug->ny, 0, ug->ny - 1);
        int max_iz = clamp(((aabb_array[i].max[2] - ug->aabb.min[2]) / ug->wz) * ug->nz, 0, ug->nz - 1);

        for(int j = min_iy; j <= max_iy; ++j)
        {
            for(int k = min_iz; k <= max_iz; ++k)
            {
                for(int p = min_ix; p <= max_ix; ++p)
                {
                    IntVector_push(&(ug->cells[j*ug->nx*ug->nz + k*ug->nx + p]), i + num_non_grid_obj);
                }
            }
        }
        /*
        AABox* aabox = (AABox*)malloc(sizeof(AABox));
        aabox->shadow = true;
        vec3_copy(aabox->min, aabb_array[i].min);
        vec3_copy(aabox->max, aabb_array[i].max);
        float tmp = 0.5f;
        vec3 color;
        vec3_assign(color, tmp, tmp, tmp);
        initDefaultPhongMat(&(aabox->mat), color);
        obj_ptrs[*num_obj] = aabox;
        obj_types[*num_obj] = AABOX;
        (*num_obj)++;
        */
    }
    free(aabb_array);
    printf("wx = %f, wy = %f, wz = %f\n", ug->wx, ug->wy, ug->wz);
    printf("nx = %d, ny = %d, nz = %d\n", ug->nx, ug->ny, ug->nz);
    return ug;
}
