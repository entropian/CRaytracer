#pragma once
#include <cmath>
#include "../util/constants.h"
#include "../util/vec.h"
#include "../util/intvector.h"
#include "../aabb.h"
#include "../shapes/shapes.h"
#include "../shapes/instanced.h"

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
    free(ug->cells);
}

void getSceneObjAABBs(AABB* aabb_array, int* num_aabb, const Object_t objects[], const int num_obj)
{   
    int aabb_index = 0;
    for(int i = 0; i < num_obj; ++i)
    {
        if(isGridObjType(objects[i]))
        {
            getObjectAABB(&(aabb_array[aabb_index++]), objects[i]);
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
UniformGrid* UniformGrid_create(Object_t* objects, int* num_obj, const int num_non_grid_obj, const int multiplier)    
{
    // Calc aabb for every obj and store in a array
    // use array to calc grid aabb
    // later use array to place obj in cells
    UniformGrid* ug = (UniformGrid*)malloc(sizeof(UniformGrid));
    AABB* aabb_array = (AABB*)malloc(sizeof(AABB) * *num_obj);
    int num_aabb;
    getSceneObjAABBs(aabb_array, &num_aabb, objects, *num_obj);
    ug->aabb = calcUniformGridAABB(aabb_array, num_aabb);

    // Grid dimensions
    ug->wx = ug->aabb.max[0] - ug->aabb.min[0];
    ug->wy = ug->aabb.max[1] - ug->aabb.min[1];
    ug->wz = ug->aabb.max[2] - ug->aabb.min[2];
    printVec3WithText("Uniform grid max", ug->aabb.max);
    printVec3WithText("Uniform grid min", ug->aabb.min);

    // Number of cells per dimension
    float s = powf(ug->wx * ug->wy * ug->wz / *num_obj, 0.33333f);    
    ug->nx = (int)(ug->wx * multiplier / s + 1);
    ug->ny = (int)(ug->wy * multiplier / s + 1);
    ug->nz = (int)(ug->wz * multiplier / s + 1);
    ug->num_cells = ug->nx * ug->ny * ug->nz;
    // TODO find out how many cells are empty and the average number of indices per cell
    printf("num_cells %d\n", ug->num_cells);

    int nx = ug->nx, ny = ug->ny, nz = ug->nz;
    float wx = ug->wx, wy = ug->wy, wz = ug->wz;
    // grid cell dimension
    float unit_x = wx / ug->nx;
    float unit_y = wy / ug->ny;
    float unit_z = wz / ug->nz;

    ug->cells = (IntVector*)malloc(ug->num_cells * sizeof(IntVector));
    for(int i = 0; i < ug->num_cells; ++i)
    {
        ug->cells[i] = IntVector_create();
    }

    for(int i = 0; i < num_aabb; ++i)
    {
        int min_ix = (int)clamp(((aabb_array[i].min[0] - ug->aabb.min[0]) / wx) * nx, 0.0f, (float)nx - 1);
        int min_iy = (int)clamp(((aabb_array[i].min[1] - ug->aabb.min[1]) / wy) * ny, 0.0f, (float)ny - 1);
        int min_iz = (int)clamp(((aabb_array[i].min[2] - ug->aabb.min[2]) / wz) * nz, 0.0f, (float)nz - 1);

        int max_ix = (int)clamp(((aabb_array[i].max[0] - ug->aabb.min[0]) / wx) * nx, 0.0f, (float)nx - 1);
        int max_iy = (int)clamp(((aabb_array[i].max[1] - ug->aabb.min[1]) / wy) * ny, 0.0f, (float)ny - 1);
        int max_iz = (int)clamp(((aabb_array[i].max[2] - ug->aabb.min[2]) / wz) * nz, 0.0f, (float)nz - 1);
        
        int num_good = 0;
        int num_bad = 0;
        for(int j = min_iy; j <= max_iy; ++j)
        {
            for(int k = min_iz; k <= max_iz; ++k)
            {
                for(int p = min_ix; p <= max_ix; ++p)
                {
                    Object_t obj = objects[i + num_non_grid_obj];
                    if((obj.type == TRIANGLE || obj.type == FLAT_TRIANGLE || obj.type == SMOOTH_TRIANGLE))
                    {
                        AABB aabb;
                        vec3 min_disp = {p * unit_x, j * unit_y, k * unit_z};
                        vec3_add(aabb.min, ug->aabb.min, min_disp);
                        vec3 unit_span = {unit_x, unit_y, unit_z};
                        vec3_add(aabb.max, aabb.min, unit_span);
                        if(triangleAABBIntersect(obj, &aabb))
                        {
                            int cell_index = j*nx*nz + k*nx + p;
                            IntVector_push(&(ug->cells[j*nx*nz + k*nx + p]), i + num_non_grid_obj);
                        }
                    }else
                    {
                        int cell_index = j*nx*nz + k*nx + p;
                        IntVector_push(&(ug->cells[j*nx*nz + k*nx + p]), i + num_non_grid_obj);
                    }
                }
            }
        }
        //printf("num_good %d num_bad %d\n", num_good, num_bad);
    }
    int empty_count = 0;
    unsigned long long total_count = 0;
    for(int i = 0; i < ug->num_cells; i++)
    {
        total_count += ug->cells[i].size;
        if(ug->cells[i].size == 0)
        {
            empty_count++;
        }
    }
    float avg = (float)total_count / (float)(ug->num_cells - empty_count);
    printf("Number of empty cells %d\n", empty_count);
    printf("Fraction of empty cells %f\n", (float)empty_count / (float)ug->num_cells);
    printf("Average number of primitives in non-empty cells %f\n", avg);
    
    free(aabb_array);
    printf("wx = %f, wy = %f, wz = %f\n", ug->wx, ug->wy, ug->wz);
    printf("nx = %d, ny = %d, nz = %d\n", nx, ny, nz);
    return ug;
}
