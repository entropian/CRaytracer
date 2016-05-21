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

typedef struct RegularGrid
{
    int num_cells;
    int nx, ny, nz;    // number of cells along each dimension
    float wx, wy, wz;    // length of the bounding box along each dimension
    IntVector* cells;
    AABB aabb;
} RegularGrid;

// TODO: add return value
void getObjectAABB(AABB* aabb, const void* obj_ptr, const ObjectType obj_type)
{
    switch(obj_type)
    {
    case SPHERE:
    {
        Sphere* sphere = (Sphere*)obj_ptr;
        aabb->min[0] = sphere->center[0] - sphere->radius;
        aabb->min[1] = sphere->center[1] - sphere->radius;
        aabb->min[2] = sphere->center[2] - sphere->radius;

        aabb->max[0] = sphere->center[0] + sphere->radius;
        aabb->max[1] = sphere->center[1] + sphere->radius;
        aabb->max[2] = sphere->center[2] + sphere->radius;
    } break;
    case RECTANGLE:
    {
        // NOTE: shit code
        Rectangle* rect = (Rectangle*)obj_ptr;
        vec3 pts[4];
        vec3 pt1, pt2, pt3, pt4;
        vec3_copy(pts[0], rect->point);
        vec3_add(pts[1], pts[0], rect->width);
        vec3_add(pts[2], pts[0], rect->height);
        vec3 tmp;
        vec3_add(tmp, rect->width, rect->height);
        vec3_add(pts[3], pts[0], tmp);
            
        vec3_copy(aabb->min, pts[0]);
        vec3_copy(aabb->max, pts[0]);

        for(int i = 1; i < 4; i++)
        {
            for(int j = 0; j < 3; j++)
            {
                if(pts[i][j] < aabb->min[j]){aabb->min[j] = pts[i][j] - K_EPSILON;}
                if(pts[i][j] > aabb->max[j]){aabb->max[j] = pts[i][j] + K_EPSILON;}
            }
        }
    } break;
    case AABOX:
    {
        AABox* aabox = (AABox*)obj_ptr;
        vec3_copy(aabb->min, aabox->min);
        vec3_copy(aabb->max, aabox->max);            
    } break;
    case TRIANGLE:
    {
        Triangle* triangle = (Triangle*)obj_ptr;
        vec3_copy(aabb->min, triangle->v0);
        vec3_copy(aabb->max, triangle->v0);

        for(int i = 0; i < 3; i++)
        {
            if(triangle->v1[i] < aabb->min[i]){aabb->min[i] = triangle->v1[i];}
            if(triangle->v1[i] > aabb->max[i]){aabb->max[i] = triangle->v1[i];}
            if(triangle->v2[i] < aabb->min[i]){aabb->min[i] = triangle->v2[i];}
            if(triangle->v2[i] > aabb->max[i]){aabb->max[i] = triangle->v2[i];}
        }
        for(int i = 0; i < 3; i++)
        {
            aabb->min[i] -= K_FLAT_AABB;
            aabb->max[i] += K_FLAT_AABB;
        }
    } break;
    case OPENCYLINDER:
    {
        OpenCylinder* oc = (OpenCylinder*)obj_ptr;
        aabb->min[0] = -(oc->radius);
        aabb->min[1] = -(oc->half_height);
        aabb->min[2] = -(oc->radius);
        aabb->max[0] = oc->radius;
        aabb->max[1] = oc->half_height;
        aabb->max[2] = oc->radius;
    } break;
    case DISK:
    {
        // TODO: improve this crap
        Disk* disk = (Disk*)obj_ptr;
        aabb->min[0] = disk->center[0] - disk->radius;
        aabb->min[1] = disk->center[1] - disk->radius;
        aabb->min[2] = disk->center[2] - disk->radius;

        aabb->max[0] = disk->center[0] + disk->radius;
        aabb->max[1] = disk->center[1] + disk->radius;
        aabb->max[2] = disk->center[2] + disk->radius;        
    } break;
    case TORUS:
    {
        Torus* torus = (Torus*)obj_ptr;
        *aabb = torus->aabb;
    } break;
    case INSTANCED:
    {
        InstancedShape* is = (InstancedShape*)obj_ptr;
        getObjectAABB(aabb, is->obj_ptr, is->obj_type);
    } break;
    case COMPOUND:
    {
        CompoundObject* co = (CompoundObject*)obj_ptr;
        *aabb = co->aabb;
    } break;
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

//AABB calcRegularGridAABB(const SceneObjects* so)
AABB calcRegularGridAABB(const AABB aabb_array[], const int num_aabb)
{
    AABB r;
    vec3_assign(r.min, FLT_MAX, FLT_MAX, FLT_MAX);
    vec3_assign(r.max, -HUGEVALUE, -HUGEVALUE, -HUGEVALUE);
    for(int i = 0; i < num_aabb; ++i)
    {
        if(aabb_array[i].min[0] < r.min[0]){r.min[0] = aabb_array[i].min[0];}
        if(aabb_array[i].min[1] < r.min[1]){r.min[1] = aabb_array[i].min[1];}
        if(aabb_array[i].min[2] < r.min[2]){r.min[2] = aabb_array[i].min[2];}

        if(aabb_array[i].max[0] > r.max[0]){r.max[0] = aabb_array[i].max[0];}
        if(aabb_array[i].max[1] > r.max[1]){r.max[1] = aabb_array[i].max[1];}
        if(aabb_array[i].max[2] > r.max[2]){r.max[2] = aabb_array[i].max[2];}                
    }
    return r;
}

// Assume non-grid objects are in the beginning of the array
//void RegularGrid_create(RegularGrid* rg, const void* obj_ptrs[], const ObjectType obj_types[],
//                        const int num_obj, const int num_non_grid_obj, const int multiplier)
void RegularGrid_create(RegularGrid* rg, void* obj_ptrs[], ObjectType obj_types[],
                        int* num_obj, const int num_non_grid_obj, const int multiplier)    
{
    // Calc aabb for every obj and store in a array
    // use array to calc grid aabb
    // later use array to place obj in cells
    AABB* aabb_array = (AABB*)malloc(sizeof(AABB) * *num_obj);
    int num_aabb;
    getSceneObjAABBs(aabb_array, &num_aabb, (const void**)obj_ptrs, obj_types, *num_obj);
    rg->aabb = calcRegularGridAABB(aabb_array, num_aabb);
    
    rg->wx = rg->aabb.max[0] - rg->aabb.min[0];
    rg->wy = rg->aabb.max[1] - rg->aabb.min[1];
    rg->wz = rg->aabb.max[2] - rg->aabb.min[2];
    printVec3WithText("aabb.max", rg->aabb.max);
    printVec3WithText("aabb.min", rg->aabb.min);    

    //float s = pow(rg->wx * rg->wy * rg->wz / *num_obj, 0.33333);
    float s = pow(rg->wx * rg->wy * rg->wz / *num_obj, 0.33333);
    rg->nx = rg->wx * multiplier / s + 1;
    rg->ny = rg->wy * multiplier / s + 1;
    rg->nz = rg->wz * multiplier / s + 1;
    rg->num_cells = rg->nx * rg->ny * rg->nz;

    rg->cells = (IntVector*)malloc(rg->num_cells * sizeof(IntVector));
    for(int i = 0; i < rg->num_cells; ++i)
    {
        rg->cells[i] = IntVector_create();
    }

    for(int i = 0; i < num_aabb; ++i)
    {
        int min_ix = clamp(((aabb_array[i].min[0] - rg->aabb.min[0]) / rg->wx) * rg->nx, 0, rg->nx - 1);
        int min_iy = clamp(((aabb_array[i].min[1] - rg->aabb.min[1]) / rg->wy) * rg->ny, 0, rg->ny - 1);
        int min_iz = clamp(((aabb_array[i].min[2] - rg->aabb.min[2]) / rg->wz) * rg->nz, 0, rg->nz - 1);

        int max_ix = clamp(((aabb_array[i].max[0] - rg->aabb.min[0]) / rg->wx) * rg->nx, 0, rg->nx - 1);
        int max_iy = clamp(((aabb_array[i].max[1] - rg->aabb.min[1]) / rg->wy) * rg->ny, 0, rg->ny - 1);
        int max_iz = clamp(((aabb_array[i].max[2] - rg->aabb.min[2]) / rg->wz) * rg->nz, 0, rg->nz - 1);

        for(int j = min_iy; j <= max_iy; ++j)
        {
            for(int k = min_iz; k <= max_iz; ++k)
            {
                for(int p = min_ix; p <= max_ix; ++p)
                {
                    IntVector_push(&(rg->cells[j*rg->nx*rg->nz + k*rg->nx + p]), i + num_non_grid_obj);
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
    printf("wx = %f, wy = %f, wz = %f\n", rg->wx, rg->wy, rg->wz);
    printf("nx = %d, ny = %d, nz = %d\n", rg->nx, rg->ny, rg->nz);
}
