#pragma once

#include "shapes/shapes.h"
#include "shapes/instanced.h"
#include "lights.h"
#include "accelerator/uniformgrid.h"
#include "accelerator/bvh.h"
#include "accelerator/bvh4.h"

float gridIntersectTest(ShadeRec* sr, const SceneObjects* so, const Ray ray)
{
    // first find where the ray intersects with the aabb or if its origin is inside the aabb
    const UniformGrid* rg = (UniformGrid*)(so->accel_ptr);
    bool hits_grid = false;
    int ix, iy, iz;
    float t0;
    if(isInsideAABB(&(rg->aabb), ray.origin))
    {
        ix = (int)clamp((ray.origin[0] - rg->aabb.min[0]) / rg->wx * rg->nx, 0.0f, (float)rg->nx - 1);
        iy = (int)clamp((ray.origin[1] - rg->aabb.min[1]) / rg->wy * rg->ny, 0.0f, (float)rg->ny - 1);
        iz = (int)clamp((ray.origin[2] - rg->aabb.min[2]) / rg->wz * rg->nz, 0.0f, (float)rg->nz - 1);
        hits_grid = true;
    }else
    {
        // Calculate where the ray hits the aabb
        t0 = rayIntersectAABB(&(rg->aabb), ray);
        if(t0 < TMAX)
        {
            vec3 point;
            getPointOnRay(point, ray, t0);
            ix = (int)clamp((point[0] - rg->aabb.min[0]) / rg->wx * rg->nx, 0.0f, (float)rg->nx - 1);
            iy = (int)clamp((point[1] - rg->aabb.min[1]) / rg->wy * rg->ny, 0.0f, (float)rg->ny - 1);
            iz = (int)clamp((point[2] - rg->aabb.min[2]) / rg->wz * rg->nz, 0.0f, (float)rg->nz - 1);
            hits_grid = true;
        }
    }

    float t = TMAX;
    if(hits_grid)
    {
        float x_comp = ray.direction[0];
        float y_comp = ray.direction[1];
        float z_comp = ray.direction[2];

        float x_min = x_comp > 0 ? rg->aabb.min[0] : rg->aabb.max[0];
        float y_min = y_comp > 0 ? rg->aabb.min[1] : rg->aabb.max[1];
        float z_min = z_comp > 0 ? rg->aabb.min[2] : rg->aabb.max[2];
        float x_max = x_comp > 0 ? rg->aabb.max[0] : rg->aabb.min[0];
        float y_max = y_comp > 0 ? rg->aabb.max[1] : rg->aabb.min[1];
        float z_max = z_comp > 0 ? rg->aabb.max[2] : rg->aabb.min[2];

        float tx_min = (x_min - ray.origin[0]) / x_comp;
        float ty_min = (y_min - ray.origin[1]) / y_comp;
        float tz_min = (z_min - ray.origin[2]) / z_comp;

        float tx_max = (x_max - ray.origin[0]) / x_comp;
        float ty_max = (y_max - ray.origin[1]) / y_comp;
        float tz_max = (z_max - ray.origin[2]) / z_comp;

        float dtx = (tx_max - tx_min) / rg->nx;
        float dty = (ty_max - ty_min) / rg->ny;
        float dtz = (tz_max - tz_min) / rg->nz;            
            
        int ix_step = x_comp > 0 ? 1 : -1;
        int iy_step = y_comp > 0 ? 1 : -1;
        int iz_step = z_comp > 0 ? 1 : -1;

        int ix_stop = x_comp > 0 ? rg->nx : -1;
        int iy_stop = y_comp > 0 ? rg->ny : -1;
        int iz_stop = z_comp > 0 ? rg->nz : -1;            

        int tmp_ix = x_comp > 0 ? ix : rg->nx - ix - 1;
        int tmp_iy = y_comp > 0 ? iy : rg->ny - iy - 1;
        int tmp_iz = z_comp > 0 ? iz : rg->nz - iz - 1;

        float tx_next = tx_min + dtx * (tmp_ix + 1);
        float ty_next = ty_min + dty * (tmp_iy + 1);
        float tz_next = tz_min + dtz * (tmp_iz + 1);

        int cell_index;
        while(ix != ix_stop && iy != iy_stop && iz != iz_stop)
        {
            cell_index = iy*rg->nz*rg->nx + iz*rg->nx + ix;                
            float obj_t = TMAX;
            ShadeRec obj_sr;
            if(rg->cells[cell_index].size > 0)
            {
                float tmp_t = TMAX;
                ShadeRec tmp_sr;
                for(int i = 0; i < rg->cells[cell_index].size; i++)
                {
                    int index = IntVector_get(&(rg->cells[cell_index]), i);
                    tmp_t = rayIntersectObject(&tmp_sr, so->objects[index], ray);
                    if(tmp_t < obj_t)
                    {
                        obj_t = tmp_t;
                        obj_sr = tmp_sr;
                    }
                }
            }
            if(min(obj_t, min(tx_next, min(ty_next, tz_next))) == obj_t)                
            {
                if(obj_t < TMAX)
                {
                    *sr = obj_sr;
                    t = obj_t;
                    break;
                }
            }
            if(tx_next < ty_next && tx_next < tz_next)
            {
                tx_next += dtx;
                ix += ix_step;
            }else if(ty_next < tz_next)
            {
                ty_next += dty;
                iy += iy_step;
            }else
            {
                tz_next += dtz;
                iz += iz_step;
            }
        }
    }
    // Test intersection with non-bounded objects
    float tmp_t = TMAX,  min_t = TMAX;
    ShadeRec tmp_sr, min_sr;
    for(int i = 0; i < so->num_non_grid_obj; i++)
    {
        tmp_t = rayIntersectObject(&tmp_sr, so->objects[i], ray);
        if(tmp_t < min_t)
        {
            min_t = tmp_t;
            min_sr = tmp_sr;
        }
    }
    if(min_t < t)
    {
        *sr = min_sr;
        t = min_t;
    }
    return t;        
}

float gridShadowIntersectTest(const SceneObjects* so, const Ray shadow_ray)
{
    // variables to make sense of: tx_max and tx_min,
    // ix_step = 1, ix_stop = nx;
    // how to find tx_max and tx_min

    // first find where the ray intersects with the aabb or if its origin is inside the aabb
    const UniformGrid* rg = (UniformGrid*)(so->accel_ptr);
    bool hits_grid = false;
    int ix, iy, iz;
    float t0;
    if(isInsideAABB(&(rg->aabb), shadow_ray.origin))
    {
        ix = (int)clamp((shadow_ray.origin[0] - rg->aabb.min[0]) / rg->wx * rg->nx, 0.0f, (float)rg->nx - 1);
        iy = (int)clamp((shadow_ray.origin[1] - rg->aabb.min[1]) / rg->wy * rg->ny, 0.0f, (float)rg->ny - 1);
        iz = (int)clamp((shadow_ray.origin[2] - rg->aabb.min[2]) / rg->wz * rg->nz, 0.0f, (float)rg->nz - 1);
        hits_grid = true;
    }else
    {
        // Calculate where the ray hits the aabb
        t0 = rayIntersectAABB(&(rg->aabb), shadow_ray);
        if(t0 < TMAX)
        {
            vec3 point;
            getPointOnRay(point, shadow_ray, t0);
            ix = (int)clamp((point[0] - rg->aabb.min[0]) / rg->wx * rg->nx, 0.0f, (float)rg->nx - 1);
            iy = (int)clamp((point[1] - rg->aabb.min[1]) / rg->wy * rg->ny, 0.0f, (float)rg->ny - 1);
            iz = (int)clamp((point[2] - rg->aabb.min[2]) / rg->wz * rg->nz, 0.0f, (float)rg->nz - 1);
            hits_grid = true;
        }
    }

    if(hits_grid)
    {
        float x_comp = shadow_ray.direction[0];
        float y_comp = shadow_ray.direction[1];
        float z_comp = shadow_ray.direction[2];

        float x_min = x_comp > 0 ? rg->aabb.min[0] : rg->aabb.max[0];
        float y_min = y_comp > 0 ? rg->aabb.min[1] : rg->aabb.max[1];
        float z_min = z_comp > 0 ? rg->aabb.min[2] : rg->aabb.max[2];
        float x_max = x_comp > 0 ? rg->aabb.max[0] : rg->aabb.min[0];
        float y_max = y_comp > 0 ? rg->aabb.max[1] : rg->aabb.min[1];
        float z_max = z_comp > 0 ? rg->aabb.max[2] : rg->aabb.min[2];

        float tx_min = (x_min - shadow_ray.origin[0]) / x_comp;
        float ty_min = (y_min - shadow_ray.origin[1]) / y_comp;
        float tz_min = (z_min - shadow_ray.origin[2]) / z_comp;

        float tx_max = (x_max - shadow_ray.origin[0]) / x_comp;
        float ty_max = (y_max - shadow_ray.origin[1]) / y_comp;
        float tz_max = (z_max - shadow_ray.origin[2]) / z_comp;

        float dtx = (tx_max - tx_min) / rg->nx;
        float dty = (ty_max - ty_min) / rg->ny;
        float dtz = (tz_max - tz_min) / rg->nz;            
            
        int ix_step = x_comp > 0 ? 1 : -1;
        int iy_step = y_comp > 0 ? 1 : -1;
        int iz_step = z_comp > 0 ? 1 : -1;

        int ix_stop = x_comp > 0 ? rg->nx : -1;
        int iy_stop = y_comp > 0 ? rg->ny : -1;
        int iz_stop = z_comp > 0 ? rg->nz : -1;            

        int tmp_ix = x_comp > 0 ? ix : rg->nx - ix - 1;
        int tmp_iy = y_comp > 0 ? iy : rg->ny - iy - 1;
        int tmp_iz = z_comp > 0 ? iz : rg->nz - iz - 1;

        float tx_next = tx_min + dtx * (tmp_ix + 1);
        float ty_next = ty_min + dty * (tmp_iy + 1);
        float tz_next = tz_min + dtz * (tmp_iz + 1);                        

        int cell_index;
        while(ix != ix_stop && iy != iy_stop && iz != iz_stop)                
        {
            cell_index = iy*rg->nz*rg->nx + iz*rg->nx + ix;                
            float obj_t = TMAX;
            if(rg->cells[cell_index].size > 0)
            {
                float t = TMAX;
                for(int i = 0; i < rg->cells[cell_index].size; i++)
                {
                    int index = IntVector_get(&(rg->cells[cell_index]), i);
                    t = shadowRayIntersectObject(so->objects[index], shadow_ray);
                    if(t < TMAX)
                    {
                        return t;
                    }
                }
            }
            if(tx_next < ty_next && tx_next < tz_next)
            {
                tx_next += dtx;
                ix += ix_step;
            }else if(ty_next < tz_next)
            {
                ty_next += dty;
                iy += iy_step;
            }else
            {
                tz_next += dtz;
                iz += iz_step;
            }
        }
    }
    // Test intersection with non-bounded objects
    float t = TMAX;
    for(int i = 0; i < so->num_non_grid_obj; i++)
    {
        t = shadowRayIntersectObject(so->objects[i], shadow_ray);
        if(t < TMAX)
        {
            return t;
        }
    }
    return t;
}

float BVHIntersectTest(ShadeRec* sr, const SceneObjects* so, const BVHNode* tree, const Ray ray)
{
    if(tree->type == LEAF)
    {
        return rayIntersectObject(sr, tree->obj, ray);
    }else
    {
        float t1 = TMAX, t2 = TMAX;
        ShadeRec sr1, sr2;
        if(rayIntersectAABB(&(tree->left->aabb), ray) < TMAX)
        {
            t1 = BVHIntersectTest(&sr1, so, tree->left, ray);
        }
        if(rayIntersectAABB(&(tree->right->aabb), ray) < TMAX)
        {
            t2 = BVHIntersectTest(&sr2, so, tree->right, ray);
        }
        if(t1 < t2)
        {
            *sr = sr1;
            return t1;
        }else
        {
			if (t2 < TMAX)
			{
				*sr = sr2;
			}
            return t2;
        }
        // TODO
        //return t1 < t2 ? t1 : t2;
    }
}

float BVHShadowIntersectTest(const BVHNode* tree, const Ray ray)
{
    if(tree->type == LEAF)
    {
        return shadowRayIntersectObject(tree->obj, ray);
    }else
    {
        float t1 = TMAX, t2 = TMAX;
        if(rayIntersectAABB(&(tree->left->aabb), ray) < TMAX)
        {
            t1 = BVHShadowIntersectTest(tree->left, ray);
        }
        if(rayIntersectAABB(&(tree->right->aabb), ray) < TMAX)
        {
            t2 = BVHShadowIntersectTest(tree->right, ray);
        }
        return t1 < t2 ? t1 : t2;
    }
}

float intersectTest(ShadeRec* sr, const SceneObjects* so, const Ray ray)
{
    if(so->accel == GRID)
    {
        return gridIntersectTest(sr, so, ray);
    }else if(so->accel == BVH)
    {
        BVHNode* tree = (BVHNode*)(so->accel_ptr);
        float t = TMAX;
        if(rayIntersectAABB(&(tree->aabb), ray) < TMAX)
        {
            t = BVHIntersectTest(sr, so, tree, ray);
        }
        float tmp_t = TMAX,  min_t = TMAX;
        ShadeRec tmp_sr, min_sr;
        for(int i = 0; i < so->num_non_grid_obj; i++)
        {
            tmp_t = rayIntersectObject(&tmp_sr, so->objects[i], ray);
            if(tmp_t < min_t)
            {
                min_t = tmp_t;
                min_sr = tmp_sr;
            }
        }
        if(min_t < t)
        {
            *sr = min_sr;
            return min_t;
        }
        return t;
    }else
    {
        float tmp_t = TMAX,  min_t = TMAX;
        ShadeRec tmp_sr, min_sr;
        for(int i = 0; i < so->num_obj; i++)
        {
            tmp_t = rayIntersectObject(&tmp_sr, so->objects[i], ray);
            if(tmp_t < min_t)
            {
                min_t = tmp_t;
                min_sr = tmp_sr;
            }
        }
        if(min_t < TMAX)
        {
            *sr = min_sr;
        }
        return min_t;
    }
}



float shadowIntersectTest(const SceneObjects *so, const Ray shadow_ray)
{
    if(so->accel == GRID)
    {
        return gridShadowIntersectTest(so, shadow_ray);
    }else if(so->accel == BVH)
    {
        BVHNode* tree = (BVHNode*)(so->accel_ptr);
        float t = TMAX;
        if(rayIntersectAABB(&(tree->aabb), shadow_ray) < TMAX)
        {
            t = BVHShadowIntersectTest(tree, shadow_ray);
        }
        if(t == TMAX)
        { 
            for(int i = 0; i < so->num_non_grid_obj; i++)
            {
                t = shadowRayIntersectObject(so->objects[i], shadow_ray);
                if(t < TMAX)
                {
                    return t;
                }
            }
        }
        return t;        
    }else
    {
        float t = TMAX;
        for(int i = 0; i < so->num_obj; i++)
        {
            t = shadowRayIntersectObject(so->objects[i], shadow_ray);
            if(t < TMAX)
            {
                return t;
            }
        }
        return t;        
    }
}
