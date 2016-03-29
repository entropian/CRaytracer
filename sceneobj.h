#pragma once

#include "shapes.h"

typedef struct SceneObjects
{
    int num_obj;
    ObjectType obj_types[MAX_SPHERES];
    void* obj_ptrs[MAX_SPHERES];
} SceneObjects2;

float intersectTest(ShadeRec *sr, SceneObjects *scene_objects, const Ray ray)
{
    float tmp_t = TMAX,  min_t = TMAX;
    ShadeRec tmp_sr, min_sr;
    for(int i = 0; i < scene_objects->num_obj; i++)
    {
        switch(scene_objects->obj_types[i])
        {
        case SPHERE:
            tmp_t = rayIntersectSphere(&tmp_sr, (Sphere*)scene_objects->obj_ptrs[i], ray);            
            break;
        case PLANE:
            tmp_t = rayIntersectPlane(&tmp_sr, (Plane*)scene_objects->obj_ptrs[i], ray);                        
            break;
        }
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

// NOTE: return after the first intersection for certain situations?
float shadowIntersectTest(SceneObjects *scene_objects, const Ray shadow_ray)
{
    float tmp_t, min_t = TMAX;

    for(int i = 0; i < scene_objects->num_obj; i++)
    {
        switch(scene_objects->obj_types[i])
        {
        case SPHERE:
            tmp_t = shadowRayIntersectSphere((Sphere*)scene_objects->obj_ptrs[i], shadow_ray);
            break;
        case PLANE:
            tmp_t = shadowRayIntersectPlane((Plane*)scene_objects->obj_ptrs[i], shadow_ray);
            break;
        }
        if(tmp_t < min_t)
        {
            min_t = tmp_t;
        }
    }
    return min_t;
}


