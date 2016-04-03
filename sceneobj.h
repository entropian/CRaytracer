#pragma once

#include "shapes.h"
#include "lights.h"

typedef struct SceneObjects
{
    int num_obj;
    ObjectType obj_types[MAX_OBJECTS];
    void* obj_ptrs[MAX_OBJECTS];
} SceneObjects;

typedef struct SceneLights
{
    int num_lights;
    bool shadow[MAX_LIGHTS];    
    LightType light_types[MAX_LIGHTS];
    void* light_ptrs[MAX_LIGHTS];
} SceneLights;

void freeSceneObjects(SceneObjects* so)
{
    for(int i = 0; i < so->num_obj; i++)
    {
        if(so->obj_ptrs[i] != NULL)
        {
            free(so->obj_ptrs[i]);
            so->obj_ptrs[i] = NULL;
        }
    }
}

void freeSceneLights(SceneLights* sl)
{
    for(int i = 0; i < sl->num_lights; i++)
    {
        if(sl->light_ptrs[i] != NULL)
        {
            if(sl->light_types[i] == AREALIGHT)
            {
                AreaLight* area_light_ptr = (AreaLight*)(sl->light_ptrs[i]);
                if(area_light_ptr->samples2D != NULL)
                {
                    freeSamples2D(area_light_ptr->samples2D);
                }
                if(area_light_ptr->samples3D != NULL)
                {
                    freeSamples3D(area_light_ptr->samples3D);
                }
            }                
            free(sl->light_ptrs[i]);
            sl->light_ptrs[i] = NULL;
        }
    }
}

float intersectTest(ShadeRec *sr, const SceneObjects *scene_objects, const Ray ray)
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
        case RECTANGLE:
            tmp_t = rayIntersectRect(&tmp_sr, (Rectangle*)scene_objects->obj_ptrs[i], ray);
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
float shadowIntersectTest(const SceneObjects *scene_objects, const Ray shadow_ray)
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
        case RECTANGLE:
            tmp_t = shadowRayIntersectRect((Rectangle*)scene_objects->obj_ptrs[i], shadow_ray);
            break;            
        }
        if(tmp_t < min_t)
        {
            min_t = tmp_t;
        }
    }
    return min_t;
}



