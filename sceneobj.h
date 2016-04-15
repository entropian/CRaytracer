#pragma once

#include "shapes/shapes.h"
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
    EnvLight* env_light_ptr;
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
                    area_light_ptr->samples2D = NULL;
                }
                if(area_light_ptr->samples3D != NULL)
                {
                    freeSamples3D(area_light_ptr->samples3D);
                    area_light_ptr->samples3D = NULL;
                }
            }else if(sl->light_types[i] == ENVLIGHT)
            {
                EnvLight* env_light_ptr = (EnvLight*)(sl->light_ptrs[i]);
                if(env_light_ptr->samples3D != NULL)
                {
                    freeSamples3D(env_light_ptr->samples3D);
                    env_light_ptr->samples3D = NULL;
                }
            }
            free(sl->light_ptrs[i]);
            sl->light_ptrs[i] = NULL;
        }
    }
}

// NOTE: try using function pointer to condense the code?
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
        case AABOX:
            tmp_t = rayIntersectAABox(&tmp_sr, (AABox*)scene_objects->obj_ptrs[i], ray);
            break;
        case TRIANGLE:
            tmp_t = rayIntersectTriangle(&tmp_sr, (Triangle*)scene_objects->obj_ptrs[i], ray);
            break;
        case OPENCYLINDER:
            tmp_t = rayIntersectOpenCylinder(&tmp_sr, (OpenCylinder*)scene_objects->obj_ptrs[i], ray);            
            break;
        case DISK:
            tmp_t = rayIntersectDisk(&tmp_sr, (Disk*)scene_objects->obj_ptrs[i], ray);
            break;
        case TORUS:
            tmp_t = rayIntersectTorus(&tmp_sr, (Torus*)scene_objects->obj_ptrs[i], ray);
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

float shadowIntersectTest(const SceneObjects *scene_objects, const Ray shadow_ray)
{
    float t = TMAX;

    for(int i = 0; i < scene_objects->num_obj; i++)
    {
        switch(scene_objects->obj_types[i])
        {
        case SPHERE:
            t = shadowRayIntersectSphere((Sphere*)scene_objects->obj_ptrs[i], shadow_ray);
            break;
        case PLANE:
            t = shadowRayIntersectPlane((Plane*)scene_objects->obj_ptrs[i], shadow_ray);
            break;
        case RECTANGLE:
            t = shadowRayIntersectRect((Rectangle*)scene_objects->obj_ptrs[i], shadow_ray);
            break;
        case AABOX:
            t = shadowRayIntersectAABox((AABox*)scene_objects->obj_ptrs[i], shadow_ray);
            break;
        case TRIANGLE:
            t = shadowRayIntersectTriangle((Triangle*)scene_objects->obj_ptrs[i], shadow_ray);
            break;
        case OPENCYLINDER:
            t = shadowRayIntersectOpenCylinder((OpenCylinder*)scene_objects->obj_ptrs[i], shadow_ray);            
            break;
        case DISK:
            t = shadowRayIntersectDisk((Disk*)scene_objects->obj_ptrs[i], shadow_ray);
            break;
        case TORUS:
            t = shadowRayIntersectTorus((Torus*)scene_objects->obj_ptrs[i], shadow_ray);
            break;
        }
        if(t < TMAX)
        {
            return t;
        }
    }
    return t;
}



