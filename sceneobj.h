#pragma once

typedef struct SceneObjects
{
    Sphere spheres[MAX_SPHERES];
    int num_spheres;
    Plane planes[MAX_PLANES];
    int num_planes;
} SceneObjects;

float intersectTest(ShadeRec *sr, SceneObjects *scene_objects, const Ray ray)
{
    float tmp_t = TMAX,  min_t = TMAX;
    ShadeRec tmp_sr, min_sr;
    // Spheres
    for(int i = 0; i < scene_objects->num_spheres; i++)
    {
        tmp_t = rayIntersectSphere(&tmp_sr, &(scene_objects->spheres[i]), ray);
        if(tmp_t < min_t)
        {
            min_t = tmp_t;
            min_sr = tmp_sr;
        }
    }

    // Planes
    for(int i = 0; i < scene_objects->num_planes; i++)
    {
        tmp_t = rayIntersectPlane(&tmp_sr, &(scene_objects->planes[i]), ray);
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
    // Spheres
    for(int i = 0; i < scene_objects->num_spheres; i++)
    {
        tmp_t = shadowRayIntersectSphere(&(scene_objects->spheres[i]), shadow_ray);
        if(tmp_t < min_t)
        {
            min_t = tmp_t;
        }
    }
    // Planes
    for(int i = 0; i < scene_objects->num_planes; i++)
    {
        tmp_t = shadowRayIntersectPlane(&(scene_objects->planes[i]), shadow_ray);
        if(tmp_t < min_t)
        {
            min_t = tmp_t;
        }
    }
    return min_t;
}
