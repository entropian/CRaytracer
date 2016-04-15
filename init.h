#pragma once

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include "util/vec.h"
#include "util/constants.h"
#include "shapes/shapes.h"
#include "gl/glcode.h"
#include "sampling.h"
#include "camera.h"
#include "lights.h"
#include "sceneobj.h"
#include "shading.h"

void initSpheres(SceneObjects *so)
{
    // TODO: tidy up material assignment
    if(so->num_obj == MAX_OBJECTS){return;}    
    Sphere* sphere = (Sphere*)malloc(sizeof(Sphere));
    //vec3_assign(sphere->center, 2.0f, 0.0f, -4.0f);
    vec3_assign(sphere->center, 0.0f, 2.0f, 0.0f);
    sphere->radius = .0f;
    sphere->shadow = true;        
    vec3_copy(sphere->mat.cd, PINK);
    vec3_copy(sphere->mat.ca, PINK);
    vec3_copy(sphere->mat.cs, PINK);    
    sphere->mat.kd = 0.6f;
    sphere->mat.ka = 1.0f;
    sphere->mat.ks = 0.4f;
    sphere->mat.exp = 10.0f;        
    sphere->mat.mat_type = PHONG;
    sphere->mat.shadow = true;
    so->obj_ptrs[so->num_obj] = sphere; 
    so->obj_types[so->num_obj] = SPHERE;
    (so->num_obj)++;

    if(so->num_obj == MAX_OBJECTS){return;}    
    sphere = (Sphere*)malloc(sizeof(Sphere));
    vec3_assign(sphere->center, 2.0f, 0.0f, -8.0f);
    sphere->radius = 1.0f;
    sphere->shadow = true;            
    vec3_copy(sphere->mat.cd, CYAN);
    vec3_copy(sphere->mat.ca, CYAN);
    vec3_copy(sphere->mat.cs, CYAN);    
    sphere->mat.kd = 0.6f;
    sphere->mat.ka = 1.0f;
    sphere->mat.ks = 0.4f;
    sphere->mat.exp = 10.0f;            
    sphere->mat.mat_type = PHONG;
    sphere->mat.shadow = true;
    so->obj_ptrs[so->num_obj] = sphere; 
    so->obj_types[so->num_obj] = SPHERE;
    (so->num_obj)++;    
}

void initPlanes(SceneObjects *so)
{
    if(so->num_obj == MAX_OBJECTS){return;}    
    Plane* plane = (Plane*)malloc(sizeof(Plane));
    plane->shadow  = true;    
    vec3_assign(plane->point, 0.0f, -1.0f, 0.0f);
    vec3_copy(plane->normal, UP);
    vec3_copy(plane->mat.cd, GREY);
    vec3_copy(plane->mat.ca, GREY);
    plane->mat.kd = 0.6f;
    plane->mat.ka = 1.0f;
    plane->mat.mat_type  = MATTE;
    plane->mat.shadow = true;    
    so->obj_ptrs[so->num_obj] = plane; 
    so->obj_types[so->num_obj] = PLANE;
    (so->num_obj)++;    
}

void initRectangles(SceneObjects *so)
{
    if(so->num_obj == MAX_OBJECTS){return;}
    Rectangle* rect = (Rectangle*)malloc(sizeof(Rectangle));
    rect->shadow = true;
    vec3_assign(rect->point, -2.0f, -1.0f, -6.0f);
    vec3_assign(rect->width, 2.0f, 0.0f, 0.0f);
    vec3_assign(rect->height, 0.0f, 4.0f, 0.0f);    
    vec3_copy(rect->normal, BACKWARD);
    vec3_copy(rect->mat.cd, YELLOW);
    vec3_copy(rect->mat.ca, YELLOW);
    rect->mat.kd = 0.6f;
    rect->mat.ka = 1.0f;
    rect->mat.mat_type  = MATTE;
    rect->mat.shadow = true;
    so->obj_ptrs[so->num_obj] = rect; 
    so->obj_types[so->num_obj] = RECTANGLE;
    (so->num_obj)++;        
}

void initTriangle(SceneObjects *so)
{
    if(so->num_obj == MAX_OBJECTS){return;}
    Triangle* triangle = (Triangle*)malloc(sizeof(Triangle));
    triangle->shadow = true;
    vec3_assign(triangle->v0, -3.0f, -1.0f, -6.0f);
    vec3_assign(triangle->v1, 0.0f, -1.0f, -6.0f);
    vec3_assign(triangle->v2, -1.0f, 2.0f, -6.0f);    
    calcTriangleNormal(triangle);
    vec3_copy(triangle->mat.cd, YELLOW);
    vec3_copy(triangle->mat.ca, YELLOW);
    triangle->mat.kd = 0.6f;
    triangle->mat.ka = 1.0f;
    triangle->mat.mat_type  = MATTE;
    triangle->mat.shadow = true;
    so->obj_ptrs[so->num_obj] = triangle; 
    so->obj_types[so->num_obj] = TRIANGLE;
    (so->num_obj)++;        
}

void initAABox(SceneObjects* so)
{
    if(so->num_obj == MAX_OBJECTS){return;}
    AABox* aabox = (AABox*)malloc(sizeof(AABox));
    aabox->shadow = true;
    vec3_assign(aabox->min, -1.0f, 1.0f, -4.0f);
    vec3_assign(aabox->max, -0.0f, 3.0f, -5.0f);        
    vec3_copy(aabox->mat.cd, YELLOW);
    vec3_copy(aabox->mat.ca, YELLOW);
    vec3_copy(aabox->mat.cs, YELLOW);    
    aabox->mat.kd = 0.6f;
    aabox->mat.ka = 1.0f;
    aabox->mat.ks = 0.4f;
    aabox->mat.exp = 10.0f;
    aabox->mat.mat_type  = PHONG;
    aabox->mat.shadow = true;
    so->obj_ptrs[so->num_obj] = aabox; 
    so->obj_types[so->num_obj] = AABOX;
    (so->num_obj)++;            
}

void initOpenCylinder(SceneObjects* so)
{
    if(so->num_obj == MAX_OBJECTS){return;}
    OpenCylinder* oc = (OpenCylinder*)malloc(sizeof(OpenCylinder));
    oc->shadow = true;
    oc->half_height = .9f;
    oc->radius = 1.0f;
    oc->phi = 2.0f * PI;
    oc->normal_type = OPEN;
    vec3_copy(oc->mat.cd, YELLOW);
    vec3_copy(oc->mat.ca, YELLOW);
    vec3_copy(oc->mat.cs, YELLOW);    
    oc->mat.kd = 0.6f;
    oc->mat.ka = 1.0f;
    oc->mat.ks = 0.4f;
    oc->mat.exp = 10.0f;
    oc->mat.mat_type  = PHONG;
    oc->mat.shadow = true;
    so->obj_ptrs[so->num_obj] = oc; 
    so->obj_types[so->num_obj] = OPENCYLINDER;
    (so->num_obj)++;      
}

void initDisk(SceneObjects* so)
{
    if(so->num_obj == MAX_OBJECTS){return;}
    Disk* disk = (Disk*)malloc(sizeof(Disk));
    disk->shadow = true;
    vec3_assign(disk->center, 0.0f, 0.0f, 0.0f);
    vec3_assign(disk->normal, 0.0f, 0.0f, -1.0f);
    disk->radius = 1.0f;
    vec3_copy(disk->mat.cd, YELLOW);
    vec3_copy(disk->mat.ca, YELLOW);
    vec3_copy(disk->mat.cs, YELLOW);    
    disk->mat.kd = 0.6f;
    disk->mat.ka = 1.0f;
    disk->mat.ks = 0.4f;
    disk->mat.exp = 10.0f;
    disk->mat.mat_type  = PHONG;
    disk->mat.shadow = true;
    so->obj_ptrs[so->num_obj] = disk; 
    so->obj_types[so->num_obj] = DISK;
    (so->num_obj)++;            
}

void initTorus(SceneObjects* so)
{
    if(so->num_obj == MAX_OBJECTS){return;}
    Torus* torus = (Torus*)malloc(sizeof(Torus));
    torus->shadow = true;
    torus->swept_radius = 2.0f;
    torus->tube_radius = 0.5f;
    vec3_copy(torus->mat.cd, YELLOW);
    vec3_copy(torus->mat.ca, YELLOW);
    vec3_copy(torus->mat.cs, YELLOW);    
    torus->mat.kd = 0.6f;
    torus->mat.ka = 1.0f;
    torus->mat.ks = 0.4f;
    torus->mat.exp = 10.0f;
    torus->mat.mat_type  = PHONG;
    torus->mat.shadow = true;
    calcAABBTorus(torus);
    so->obj_ptrs[so->num_obj] = torus; 
    so->obj_types[so->num_obj] = TORUS;
    (so->num_obj)++;            
}

void initSceneObjects(SceneObjects *so, const SceneLights *sl)
{
    for(int i = 0; i < MAX_OBJECTS; i++)
    {
        so->obj_ptrs[i] = NULL;
    }
    so->num_obj = 0;
    //initSpheres(so);
    initPlanes(so);
    //initRectangles(so);
    //initAABox(so);
    //initTriangle(so);
    //initOpenCylinder(so);
    //initDisk(so);
    initTorus(so);
    
    for(int i = 0; i < sl->num_lights; i++)
    {
        if(sl->light_types[i] == AREALIGHT)
        {
            if(so->num_obj < MAX_OBJECTS)
            {
                so->obj_ptrs[so->num_obj] = ((AreaLight*)(sl->light_ptrs[i]))->obj_ptr;
                so->obj_types[so->num_obj] = ((AreaLight*)(sl->light_ptrs[i]))->obj_type;
                (so->num_obj)++;
            }
        }
    }
}

void initAreaLights(SceneLights* sl, const int num_samples, const int num_sets)
{
    // Area light
    if(sl->num_lights == MAX_LIGHTS){return;}
    AreaLight* area_light_ptr = (AreaLight*)malloc(sizeof(AreaLight));
    area_light_ptr->intensity = 5.0f;
    vec3_copy(area_light_ptr->color, WHITE);
    vec3_assign(area_light_ptr->sample_point, 0.0f, 0.0f, 0.0f);

    /*
    // Rectangle
    Rectangle* rect = (Rectangle*)malloc(sizeof(Rectangle));
    rect->shadow = false;
    vec3_assign(rect->point, -2.0f, 0.0f, -6.0f);
    vec3_assign(rect->width, 2.0f, 0.0f, 0.0f);
    vec3_assign(rect->height, 0.0f, 4.0f, 0.0f);    
    vec3_copy(rect->normal, BACKWARD);
    vec3_copy(rect->mat.ce, area_light_ptr->color);
    rect->mat.ke = area_light_ptr->intensity;    
    rect->mat.mat_type = EMISSIVE;
    float width = sqrt(vec3_dot(rect->width, rect->width));
    float height = sqrt(vec3_dot(rect->height, rect->height));        
    
    Samples2D* samples = (Samples2D*)malloc(sizeof(Samples2D));
    samples->samples = NULL;
    samples->shuffled_indices = NULL;
    prepSample2DStruct(samples, num_samples, num_sets);
    genMultijitteredSamples(samples, num_samples, num_sets);

    area_light_ptr->pdf = 1.0f/(width * height);
    area_light_ptr->samples2D = samples;
    area_light_ptr->samples3D = NULL;
    area_light_ptr->obj_ptr = rect;
    area_light_ptr->obj_type = RECTANGLE;
    sl->shadow[sl->num_lights] = true;    
    sl->light_ptrs[sl->num_lights] = area_light_ptr;
    sl->light_types[sl->num_lights] = AREALIGHT;
    (sl->num_lights)++;
    */
    // Sphere
    Sphere* sphere = (Sphere*)malloc(sizeof(Sphere));
    sphere->shadow = false;
    vec3_assign(sphere->center, -2.0f, 2.0f, -3.0f);
    sphere->radius = 1.0f;
    vec3_copy(sphere->mat.ce, area_light_ptr->color);    
    sphere->mat.ke = area_light_ptr->intensity;    
    sphere->mat.mat_type = EMISSIVE;

    Samples3D* samples = genHemisphereSamples(MULTIJITTERED, num_samples, num_sets);    

    area_light_ptr->pdf = 1.0f / (4.0f * (float)PI * sphere->radius * sphere->radius);
    area_light_ptr->samples2D = NULL;
    area_light_ptr->samples3D = samples;
    area_light_ptr->obj_ptr = sphere;
    area_light_ptr->obj_type = SPHERE;
    sl->shadow[sl->num_lights] = true;    
    sl->light_ptrs[sl->num_lights] = area_light_ptr;
    sl->light_types[sl->num_lights] = AREALIGHT;
    (sl->num_lights)++;
}

void initEnvLight(SceneLights*sl, const int num_samples, const int num_sets)
{
    if(sl->num_lights == MAX_LIGHTS){return;}
    EnvLight* env_light_ptr = (EnvLight*)malloc(sizeof(EnvLight));
    env_light_ptr->intensity = 1.0f;
    vec3_copy(env_light_ptr->color, WHITE);

    Samples3D* samples = genHemisphereSamples(MULTIJITTERED, num_samples, num_sets);

    env_light_ptr->samples3D = samples;
    sl->shadow[sl->num_lights] = true;
    sl->light_ptrs[sl->num_lights] = env_light_ptr;
    sl->light_types[sl->num_lights] = ENVLIGHT;
    (sl->num_lights)++;
    sl->env_light_ptr = env_light_ptr;
}

void initSceneLights(SceneLights* sl, const int num_samples, const int num_sets)
{
    for(int i = 0; i < MAX_LIGHTS; i++)
    {
        sl->light_ptrs[i] = NULL;
    }
    sl->num_lights = 0;
    // Directional light
    if(sl->num_lights == MAX_LIGHTS){return;}
    DirLight* dir_light_ptr = (DirLight*)malloc(sizeof(DirLight));
    float intensity = 0.0f;
    vec3 direction = {10.0f, 10.0f, 10.0f};
    vec3_normalize(direction, direction);
    assignDirLight(dir_light_ptr, intensity, WHITE, direction);
    sl->shadow[sl->num_lights] = true;
    sl->light_ptrs[sl->num_lights] = dir_light_ptr;
    (sl->num_lights)++;
    
    // Point light
    if(sl->num_lights == MAX_LIGHTS){return;}
    PointLight* point_light_ptr = (PointLight*)malloc(sizeof(PointLight));
    intensity = 0.1f;
    vec3 point = {-3.0f, -5.0f, 0.0f};
    assignPointLight(point_light_ptr, intensity, WHITE, point);
    sl->shadow[sl->num_lights] = true;    
    sl->light_ptrs[sl->num_lights] = point_light_ptr;
    (sl->num_lights)++;

    sl->env_light_ptr = NULL;
    //initAreaLights(sl, num_samples, num_sets);
    initEnvLight(sl, num_samples, num_sets);
}

void initScene(SceneObjects* so, SceneLights* sl, const int num_samples, const int num_sets)
{
    // NOTE: initSceneObjects must be called after after initSceneLights if there are area lights
    initSceneLights(sl, num_samples, num_sets);    
    initSceneObjects(so, sl);        
}

void initBackgroundColor(vec3 r, const SceneLights* sl, const vec3 default_color)
{
    if(sl->env_light_ptr != NULL)
    {
        getIncRadiance(r, ENVLIGHT, sl->env_light_ptr);
    }else
    {
        vec3_copy(r, default_color);
    }
}
