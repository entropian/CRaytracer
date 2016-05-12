#pragma once

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <ctime>
#include "util/vec.h"
#include "util/constants.h"
#include "shapes/shapes.h"
#include "shapes/instanced.h"
#include "gl/glcode.h"
#include "sampling.h"
#include "camera.h"
#include "lights.h"
#include "sceneobj.h"
#include "scenefile.h"
#include "shading.h"

void initSolidCylinder(SceneObjects* so)
{
    if(so->num_obj == MAX_OBJECTS){return;}
    Material mat;
    initDefaultPhongMat(&mat, YELLOW);
    float radius = 1.0f, half_height = 0.5f, phi = (float)PI;
    SolidCylinder* sc = initSolidCylinder(radius, half_height, phi, &mat, true);
    so->obj_ptrs[so->num_obj] = sc;
    so->obj_types[so->num_obj] = COMPOUND;
    (so->num_obj)++;    
}

void initInstanced(SceneObjects* so)
{
    if(so->num_obj == MAX_OBJECTS){return;}
    Material mat;
    InstancedShape* is = (InstancedShape*)malloc(sizeof(InstancedShape));

    /*
    initDefaultPhongMat(&mat, YELLOW);
    float radius = 1.0f, half_height = 0.5f, phi = (float)PI;
    SolidCylinder* sc = initSolidCylinder(radius, half_height, phi, &mat, true);
    is->obj_ptr = sc;
    is->obj_type = SOLIDCYLINDER;
    */


    Torus* torus = (Torus*)malloc(sizeof(Torus));
    torus->shadow = true;
    torus->swept_radius = 2.0f;
    torus->tube_radius = 0.5f;
    torus->phi = (float)PI;
    initDefaultPhongMat(&(torus->mat), YELLOW);    
    calcAABBTorus(torus);    
    is->obj_ptr = torus;
    is->obj_type = TORUS;    

    vec3 translation = {2.0f, 0.0f, 0.0f};
    vec3 rotate_axis = {1.0f, 1.0f, -1.0f};
    vec3_normalize(rotate_axis, rotate_axis);
    float theta = PI / 4.0f;    
    vec3 scaling = {2.0f, 2.0f, 2.0f};    
    defaultInvTransform(is->inv_transform, scaling, rotate_axis, theta, translation);

    so->obj_ptrs[so->num_obj] = is;
    so->obj_types[so->num_obj] = INSTANCED;
    (so->num_obj)++;
}

void initSceneObjects(SceneObjects *so, const SceneLights *sl)
{
    for(int i = 0; i < MAX_OBJECTS; i++)
    {
        so->obj_ptrs[i] = NULL;
    }
    so->num_obj = 0;    
    FILE* fp = fopen("scene.txt", "r");
    char buffer[128];
    while(getNextTokenInFile(buffer, fp))
    {
        if(strcmp(buffer, "OBJECT") == 0)
        {
            getNextTokenInFile(buffer, fp);
            if(strcmp(buffer, "SPHERE") == 0)
            {
                Sphere* sphere_ptr;
                if(parseSphereEntry(&sphere_ptr, fp))
                {
                    so->obj_ptrs[so->num_obj] = sphere_ptr;
                    so->obj_types[so->num_obj] = SPHERE;
                    (so->num_obj)++;
                }
            }else if(strcmp(buffer, "PLANE") == 0)
            {
                Plane* plane_ptr;
                if(parsePlaneEntry(&plane_ptr, fp))
                {
                    so->obj_ptrs[so->num_obj] = plane_ptr;
                    so->obj_types[so->num_obj] = PLANE;
                    (so->num_obj)++;
                }      
            }else if(strcmp(buffer, "RECTANGLE") == 0)
            {
                Rectangle* rect_ptr;
                if(parseRectEntry(&rect_ptr, fp))
                {
                    so->obj_ptrs[so->num_obj] = rect_ptr;
                    so->obj_types[so->num_obj] = RECTANGLE;
                    (so->num_obj)++;
                }      
            }else if(strcmp(buffer, "TRIANGLE") == 0)
            {
                Triangle* tri_ptr;
                if(parseTriangleEntry(&tri_ptr, fp))
                {
                    so->obj_ptrs[so->num_obj] = tri_ptr;
                    so->obj_types[so->num_obj] = TRIANGLE;
                    (so->num_obj)++;
                }
            }else if(strcmp(buffer, "OPENCYLINDER") == 0)
            {
                OpenCylinder* cyl_ptr;
                if(parseOpenCylEntry(&cyl_ptr, fp))
                {
                    so->obj_ptrs[so->num_obj] = cyl_ptr;
                    so->obj_types[so->num_obj] = OPENCYLINDER;
                    (so->num_obj)++;
                }
            }else if(strcmp(buffer, "DISK") == 0)
            {
                Disk* disk_ptr;
                if(parseDiskEntry(&disk_ptr, fp))
                {
                    so->obj_ptrs[so->num_obj] = disk_ptr;
                    so->obj_types[so->num_obj] = DISK;
                    (so->num_obj)++;
                }
            }else if(strcmp(buffer, "TORUS") == 0)
            {
                Torus* torus_ptr;
                if(parseTorusEntry(&torus_ptr, fp))
                {
                    so->obj_ptrs[so->num_obj] = torus_ptr;
                    so->obj_types[so->num_obj] = TORUS;
                    (so->num_obj)++;
                }else
                {
                    printf("falsed\n");
                }
            }
        }   
    }
    
    //initTorus(so);
    //initSolidCylinder(so);
    //initInstanced(so);
    //initSolidCylinder(so);
    
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
    vec3_assign(sphere->center, -2.0f, 3.0f, -2.0f);
    sphere->radius = 1.0f;
    sphere->min_theta = 0.0f;
    sphere->max_theta = (float)PI;
    sphere->phi = (float)PI;
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
    /*
    if(sl->num_lights == MAX_LIGHTS){return;}
    PointLight* point_light_ptr = (PointLight*)malloc(sizeof(PointLight));
    intensity = 0.1f;
    vec3 point = {-3.0f, -5.0f, 0.0f};
    assignPointLight(point_light_ptr, intensity, WHITE, point);
    sl->shadow[sl->num_lights] = true;    
    sl->light_ptrs[sl->num_lights] = point_light_ptr;
    (sl->num_lights)++;
    */

    sl->env_light_ptr = NULL;
    initAreaLights(sl, num_samples, num_sets);
    //initEnvLight(sl, num_samples, num_sets);
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
