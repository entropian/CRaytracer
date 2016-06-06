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
#include "objloader/objloader.h"

void initSolidCylinder(SceneObjects* so)
{
    if(so->num_obj == MAX_OBJECTS){return;}
    Material mat;
    initDefaultPhongMat(&mat, YELLOW);
    float radius = 1.0f, half_height = 0.5f, phi = (float)PI;
    SolidCylinder* sc = initSolidCylinder(radius, half_height, phi, &mat, true);
    Object_t obj = {sc, COMPOUND};
    SceneObjects_push_obj(so, obj);
}

void initInstanced(SceneObjects* so, SceneMaterials* sm)
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
    initDefaultPhongMat(&mat, YELLOW);
    torus->mat = SceneMaterials_push(sm, mat, "torus_mat1");
    calcAABBTorus(torus);    
    is->obj.ptr = torus;
    is->obj.type = TORUS;    

    vec3 translation = {2.0f, 0.0f, 0.0f};
    vec3 rotate_axis = {1.0f, 1.0f, -1.0f};
    vec3_normalize(rotate_axis, rotate_axis);
    float theta = PI / 4.0f;    
    vec3 scaling = {2.0f, 2.0f, 2.0f};    
    defaultInvTransform(is->inv_transform, scaling, rotate_axis, theta, translation);

    Object_t obj = {is, INSTANCED};
    SceneObjects_push_obj(so, obj);
}

void initSceneObjects(SceneObjects *so, const SceneLights *sl, SceneMaterials *sm, const char* scenefile)
{
#ifdef _MSC_VER

#else
    FILE* fp = fopen(scenefile, "r");
#endif
    if(fp)
    {
        parseMaterials(sm, fp);
        int num_mat = sm->size;
        char buffer[128];
        while(getNextTokenInFile(buffer, fp))
        {
            if(buffer[0] == '#')
            {
                while(fgetc(fp) != '\n'){}
            }
            if(strcmp(buffer, "OBJECT") == 0)
            {
                Object_t obj;
                bool parse_status = false;
                getNextTokenInFile(buffer, fp);
                if(strcmp(buffer, "SPHERE") == 0)
                {
                    parse_status = parseSphereEntry(&obj, fp, sm, num_mat);
                }else if(strcmp(buffer, "PLANE") == 0)
                {
                    parse_status = parsePlaneEntry(&obj, fp, sm, num_mat);
                }else if(strcmp(buffer, "RECTANGLE") == 0)
                {
                    parse_status = parseRectEntry(&obj, fp, sm, num_mat);
                }else if(strcmp(buffer, "TRIANGLE") == 0)
                {
                    parse_status = parseTriangleEntry(&obj, fp, sm, num_mat);
                }else if(strcmp(buffer, "AABOX") == 0)
                {
                    parse_status = parseAABoxEntry(&obj, fp, sm, num_mat);
                }else if(strcmp(buffer, "OPENCYLINDER") == 0)
                {
                    parse_status = parseOpenCylEntry(&obj, fp, sm, num_mat);
                }else if(strcmp(buffer, "DISK") == 0)
                {
                    parse_status = parseDiskEntry(&obj, fp, sm, num_mat);
                }else if(strcmp(buffer, "TORUS") == 0)
                {
                    parse_status = parseTorusEntry(&obj, fp, sm, num_mat);
                }else if(strcmp(buffer, "MESH") == 0)
                {
                    /*
                      loadOBJ returns an array of OBJShapes
                      the shapes should be thrown into an array
                      an object from an obj file can have many groups
                      
                     */
                }

                if(parse_status)
                {
                    SceneObjects_push_obj(so, obj);
                }
            }   
        }
    }
    
    //initSolidCylinder(so);
    //initInstanced(so);
    /*
    for(int i = 0; i < 50; ++i)
    {
        Sphere* sphere = (Sphere*)malloc(sizeof(Sphere));
        sphere->shadow = true;
        sphere->min_theta = 0;
        sphere->max_theta = PI;
        sphere->phi = PI;
        sphere->radius = 0.5;
        float r = (float)rand() / (float)RAND_MAX;
        float g = (float)rand() / (float)RAND_MAX;
        float b = (float)rand() / (float)RAND_MAX;
        vec3 color = {r, g, b};
        initDefaultMatteMat(&(sphere->mat), color);
        float x = (float)rand() / (float)RAND_MAX * 18.0f - 9.0f;
        float y = (float)rand() / (float)RAND_MAX * 18.0f - 10.0f;
        float z = (float)rand() / (float)RAND_MAX * -10.0f - 2.0f;
        vec3_assign(sphere->center, x, y, z);
        //so->obj_ptrs[so->num_obj] = sphere;
        //so->obj_types[so->num_obj] = SPHERE;
        //(so->num_obj)++;
        Object_t obj = {sphere, SPHERE};
        SceneObjects_push_obj(so, obj);
    }
    */
    for(int i = 0; i < sl->num_lights; i++)
    {
        if(sl->light_types[i] == AREALIGHT)
        {
            if(so->num_obj < MAX_OBJECTS)
            {
                AreaLight* area_light_ptr = (AreaLight*)(sl->light_ptrs[i]);
                Object_t obj = {area_light_ptr->obj_ptr, area_light_ptr->obj_type};
                
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
    // TODO: fix the material pointer situation
    sphere->mat = (Material*)malloc(sizeof(Material));
    vec3_copy(sphere->mat->ce, area_light_ptr->color);    
    sphere->mat->ke = area_light_ptr->intensity;    
    sphere->mat->mat_type = EMISSIVE;

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
    //initAreaLights(sl, num_samples, num_sets);
    initEnvLight(sl, num_samples, num_sets);
}

void initScene(SceneObjects* so, SceneLights* sl, SceneMaterials* sm, const int num_samples, const int num_sets,
               const char* scenefile)
{
    // NOTE: initSceneObjects must be called after after initSceneLights if there are area lights
    initSceneLights(sl, num_samples, num_sets);    
    initSceneObjects(so, sl, sm, scenefile);
    mvNonGridObjToStart(so);

    so->accel = BVH;
    if(so->accel == GRID)
    {
        //UniformGrid* rg = UniformGrid_create(so->obj_ptrs, so->obj_types, &(so->num_obj), so->num_non_grid_obj, 2);
        UniformGrid* rg = UniformGrid_create(so->objects, &(so->num_obj), so->num_non_grid_obj, 2);
        so->accel_ptr = rg;
    }else if(so->accel == BVH)
    {
        BVHNode* tree;
        BVH_build(&tree, &(so->objects[so->num_non_grid_obj]), so->num_obj - so->num_non_grid_obj, 0);        
        int leaf_count = 0;
        printBVH(tree, &leaf_count, 0);        
        so->accel_ptr = tree;
    }
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
