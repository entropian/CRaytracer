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
    Object_t obj = {COMPOUND, sc};
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
    float theta = (float)PI / 4.0f;    
    vec3 scaling = {2.0f, 2.0f, 2.0f};    
    defaultInvTransform(is->inv_transform, scaling, rotate_axis, theta, translation);

    Object_t obj = {INSTANCED, is};
    SceneObjects_push_obj(so, obj);
}

void generateMeshTriangles(SceneObjects*so, const MeshEntry mesh_entry, const SceneMaterials *sm,
                           const SceneMeshes* s_meshes)
{
    mat3 rotation;
    eulerAngToMat3(rotation, mesh_entry.orientation);
    
    for(int i = 0; i < s_meshes->size; i++)
    {
        if(strcmp(mesh_entry.mesh_name, s_meshes->meshes[i].mesh_name) == 0)
        {
            OBJShape* mesh = &(s_meshes->meshes[i]);
            Material* mat = findMaterial(mesh_entry.mat_name, sm);
            for(int j = 0; j < mesh->num_indices; j += 3)
            {
                MeshTriangle* mesh_tri = (MeshTriangle*)malloc(sizeof(MeshTriangle));
                mesh_tri->i0 = mesh->indices[j];
                mesh_tri->i1 = mesh->indices[j + 1];
                mesh_tri->i2 = mesh->indices[j + 2];
                mesh_tri->shadow = mesh_entry.shadow;
                mesh_tri->mesh_ptr = mesh;

                vec3 v0, v1, v2;
                int index;
                index = mesh_tri->i0 * 3;
                vec3_assign(v0, mesh->positions[index], mesh->positions[index+1],
                            mesh->positions[index+2]);
                index = mesh_tri->i1 * 3;    
                vec3_assign(v1, mesh->positions[index], mesh->positions[index+1],
                            mesh->positions[index+2]);
                index = mesh_tri->i2 * 3;
                vec3_assign(v2, mesh->positions[index], mesh->positions[index+1],
                            mesh->positions[index+2]);

                vec3_mult(v0, mesh_entry.scaling, v0);
                vec3_mult(v1, mesh_entry.scaling, v1);
                vec3_mult(v2, mesh_entry.scaling, v2);

                vec3 new_v0, new_v1, new_v2;
                mat3_mult_vec3(new_v0, rotation, v0);
                mat3_mult_vec3(new_v1, rotation, v1);
                mat3_mult_vec3(new_v2, rotation, v2);                
                
                vec3_add(new_v0, mesh_entry.location, new_v0);
                vec3_add(new_v1, mesh_entry.location, new_v1);
                vec3_add(new_v2, mesh_entry.location, new_v2);
                
                vec3_copy(mesh_tri->v0, new_v0);
                vec3_copy(mesh_tri->v1, new_v1);
                vec3_copy(mesh_tri->v2, new_v2);

                if(mesh->num_normals == 0)
                {
                    calcTriangleNormal(mesh_tri->normal, new_v0, new_v1, new_v2);
                }else
                {
                    // TODO: fix this for rotated mesh
                    int index = mesh->indices[j] * 3;
                    vec3_assign(mesh_tri->normal, mesh->normals[index],
                                mesh->normals[index+1], mesh->normals[index+2]);
                }
                

                mesh_tri->mat = mat;                                    
                Object_t obj = {MESH_TRIANGLE, mesh_tri};
                SceneObjects_push_obj(so, obj);
            }
        }
    }
}

void initSceneObjects(SceneObjects *so, SceneMaterials *sm, SceneMeshes* s_meshes,
                      const SceneLights* sl, const char* scenefile)
{
    FILE* fp;
#ifdef _MSC_VER
    fopen_s(&fp, scenefile, "r");
#else
    fp = fopen(scenefile, "r");
#endif
    if(fp)
    {
        int num_mat = parseMaterials(sm, fp);
        char buffer[128];
        while(getNextTokenInFile(buffer, fp))
        {
            if(buffer[0] == '#')
            {
                while(fgetc(fp) != '\n'){}
            }
            if(strcmp(buffer, "OBJECT") == 0)
            {
                getNextTokenInFile(buffer, fp);
                if(strcmp(buffer, "MESH") == 0)
                {
                    OBJShape* shapes;
                    int num_mesh;
                    char mesh_file_names[MAX_MESH][NAME_LENGTH];
                    int num_file_names = 0;
                    MeshEntry mesh_entry;
                    num_mesh = parseMesh(&mesh_entry, &shapes, mesh_file_names, &num_file_names, fp);
                    if(num_mesh > 0)
                    {
                        for(int i = 0; i < num_mesh; i++)
                        {
                            SceneMeshes_push(s_meshes, shapes[i]);
                        }
                    }
                    free(shapes);
                    if(num_mesh != -1)
                    {
                        generateMeshTriangles(so, mesh_entry, sm, s_meshes);
                    }
                }else
                {
                    Object_t obj;                    
                    if(parsePrimitive(&obj, fp, sm, buffer))
                    {
                        SceneObjects_push_obj(so, obj);
                    }
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
        Object_t obj = {SPHERE, sphere};
        SceneObjects_push_obj(so, obj);
    }
    */
    for(int i = 0; i < sl->num_lights; i++)
    {
        if(sl->light_types[i] == AREALIGHT)
        {
            AreaLight* area_light_ptr = (AreaLight*)(sl->light_ptrs[i]);
            Object_t obj = {area_light_ptr->obj_type, area_light_ptr->obj_ptr};
            SceneObjects_push_obj(so, obj);
        }
    }
}

void initAreaLights(SceneLights* sl)
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
    prepSample2DStruct(samples);
    genMultijitteredSamples(samples);

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
    vec3_assign(sphere->center, -4.0f, 5.0f, 1.0f);
    sphere->radius = 1.0f;
    sphere->min_theta = 0.0f;
    sphere->max_theta = (float)PI;
    sphere->phi = (float)PI;
    // TODO: fix the material pointer situation
    sphere->mat = (Material*)malloc(sizeof(Material));
    vec3_copy(sphere->mat->ce, area_light_ptr->color);    
    sphere->mat->ke = area_light_ptr->intensity;    
    sphere->mat->mat_type = EMISSIVE;

    Samples3D* samples = genHemisphereSamples(MULTIJITTERED, 1.0f);    

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

void initEnvLight(SceneLights* sl)
{
    if(sl->num_lights == MAX_LIGHTS){return;}
    EnvLight* env_light = (EnvLight*)malloc(sizeof(EnvLight));
    env_light->intensity = 1.0f;
    vec3_copy(env_light->color, WHITE);

    Samples3D* samples = genHemisphereSamples(MULTIJITTERED, 1.0f);

    env_light->samples3D = samples;
    sl->shadow[sl->num_lights] = true;
    sl->light_ptrs[sl->num_lights] = env_light;
    sl->light_types[sl->num_lights] = ENVLIGHT;
    (sl->num_lights)++;
    sl->env_light = env_light;
}

void initAmbLight(SceneLights *sl)
{
    sl->amb_light = (AmbientLight*)malloc(sizeof(AmbientLight));
    vec3_copy(sl->amb_light->color, WHITE);
    sl->amb_light->intensity = 0.2f;
    sl->amb_light->amb_occlusion = false;
}

void initBackgroundColor(SceneLights* sl)
{
    if(sl->env_light != NULL)
    {
        getIncRadiance(sl->bg_color, ENVLIGHT, sl->env_light);
    }else
    {
        vec3_copy(sl->bg_color, WHITE);
    }
}

void initSceneLights(SceneLights* sl)
{
    sl->env_light = NULL;    
    for(int i = 0; i < MAX_LIGHTS; i++)
    {
        sl->light_ptrs[i] = NULL;
    }
    sl->num_lights = 0;
    // Directional light

    if(sl->num_lights == MAX_LIGHTS){return;}
    DirLight* dir_light_ptr = (DirLight*)malloc(sizeof(DirLight));
    float intensity = 2.0f;
    vec3 direction = {10.0f, 10.0f, 10.0f};
    vec3_normalize(direction, direction);
    assignDirLight(dir_light_ptr, intensity, WHITE, direction);
    sl->shadow[sl->num_lights] = true;
    sl->light_ptrs[sl->num_lights] = dir_light_ptr;
    sl->light_types[sl->num_lights] = DIRECTIONAL;
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
    sl->light_ptrs[sl->num_lights] = POINTLIGHT;    
    (sl->num_lights)++;
    */


    //initAreaLights(sl);
    //initEnvLight(sl);
    initAmbLight(sl);
    initBackgroundColor(sl);
}

void initScene(SceneObjects* so, SceneLights* sl, SceneMaterials* sm, SceneMeshes* s_meshes,
               const char* scenefile)
{
    // NOTE: initSceneObjects must be called after after initSceneLights if there are area lights
    initSceneLights(sl);    
    initSceneObjects(so, sm, s_meshes, sl, scenefile);
    mvNonGridObjToStart(so);
    printf("num_obj %d\n", so->num_obj);
    printf("non grid obj %d\n", so->num_non_grid_obj);

    so->accel = GRID;
    if(so->accel == GRID)
    {
        UniformGrid* rg = UniformGrid_create(so->objects, &(so->num_obj), so->num_non_grid_obj, 2);
        so->accel_ptr = rg;
    }else if(so->accel == BVH)
    {
        BVHNode* tree;
        BVH_build(&tree, &(so->objects[so->num_non_grid_obj]), so->num_obj - so->num_non_grid_obj, 0);
        //int leaf_count = 0;
        //printBVH(tree, &leaf_count, 0);        
        so->accel_ptr = tree;
    }
}
