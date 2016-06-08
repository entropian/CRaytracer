#pragma once

#include "shapes/shapes.h"
#include "shapes/instanced.h"
#include "lights.h"
#include "accelerator/uniformgrid.h"
#include "accelerator/bvh.h"
#include "objloader/objloader.h"

enum AccelType
{
    NONE,
    GRID,
    BVH
};

typedef struct 
{
    AccelType accel;    
    int num_obj;
    int num_non_grid_obj;
    int max;
    void* accel_ptr;    
    Object_t* objects;
} SceneObjects;

SceneObjects SceneObjects_create()
{
    SceneObjects so;
    so.accel = NONE;
    so.num_obj = so.num_non_grid_obj = 0;
    so.objects = (Object_t*)malloc(sizeof(Object_t) * INITIAL_NUM_OBJECTS);
    so.max = INITIAL_NUM_OBJECTS;
    return so;
}

void SceneObjects_push_obj(SceneObjects* so, Object_t obj)
{
    DBuffer obj_buffer;
    DBuffer_assume(&obj_buffer, (char*)so->objects, so->num_obj, so->max, sizeof(Object_t));    
    DBuffer_push(obj_buffer, obj);
    so->objects = (Object_t*)(obj_buffer.data);
    so->num_obj = DBuffer_size(obj_buffer);
    so->max = DBuffer_max_elements(obj_buffer);
}

typedef struct
{
    int size;
    int max;
    Material* materials;
    char** names;
}SceneMaterials;

SceneMaterials SceneMaterials_create()
{
    SceneMaterials sm;
    sm.materials = (Material*)malloc(sizeof(Material) * DEFAULT_MATERIAL);
    sm.names = (char**)malloc(sizeof(char*) * DEFAULT_MATERIAL);
    sm.size = 0;
    sm.max = DEFAULT_MATERIAL;
    return sm;
}

Material* SceneMaterials_push(SceneMaterials* sm, const Material mat, const char* name)
{
    DBuffer mat_buffer;
    DBuffer_assume(&mat_buffer, (char*)sm->materials, sm->size, sm->max, sizeof(Material));
    DBuffer_push(mat_buffer, mat);    

    DBuffer name_buffer;
    DBuffer_assume(&name_buffer, (char*)sm->names, sm->size, sm->max, sizeof(char*));
    char* mat_name = (char*)malloc(sizeof(char) * MAX_NAME_LENGTH);
    strcpy_s(mat_name, NAME_LENGTH, name);
    DBuffer_push(name_buffer, mat_name);
    
    sm->materials = (Material*)(mat_buffer.data);
    sm->names = (char**)(name_buffer.data);
    sm->size = DBuffer_size(mat_buffer);
    sm->max = DBuffer_max_elements(mat_buffer);
    return &(sm->materials[sm->size - 1]);
}

Material* findMaterial(const char* mat_name, const SceneMaterials* sm)
{
    for(int i = 0; i < sm->size; i++)
    {
        if(strcmp(mat_name, sm->names[i]) == 0)
        {
            return &(sm->materials[i]);
        }
    }
    // TODO: return a default material
    fprintf(stderr, "Material not found\n");
    return NULL;    
}

typedef struct
{
    int size;
    int max;    
    OBJShape* meshes;
}SceneMeshes;

SceneMeshes SceneMeshes_create()
{
    SceneMeshes s_meshes;
    s_meshes.meshes = (OBJShape*)malloc(sizeof(OBJShape) * MAX_MESH);
    s_meshes.size = 0;
    s_meshes.max = MAX_MESH;
    return s_meshes;
}

OBJShape* SceneMeshes_push(SceneMeshes* s_meshes, const OBJShape obj_shape)
{
    DBuffer mesh_buffer;
    DBuffer_assume(&mesh_buffer, (char*)s_meshes->meshes, s_meshes->size, s_meshes->max, sizeof(OBJShape));
    DBuffer_push(mesh_buffer, obj_shape);

    s_meshes->meshes = (OBJShape*)(mesh_buffer.data);
    s_meshes->size = DBuffer_size(mesh_buffer);
    s_meshes->max = DBuffer_max_elements(mesh_buffer);
    return &(s_meshes->meshes[s_meshes->size - 1]);
}

typedef struct SceneLights
{
    int num_lights;
    bool shadow[MAX_LIGHTS];    
    LightType light_types[MAX_LIGHTS];
    void* light_ptrs[MAX_LIGHTS];
    EnvLight* env_light_ptr;
} SceneLights;

Material* tmp_mat = (Material*)malloc(sizeof(Material));
vec3 color = {0.4f, 0.4f, 0.4f};

void freeSceneObjects(SceneObjects* so)
{

    for(int i = 0; i < so->num_obj; i++)
    {
        if(so->objects[i].ptr != NULL)
        {
            if(so->objects[i].type == COMPOUND)
            {
                freeCompoundObject((CompoundObject*)(so->objects[i].ptr));
            }else if(so->objects[i].type == INSTANCED)
            {
                // TODO
                freeInstancedShape((InstancedShape*)(so->objects[i].ptr));                
            }else
            {
                free(so->objects[i].ptr);
                so->objects[i].ptr = NULL;
            }
        }
    }

    free(so->objects);
    switch(so->accel)
    {
    case GRID:
        UniformGrid_destroy((UniformGrid*)(so->accel_ptr));
        break;
    case BVH:
        // TODO
        break;
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

void mvNonGridObjToStart(SceneObjects *so)
{
    int j = 0;
    for(int i = 0; i < so->num_obj; ++i)
    {
        if(!isGridObjType(so->objects[i].type))
        {
            Object_t tmp_obj = so->objects[i];
            so->objects[i] = so->objects[j];
            so->objects[j] = tmp_obj;
            ++j;
        }
    }
    so->num_non_grid_obj = j;
}
