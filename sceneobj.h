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
    DBuffer objects;
} scene_obj_data;

/*
  Rethink SceneObjects
  other relevant data: mesh array
 */

typedef struct SceneObjects
{
    AccelType accel;    
    int num_obj;
    int num_mat;    
    int num_non_grid_obj;
    void* accel_ptr;    
    Object_t* objects;
    scene_obj_data* data;
} SceneObjects;

SceneObjects SceneObjects_create()
{
    SceneObjects so;
    so.accel = NONE;
    so.num_obj = so.num_mat = so.num_non_grid_obj = 0;
    so.data = (scene_obj_data*)malloc(sizeof(scene_obj_data));
    so.data->objects = DBuffer_create(Object_t);
    so.objects = (Object_t*)(so.data->objects.data);
    return so;
}

void SceneObjects_push_obj(SceneObjects* so, Object_t obj)
{
    DBuffer_push(so->data->objects, obj);
    so->objects = (Object_t*)(so->data->objects.data);
    so->num_obj = DBuffer_size(so->data->objects);
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
    DBuffer name_buffer;
    mat_buffer.data = (char*)(sm->materials);
    mat_buffer.size = sm->size * sizeof(Material);
    mat_buffer.max = sm->max * sizeof(Material);
    mat_buffer.element_size = sizeof(Material);
    DBuffer_push(mat_buffer, mat);    
    
    name_buffer.data = (char*)(sm->names);
    name_buffer.size = sm->size * sizeof(char*);
    name_buffer.max = sm->max * sizeof(char*);
    name_buffer.element_size = sizeof(char*);

    char* mat_name = (char*)malloc(sizeof(char) * MAX_NAME_LENGTH);
    strcpy(mat_name, name);
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
    OBJShape* meshes;
    int size;
    int max;
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
    mesh_buffer.data = (char*)s_meshes->meshes;
    mesh_buffer.size = s_meshes->size * sizeof(OBJShape);
    mesh_buffer.max = s_meshes->max * sizeof(OBJShape);
    mesh_buffer.element_size = sizeof(OBJShape);

    DBuffer_push(mesh_buffer, obj_shape);

    s_meshes->meshes = (OBJShape*)(mesh_buffer.data);
    s_meshes->size = DBuffer_size(mesh_buffer);
    s_meshes->max = DBuffer_max_elements(mesh_buffer);

    return &(s_meshes->meshes[s_meshes->size - 1]);
}

/* old
typedef struct SceneObjects
{
    AccelType accel;    
    int num_obj;
    int num_mat;    
    int num_non_grid_obj;
    ObjectType obj_types[MAX_OBJECTS];    // TODO: use dynamic buffer 
    void* obj_ptrs[MAX_OBJECTS];
    Material materials[MAX_MATERIAL];    
    void* accel_ptr;
} SceneObjects;
*/

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

    DBuffer_destroy(&(so->data->objects));
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
