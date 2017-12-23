#pragma once

#include "../shapes/shapes.h"
#include "../shapes/instanced.h"
#include "../lights.h"
#include "../accelerator/uniformgrid.h"
#include "../accelerator/bvh.h"
#include "../accelerator/bvh4.h"
#include "../objloader/objloader.h"
#include "../texture.h"

enum AccelType
{
    NONE,
    GRID,
    BVH,
    BVH4
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

void SceneObjects_push_obj(SceneObjects* so, const Object_t* obj)
{
    DBuffer obj_buffer;
    DBuffer_assume(&obj_buffer, (char*)so->objects, so->num_obj, so->max, sizeof(Object_t));
    DBuffer_push(obj_buffer, *obj);
    so->objects = (Object_t*)(obj_buffer.data);
    so->num_obj = DBuffer_size(obj_buffer);
    so->max = DBuffer_max_elements(obj_buffer);
}

typedef struct SceneTextures
{
    int size;
    int max;
    Texture* textures;
    char** names;
}SceneTextures;

SceneTextures SceneTextures_create()
{
    SceneTextures st;
    st.textures = (Texture*)malloc(sizeof(Texture) * DEFAULT_NUM_TEXTURES);
    st.names = (char**)malloc(sizeof(char*) * DEFAULT_NUM_TEXTURES);
    st.size = 0;
    st.max = DEFAULT_NUM_TEXTURES;
    return st;
}

void SceneTextures_destroy(SceneTextures* st)
{
    for(int i = 0; i < st->size; i++)
    {
        freeTexture(&(st->textures[i]));
        free(st->names[i]);
    }
    free(st->textures);
    free(st->names);
    st->size = st->max = 0;    
}

Texture* SceneTextures_push(SceneTextures* st, const Texture* tex, const char* name)
{
    DBuffer tex_buffer;
    DBuffer_assume(&tex_buffer, (char*)st->textures, st->size, st->max, sizeof(Texture));
    DBuffer_push(tex_buffer, *tex);

    DBuffer name_buffer;
    DBuffer_assume(&name_buffer, (char*)st->names, st->size, st->max, sizeof(char*));
    char* tex_name = (char*)malloc(sizeof(char) * MAX_NAME_LENGTH);
    stringCopy(tex_name, MAX_NAME_LENGTH, name);
    DBuffer_push(name_buffer, tex_name);

    st->textures = (Texture*)(tex_buffer.data);
    st->names = (char**)(name_buffer.data);
    st->size = DBuffer_size(tex_buffer);
    st->max = DBuffer_max_elements(tex_buffer);
    return &(st->textures[st->size - 1]);
}

Texture* findTexture(const char* tex_name, const SceneTextures* st)
{
    for(int i = 0; i < st->size; i++)
    {
        if(strcmp(tex_name, st->names[i]) == 0)
        {
            return &(st->textures[i]);
        }
    }    
    //fprintf(stderr, "Texture not found.\n");
    return NULL;
}

typedef struct SceneMaterials_s
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
    initDefaultMatteMat(&(sm.materials[0]), GREY);
    sm.names[0] = (char*)malloc(sizeof(char) * NAME_LENGTH);
    //strcpy_s(sm.names[0], NAME_LENGTH, "DEFAULT");
    stringCopy(sm.names[0], NAME_LENGTH, "DEFAULT");
    sm.size = 1;
    sm.max = DEFAULT_MATERIAL;
    return sm;
}

void SceneMaterials_destroy(SceneMaterials* sm)
{
    for(int i = 0; i < sm->size; i++)
    {
        free(sm->names[i]);
    }    
    free(sm->materials);
    free(sm->names);
    sm->size = sm->max = 0;
}

Material* SceneMaterials_push(SceneMaterials* sm, const Material* mat, const char* name)
{
    DBuffer mat_buffer;
    DBuffer_assume(&mat_buffer, (char*)sm->materials, sm->size, sm->max, sizeof(Material));
    DBuffer_push(mat_buffer, *mat);    

    DBuffer name_buffer;
    DBuffer_assume(&name_buffer, (char*)sm->names, sm->size, sm->max, sizeof(char*));
    char* mat_name = (char*)malloc(sizeof(char) * MAX_NAME_LENGTH);
    //strcpy_s(mat_name, NAME_LENGTH, name);
    stringCopy(mat_name, MAX_NAME_LENGTH, name);
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
    
    //fprintf(stderr, "Material not found. Default material returned\n");
    //return &(sm->materials[0]);
    fprintf(stderr, "Material not found. Returning NULL\n");
    return NULL;
}

typedef struct SceneMeshes_s
{
    int size;
    int max;    
    Mesh* meshes;
}SceneMeshes;

SceneMeshes SceneMeshes_create()
{
    SceneMeshes s_meshes;
    s_meshes.meshes = (Mesh*)malloc(sizeof(Mesh) * MAX_MESH);
    s_meshes.size = 0;
    s_meshes.max = MAX_MESH;
    return s_meshes;
}

Mesh* SceneMeshes_push(SceneMeshes* s_meshes, const Mesh* mesh)
{
    DBuffer mesh_buffer;
    DBuffer_assume(&mesh_buffer, (char*)s_meshes->meshes, s_meshes->size, s_meshes->max, sizeof(Mesh));
    DBuffer_push(mesh_buffer, *mesh);

    s_meshes->meshes = (Mesh*)(mesh_buffer.data);
    s_meshes->size = DBuffer_size(mesh_buffer);
    s_meshes->max = DBuffer_max_elements(mesh_buffer);
    return &(s_meshes->meshes[s_meshes->size - 1]);
}

void SceneMeshes_destroy(SceneMeshes* s_meshes)
{
    for(int i = 0; i < s_meshes->size; i++)
    {
        Mesh_destroy(&(s_meshes->meshes[i]));
    }
    free(s_meshes->meshes);
    s_meshes->size = 0;
    s_meshes->max = 0;    
}

Mesh** findMeshes(int* num_meshes, const SceneMeshes* s_meshes, const char* name)
{
    DBuffer mesh_buffer = DBuffer_create(Mesh*);
    for(int i = 0; i < s_meshes->size; i++)
    {
        Mesh* mesh_ptr = &(s_meshes->meshes[i]);
        if(strcmp(name, mesh_ptr->mesh_name) == 0)
        {
            DBuffer_push(mesh_buffer, mesh_ptr);
        }
    }
    *num_meshes = DBuffer_size(mesh_buffer);
    return (Mesh**)(mesh_buffer.data);
}

typedef struct SceneTransforms_s
{
    int size;
    int max;
    mat3* matrices;
}SceneTransforms;

SceneTransforms SceneTransforms_create()
{
    SceneTransforms transforms;
    transforms.matrices = (mat3*)malloc(sizeof(mat3) * 128);
    transforms.size = 0;
    transforms.max = 128;
    return transforms;
}

mat3* SceneTransforms_push(SceneTransforms* transforms, const mat3* mat)
{
    DBuffer matrix_buffer;
    DBuffer_assume(&matrix_buffer, (char*)transforms->matrices, transforms->size, transforms->max, sizeof(mat3));
    DBuffer_push(matrix_buffer, *mat);

    transforms->matrices = (mat3*)(matrix_buffer.data);
    transforms->size = DBuffer_size(matrix_buffer);
    transforms->max = DBuffer_max_elements(matrix_buffer);
    return &(transforms->matrices[transforms->size - 1]);
}

void SceneTransforms_destroy(SceneTransforms* transforms)
{
    free(transforms->matrices);
    transforms->size = 0;
    transforms->max = 0;
}

typedef struct SceneLights
{
    int num_lights;
    bool shadow[MAX_LIGHTS];    
    LightType light_types[MAX_LIGHTS];
    void* light_ptrs[MAX_LIGHTS];
    float power[MAX_LIGHTS];
    EnvLight* env_light;
    AmbientLight* amb_light;
    vec3 bg_color;
} SceneLights;

SceneLights SceneLights_create()
{
    SceneLights sl;
    sl.num_lights = 0;
    sl.env_light = NULL;
    sl.amb_light = NULL;
    vec3_copy(sl.bg_color, WHITE);
    return sl;
}

int SceneLights_addLight(void* light_ptr, LightType light_type, SceneLights* sl)
{
    if(sl->num_lights == MAX_LIGHTS)
    {
        return 0;
    }
    sl->shadow[sl->num_lights] = true;
    sl->light_types[sl->num_lights] = light_type;
    sl->light_ptrs[sl->num_lights] = light_ptr;
    sl->num_lights++;
    return 1;
}

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
        BVH_destroy((BVHNode*)(so->accel_ptr));
        break;
    }
    
}

void freeSceneLights(SceneLights* sl)
{
    EnvLight_destroy(sl->env_light);
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
            }else if(sl->light_types[i] == POINTLIGHT)
            {
                PointLight *point_light = (PointLight*)(sl->light_ptrs[i]);
                if(point_light->proj_map)
                {
                    free(point_light->proj_map);
                }
            }
            free(sl->light_ptrs[i]);
            sl->light_ptrs[i] = NULL;
        }
    }
    free(sl->amb_light);
}

void mvNonGridObjToStart(SceneObjects *so)
{
    int j = 0;
    for(int i = 0; i < so->num_obj; ++i)
    {
        if(!isGridObjType(so->objects[i]))
        {
            Object_t tmp_obj = so->objects[i];
            so->objects[i] = so->objects[j];
            so->objects[j] = tmp_obj;
            ++j;
        }
    }
    so->num_non_grid_obj = j;
}
