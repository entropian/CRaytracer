#pragma once

#include "scenedata.h"
#include "../camera.h"

typedef struct Scene_s
{
    SceneObjects objects;
    SceneMaterials materials;
    SceneTextures textures;
    SceneMeshes meshes;
    SceneTransforms transforms;
    SceneLights lights;
    Camera camera;
}Scene;

Scene Scene_create()
{
    Scene scene;
    scene.objects = SceneObjects_create();
    scene.materials = SceneMaterials_create();
    scene.textures = SceneTextures_create();
    scene.meshes = SceneMeshes_create();
    scene.transforms = SceneTransforms_create();
    scene.lights = SceneLights_create();
    return scene;
}

void Scene_destroy(Scene* scene)
{
    freeSceneObjects(&(scene->objects));
    SceneMaterials_destroy(&(scene->materials));
    SceneTextures_destroy(&(scene->textures));
    SceneMeshes_destroy(&(scene->meshes));
    SceneTransforms_destroy(&(scene->transforms));
    freeSceneLights(&(scene->lights));
    Camera_destroy(&(scene->camera));
}


void Scene_addObject(Scene* scene, const Object_t* obj)
{
    SceneObjects_push_obj(&(scene->objects), obj);
}

Material* Scene_addMaterial(Scene* scene, const Material* mat)
{
    return SceneMaterials_push(&(scene->materials), mat);
}

Material* Scene_findMaterial(Scene* scene, const char* name)
{
    return findMaterial(name, &(scene->materials));
}

Texture* Scene_addTexture(Scene* scene, const Texture* tex, const char* name)
{
    return SceneTextures_push(&(scene->textures), tex, name);
}

Texture* Scene_findTexture(Scene* scene, const char *name)
{
    return findTexture(name, &(scene->textures));
}

Mesh* Scene_addMesh(Scene* scene, const Mesh* mesh)
{
    return SceneMeshes_push(&(scene->meshes), mesh);
}

mat3* Scene_addTransform(Scene* scene, const mat3* mat)
{
    return SceneTransforms_push(&(scene->transforms), mat);
}

Mesh** Scene_findMeshes(int* num_meshes, const Scene* scene, const char* name)
{
    return findMeshes(num_meshes, &(scene->meshes), name);
}

int Scene_getNumMaterials(Scene* scene)
{
    return scene->materials.size;
}

void Scene_printMaterials(Scene *scene)
{
    SceneMaterials *sm = &(scene->materials);
    for(int i = 0; i < sm->size; i++)
    {
        printf("%s\n", sm->materials[i].name);
        //printMaterial(&(sm->materials[i]));
    }
}
