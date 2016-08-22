#pragma once

#include "scenedata.h"

typedef struct Scene_s
{
    SceneObjects objects;
    SceneMaterials materials;
    SceneTextures textures;
    SceneMeshes meshes;
    SceneTransforms transforms;
    SceneLights lights;
}Scene;

Scene Scene_create()
{
    Scene scene;
    scene.objects = SceneObjects_create();
    scene.materials = SceneMaterials_create();
    scene.textures = SceneTextures_create();
    scene.meshes = SceneMeshes_create();
    scene.transforms = SceneTransforms_create();
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
}


void Scene_addObject(Scene* scene, const Object_t* obj)
{
    SceneObjects_push_obj(&(scene->objects), obj);
}

Material* Scene_addMaterial(Scene* scene, const Material* mat, const char* name)
{
    return SceneMaterials_push(&(scene->materials), mat, name);
}

Texture* Scene_addTexture(Scene* scene, const Texture* tex, const char* name)
{
    return SceneTextures_push(&(scene->textures), tex, name);
}

Mesh* Scene_addMesh(Scene* scene, const Mesh* mesh)
{
    return SceneMeshes_push(&(scene->meshes), mesh);
}

mat3* Scene_addTransfrom(Scene* scene, const mat3* mat)
{
    return SceneTransforms_push(&(scene->transforms), mat);
}

Material* Scene_findMaterial(Scene* scene, const char* name)
{
    return findMaterial(name, &(scene->materials));
}

Mesh** Scene_findMeshes(int* num_meshes, const Scene* scene, const char* name)
{
    return findMeshes(num_meshes, &(scene->meshes), name);
}

int Scene_getNumMaterials(Scene* scene)
{
    return scene->materials.size;
}
