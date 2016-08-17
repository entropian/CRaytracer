#pragma once

#include "scenedata.h"

typedef struct Scene_s
{
    SceneObjects objects;
    SceneMaterials materials;
    SceneTextures textures;
    SceneMeshes meshes;
    SceneLights lights;
}Scene;

Scene Scene_create()
{
    Scene scene;
    scene.objects = SceneObjects_create();
    scene.materials = SceneMaterials_create();
    scene.textures = SceneTextures_create();
    scene.meshes = SceneMeshes_create();
    return scene;
}

void Scene_destroy(Scene* scene)
{
    freeSceneObjects(&(scene->objects));
    SceneMaterials_destroy(&(scene->materials));
    SceneTextures_destroy(&(scene->textures));
    SceneMeshes_destroy(&(scene->meshes));    
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

OBJShape* Scene_addMesh(Scene* scene, const OBJShape* obj_shape)
{
    return SceneMeshes_push(&(scene->meshes), obj_shape);
}

Material* Scene_findMaterial(Scene* scene, const char* name)
{
    return findMaterial(name, &(scene->materials));
}

OBJShape** Scene_findMeshes(int* num_meshes, const Scene* scene, const char* name)
{
    return findMeshes(num_meshes, &(scene->meshes), name);
}

int Scene_getNumMaterials(Scene* scene)
{
    return scene->materials.size;
}
