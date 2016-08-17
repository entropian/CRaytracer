#pragma once

#include "scenedata.h"

typedef struct Scene_s
{
    SceneObjects objects;
    SceneMaterials materials;
    SceneMeshes meshes;
    SceneLights lights;
}Scene;

Scene Scene_create()
{
    Scene scene;
    scene.objects = SceneObjects_create();
    scene.materials = SceneMaterials_create();
    scene.meshes = SceneMeshes_create();
    return scene;
}

void Scene_addObject(Scene* scene, const Object_t* obj)
{
    SceneObjects_push_obj(&(scene->objects), obj);
}

Material* Scene_addMaterial(Scene* scene, const Material* mat, const char* name)
{
    return SceneMaterials_push(&(scene->materials), mat, name);
}

OBJShape* Scene_addMesh(Scene* scene, const OBJShape* obj_shape)
{
    return SceneMeshes_push(&(scene->meshes), obj_shape);
}

Material* Scene_findMaterial(Scene* scene, const char* name)
{
    return findMaterial(name, &(scene->materials));
}

void Scene_destroy(Scene* scene)
{
    freeSceneObjects(&(scene->objects));
    SceneMaterials_destroy(&(scene->materials));
    SceneMeshes_destroy(&(scene->meshes));
    freeSceneLights(&(scene->lights));    
}
