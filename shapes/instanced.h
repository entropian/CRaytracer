#pragma once

#include "../util/mat.h"
#include "objecttype.h"
#include "generic.h"

typedef struct CompoundObject
{
    bool shadow;
    int num_obj;
    Object_t* objects;
    AABB aabb;    
} CompoundObject;

typedef struct InstancedShape
{
    Object_t obj;
    mat4 inv_transform;
    Material* mat;        
} InstancedShape;

float rayIntersectCompound(ShadeRec* sr, CompoundObject* co, const Ray ray);
float shadowRayIntersectCompound(CompoundObject* co, const Ray shadow_ray);
void freeInstancedShape(InstancedShape* is);

void freeCompoundObject(CompoundObject* co);
void freeInstancedShape(InstancedShape* is);

float rayIntersectInstanced(ShadeRec* sr, InstancedShape* is, const Ray src_ray);
float shadowRayIntersectInstanced(InstancedShape* is, const Ray src_ray);

float rayIntersectCompound(ShadeRec* sr, CompoundObject* co, const Ray ray);
float shadowRayIntersectCompound(CompoundObject* co, const Ray shadow_ray);
