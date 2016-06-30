#pragma once

#include "../util/util.h"
#include "../util/vec.h"
#include "../util/ray.h"
#include "../util/shaderec.h"
#include "../util/constants.h"
#include "../materials.h"
#include "generic.h"
#include "instanced.h"

typedef InstancedShape OpenCylinder;

OpenCylinder* initOpenCylinder(const mat4 inv_transform, const float phi, Material* mat,
                               const NormalType normal_type, const bool shadow)
{
    GenericOpenCylinder* goc = (GenericOpenCylinder*)malloc(sizeof(GenericOpenCylinder));
    goc->shadow = shadow;
    goc->radius = 1.0f;
    goc->half_height = 1.0f;
    goc->phi = phi;
    goc->mat = mat;
    goc->normal_type = normal_type;
    
    OpenCylinder* open_cyl = (OpenCylinder*)malloc(sizeof(OpenCylinder));
    open_cyl->obj.ptr = goc;
    open_cyl->obj.type = GENERICOPENCYLINDER;
    open_cyl->mat = mat;
    mat4_copy(open_cyl->inv_transform, inv_transform);
    return open_cyl;
}

typedef CompoundObject CompoundSolidCylinder;

void calcAABBSolidCylinder(AABB* aabb, const Disk* top, const Disk* bottom, const GenericOpenCylinder* tube)
{
    aabb->min[0] = -(tube->radius);
    aabb->min[1] = -(tube->half_height);
    aabb->min[2] = -(tube->radius);

    aabb->max[0] = tube->radius;
    aabb->max[1] = tube->half_height;
    aabb->max[2] = (tube->radius);    
}

CompoundSolidCylinder* initCompoundSolidCylinder(const float radius, const float half_height, const float phi,
                    Material* mat, const bool shadow)
{
    int num_obj = 3;
    CompoundSolidCylinder* sc = (CompoundSolidCylinder*)malloc(sizeof(CompoundSolidCylinder));
    sc->objects = (Object_t*)malloc(sizeof(Object_t) * num_obj);
    sc->shadow = shadow;
    Disk* top = (Disk*)malloc(sizeof(Disk));
    Disk* bottom = (Disk*)malloc(sizeof(Disk));
    GenericOpenCylinder* tube = (GenericOpenCylinder*)malloc(sizeof(GenericOpenCylinder));
    top->radius = bottom->radius = tube->radius = radius;
    tube->half_height = half_height;
    vec3_assign(top->center, 0.0f, half_height, 0.0f);
    vec3_assign(bottom->center, 0.0f, -half_height, 0.0f);
    vec3_assign(top->normal, 0.0f, 1.0f, 0.0f);
    vec3_assign(bottom->normal, 0.0f, -1.0f, 0.0f);
    tube->phi = phi;
    tube->normal_type = OPEN;
    top->mat = bottom->mat = tube->mat = mat;

    sc->objects[0].ptr = top;
    sc->objects[0].type = DISK;
    sc->objects[1].ptr = bottom;
    sc->objects[1].type = DISK;
    sc->objects[2].ptr = tube;
    sc->objects[2].type = GENERICOPENCYLINDER;
    sc->num_obj = num_obj;
    calcAABBSolidCylinder(&(sc->aabb), top, bottom, tube);
    return sc;
}

typedef InstancedShape SolidCylinder;

SolidCylinder* initSolidCylinder(const mat4 inv_transform, const float radius, const float half_height,
                                 const float phi, Material* mat, const bool shadow)
{
    CompoundSolidCylinder* csc = initCompoundSolidCylinder(radius, half_height, phi, mat, shadow);
    SolidCylinder* sc = (SolidCylinder*)malloc(sizeof(SolidCylinder));
    sc->obj.ptr = csc;
    sc->obj.type = COMPOUND;
    sc->mat = mat;
    mat4_copy(sc->inv_transform, inv_transform);
    return sc;
}
