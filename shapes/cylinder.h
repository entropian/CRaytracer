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

