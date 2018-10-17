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
typedef CompoundObject CompoundSolidCylinder;
typedef InstancedShape SolidCylinder;
OpenCylinder* initOpenCylinder(const mat4 inv_transform, const float phi, Material* mat,
                               const NormalType normal_type);
void calcAABBSolidCylinder(AABB* aabb, const Disk* top, const Disk* bottom, const GenericOpenCylinder* tube);
CompoundSolidCylinder* initCompoundSolidCylinder(const float radius, const float half_height, const float phi,
                                                 Material* mat);
SolidCylinder* initSolidCylinder(const mat4 inv_transform, const float radius, const float half_height,
                                 const float phi, Material* mat);

