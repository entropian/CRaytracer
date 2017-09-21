#pragma once

#include <cstdlib>
#include "../util/util.h"
#include "../util/vec.h"
#include "../util/ray.h"
#include "../util/math.h"
#include "../util/shaderec.h"
#include "../util/constants.h"
#include "../materials.h"
#include "../aabb.h"
#include "instanced.h"

typedef InstancedShape Torus;

Torus* initTorus(const mat4 inv_transform, const float swept_radius, const float tube_radius, const float phi, 
                 Material* mat);
