#include "generic.h"
#include "torus.h"

Torus* initTorus(const mat4 inv_transform, const float swept_radius, const float tube_radius, const float phi, 
                 Material* mat, const bool shadow)
{
    GenericTorus* gt = (GenericTorus*)malloc(sizeof(GenericTorus));
    gt->shadow = shadow;
    gt->swept_radius = swept_radius;
    gt->tube_radius = tube_radius;
    gt->phi = phi;
    gt->mat = mat;

    Torus* torus = (Torus*)malloc(sizeof(Torus));
    torus->obj.ptr = gt;
    torus->obj.type = GENERICTORUS;
    torus->mat = mat;
    mat4_copy(torus->inv_transform, inv_transform);
    return torus;
}
