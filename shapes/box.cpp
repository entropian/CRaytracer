#include "box.h"
#include "generic.h"

Box* initBox(const mat4 inv_transform, const float x_span, const float y_span, const float z_span, Material* mat, const bool shadow)
{
    float half_x = x_span / 2.0f;
    float half_y = y_span / 2.0f;
    float half_z = z_span / 2.0f;    
    AABox* aabox = (AABox*)malloc(sizeof(AABox));
    vec3_assign(aabox->min, -half_x, -half_y, -half_z);
    vec3_assign(aabox->max, half_x, half_y, half_z);
    aabox->mat = mat;

    Box* box = (Box*)malloc(sizeof(Box));
    box->obj.ptr = aabox;
    box->obj.type = AABOX;
    box->mat = mat;
    mat4_copy(box->inv_transform, inv_transform);
    return box;
}
