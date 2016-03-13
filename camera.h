#pragma once

#include "vec.h"

typedef struct Camera
{
    vec3 position;
    vec3 x_axis;
    vec3 y_axis;
    vec3 z_axis;
    vec3 focal_point;
    float focal_dist;
} Camera;

void initCameraDefault(Camera *camera)
{
    vec3 position = {0.0f, 0.0f, 0.0f};
    vec3 x_axis = {1.0f, 0.0f, 0.0f};
    vec3 y_axis = {0.0f, 1.0f, 0.0f};
    vec3 z_axis = {0.0f, 0.0f, 1.0f};
    vec3_copy(camera->position, position);
    vec3_copy(camera->x_axis, x_axis);
    vec3_copy(camera->y_axis, y_axis);
    vec3_copy(camera->z_axis, z_axis);
    camera->focal_dist = 1.0f;
    vec3 focal_point = {0.0f, 0.0f, 1.0f};
}

// Sets up orthonormal basis, position, and focal point
void cameraLookAt(Camera *camera, const vec3 position, const vec3 look_point, const vec3 up_vec)
{
    vec3 cam_upright_look_point;
    vec3_sub(cam_upright_look_point, look_point, position);
    vec3 look_vec;
    vec3_normalize(look_vec, cam_upright_look_point);
    vec3_negate(camera->z_axis, look_vec);
    vec3_cross(camera->x_axis, camera->z_axis, up_vec);
    vec3_normalize(camera->x_axis, camera->x_axis);
    vec3_negate(camera->x_axis, camera->x_axis);
    vec3_cross(camera->y_axis, camera->z_axis, camera->x_axis);
    vec3 displacement;
    vec3_scale(displacement, camera->z_axis, camera->focal_dist);
    vec3_add(camera->focal_point, position, displacement);
    vec3_copy(camera->position, position);
}
