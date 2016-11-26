#pragma once

#include "sampling.h"
#include "util/vec.h"
#include "util/constants.h"
#include "util/ray.h"

enum CameraType
{
    Pinhole,
    ThinLens
};

typedef struct Camera
{
    CameraType camera_type;    
    float focal_pt_dist;                // Distance between focal_point and the view plane
    float focal_length;                 // Distance between the view plane and focal plane
    float lens_radius;    
    vec3 position;
    vec3 x_axis;
    vec3 y_axis;
    vec3 z_axis;
    vec3 focal_point;
    Samples2D* samples;
} Camera;

void Camera_destroy(Camera* camera)
{
    freeSamples2D(camera->samples);
    free(camera->samples);
    camera->focal_pt_dist = 0.0f;
    camera->focal_length = 0.0f;
    camera->lens_radius = 0.0f;
    vec3_copy(camera->position, ORIGIN);
    vec3_copy(camera->x_axis, ORIGIN);
    vec3_copy(camera->y_axis, ORIGIN);
    vec3_copy(camera->z_axis, ORIGIN);
    vec3_copy(camera->focal_point, ORIGIN);        
}

void initPinholeCameraDefault(Camera *camera)
{
    vec3 position = {0.0f, 0.0f, 0.0f};
    vec3 x_axis = {1.0f, 0.0f, 0.0f};
    vec3 y_axis = {0.0f, 1.0f, 0.0f};
    vec3 z_axis = {0.0f, 0.0f, 1.0f};
    vec3_copy(camera->position, position);
    vec3_copy(camera->x_axis, x_axis);
    vec3_copy(camera->y_axis, y_axis);
    vec3_copy(camera->z_axis, z_axis);
    camera->focal_pt_dist = 0.035f;
    vec3 focal_point = {0.0f, 0.0f, camera->focal_pt_dist};
    vec3_copy(camera->focal_point, focal_point);
    camera->camera_type = Pinhole;
    camera->focal_length = 0.0f;
    camera->lens_radius = 0.0f;
    // TODO: move below into initThinLensCameraDefault
    camera->samples = (Samples2D*)malloc(sizeof(Samples2D));
    *(camera->samples) = getDefaultSamples2D();
    genMultijitteredSamples(camera->samples);
    mapSamplesToDisk(camera->samples, camera->samples);
}

void initThinLensCameraDefault(Camera *camera, const float focal_length, const float lens_radius, const int num_samples,
    const int num_sets)
{
    initPinholeCameraDefault(camera);
    camera->camera_type = ThinLens;
    camera->focal_length = focal_length;
    camera->lens_radius = lens_radius;
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
    vec3_scale(displacement, camera->z_axis, camera->focal_pt_dist);
    vec3_add(camera->focal_point, position, displacement);
    vec3_copy(camera->position, position);
}

// TODO: change the parameter to fov
void setFocalDist(Camera *camera, const float focal_pt_dist)
{
    camera->focal_pt_dist = focal_pt_dist;
    vec3 displacement;
    vec3_scale(displacement, camera->z_axis, focal_pt_dist);
    vec3_add(camera->focal_point, camera->position, displacement);
}

// Compute a ray for pinhole camera
void calcRayPinhole(Ray *ray, const vec2 viewplane_coord, const Camera *camera)
{
    vec3 sample_loc_cam, tmp_vec_1, tmp_vec_2, focal_point;
    vec3_scale(tmp_vec_1, camera->x_axis, viewplane_coord[0]);
    vec3_scale(tmp_vec_2, camera->y_axis, viewplane_coord[1]);
    vec3_add(sample_loc_cam, tmp_vec_1, tmp_vec_2);
    vec3_scale(focal_point, camera->z_axis, camera->focal_pt_dist);    // focal_point is rotated but not translated
    vec3_sub(ray->direction, sample_loc_cam, focal_point);
    vec3_normalize(ray->direction, ray->direction);
    vec3_add(ray->origin, focal_point, camera->position);
}

// Compute a ray for thin lens camera
void calcRayThinLens(Ray *ray, const vec2 viewplane_coord, const Camera *camera, const int sample_index)
{
    vec2 disk_sample;
    //getNextSample2D(disk_sample, camera->samples);
    getSample2D(disk_sample, camera->samples, sample_index);
    // NOTE: consider doing the calculations in camera space instead of transforming the resulting ray to camera space
    vec3 point_viewplane = {viewplane_coord[0], viewplane_coord[1], 0.0f};
    vec3_copy(ray->origin, camera->focal_point);
    vec3_sub(ray->direction, point_viewplane, ray->origin);
    vec3_normalize(ray->direction, ray->direction);
    vec3 focal_plane_point = {point_viewplane[0] * (camera->focal_length/camera->focal_pt_dist),
                              point_viewplane[1] * (camera->focal_length/camera->focal_pt_dist),
                              -camera->focal_length};
    vec3_assign(ray->origin, disk_sample[0] * camera->lens_radius,
                disk_sample[1] * camera->lens_radius, camera->focal_pt_dist);
    vec3_sub(ray->direction, focal_plane_point, ray->origin);
    vec3_normalize(ray->direction, ray->direction);

    // Transform ray to camera space
    // NOTE: shrink below into two function calls
    vec3 tmp_vec_1, tmp_vec_2, tmp_vec_3;
    vec3_scale(tmp_vec_1, camera->x_axis, ray->direction[0]);
    vec3_scale(tmp_vec_2, camera->y_axis, ray->direction[1]);
    vec3_scale(tmp_vec_3, camera->z_axis, ray->direction[2]);
    vec3_add(tmp_vec_1, tmp_vec_1, tmp_vec_2);
    vec3_add(ray->direction, tmp_vec_1, tmp_vec_3);
    vec3_scale(tmp_vec_1, camera->x_axis, ray->origin[0]);
    vec3_scale(tmp_vec_2, camera->y_axis, ray->origin[1]);
    vec3_scale(tmp_vec_3, camera->z_axis, ray->origin[2]);
    vec3_add(tmp_vec_1, tmp_vec_1, tmp_vec_2);
    vec3_add(tmp_vec_1, tmp_vec_1, tmp_vec_3);
    vec3_add(ray->origin, tmp_vec_1, camera->position);
}

void calcCameraRay(Ray* ray, const vec2 imageplane_coord, const Camera* camera, const int sample_index)
{
    switch(camera->camera_type)
    {
    case Pinhole:
    {
        calcRayPinhole(ray, imageplane_coord, camera);
    } break;
    case ThinLens:
    {
        calcRayThinLens(ray, imageplane_coord, camera, sample_index);
    } break;
    }
}
