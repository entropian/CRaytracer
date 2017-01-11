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

void Camera_destroy(Camera* camera);
void initPinholeCameraDefault(Camera *camera);
void initThinLensCameraDefault(Camera *camera, const float focal_length, const float lens_radius, const int num_samples,
                               const int num_sets);
// Sets up orthonormal basis, position, and focal point
void cameraLookAt(Camera *camera, const vec3 position, const vec3 look_point, const vec3 up_vec);
// TODO: change the parameter to fov
void setFocalDist(Camera *camera, const float focal_pt_dist);

// Compute a ray for pinhole camera
void calcRayPinhole(Ray *ray, const vec2 viewplane_coord, const Camera *camera);
// Compute a ray for thin lens camera
void calcRayThinLens(Ray *ray, const vec2 viewplane_coord, const Camera *camera, const int sample_index);
void calcCameraRay(Ray* ray, const vec2 imageplane_coord, const Camera* camera, const int sample_index);

typedef struct Film_s
{
    float fov;
    float frame_length, frame_height; // Dimension of the film inside the scene
    int frame_res_width, frame_res_height; // Actual image resolution
    float pixel_length;
    int num_pixels;
    Samples2D samples;
}Film;

void calcFilmDimension(Film *film, const Camera *camera);
void calcImageCoord(vec2 image_coord, Film *film, const int sample_index, const int pixel_index);

