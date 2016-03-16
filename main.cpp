/*
  Things that can be swapped:
  camera --  perspective vs orthographic
  objects -- implicit surfaces vs explicit surfaces
  sampling -- different sampling algorithms

  TODO:
  gamma correction
  consider how variables should be grouped into objects?
 */

#include <GL/glew.h>

#if __GNUG__
#    include <GLFW/glfw3.h>
#else
#    include <GL/glfw3.h>
#    include <SOIL.h>
#endif

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include "vec.h"
#include "shaders.h"
#include "shapes.h"
#include "glcode.h"
#include "sampling.h"
#include "camera.h"

static const int MAX_SPHERES = 1000;

vec3 cam_position = {0.0f, 0.0f, 0.0f};

// TODO: not finished
void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    if(action == GLFW_PRESS)
    {
        switch(key)
        {
        case GLFW_KEY_A:
            cam_position[0] -= 0.01f;
            break;
        case GLFW_KEY_D:
            cam_position[0] += 0.01f;
            break;
        case GLFW_KEY_W:
            cam_position[2] -= 0.01f;
            break;
        case GLFW_KEY_S:
            cam_position[2] += 0.01f;
            break;            
        }
    }
}

int initSpheres(Sphere spheres[])
{
    int sphere_count = 0;
    vec3_assign(spheres[sphere_count].center, 0.0f, 0.0f, -4.0f);
    spheres[sphere_count].radius = 1.0f;
    vec3_assign(spheres[sphere_count].color, 255.0f, 0.0f, 0.0f);
    sphere_count++;

    vec3_assign(spheres[sphere_count].center, 0.75f, 0.0f, -4.0f);
    spheres[sphere_count].radius = 1.0f;
    vec3_assign(spheres[sphere_count].color, 0.0f, 255.0f, 0.0f);
    sphere_count++;

    return sphere_count;
}

void transformPoint(vec3 r, const vec3 a, const vec3 x, const vec3 y, const vec3 z, const vec3 cam_position)
{
    vec3 result_vec;
    result_vec[0] = vec3_dot(a, x[0], y[0], z[0]);
    result_vec[1] = vec3_dot(a, x[1], y[1], z[1]);
    result_vec[2] = vec3_dot(a, x[2], y[2], z[2]);    
    vec3_add(result_vec, result_vec, cam_position);
    vec3_copy(r, result_vec);
}

int main()
{
    int window_width = 1920, window_height = 1080;
    GLFWwindow* window = initWindow(1920, 1080);
    glfwSetKeyCallback(window, keyCallback);
    
    GlViewport viewport;
    initViewport(&viewport);
    
    unsigned char *image;
    int frame_res_width = 720, frame_res_height = 405;
    int num_pixels = frame_res_width * frame_res_height;
    image = (unsigned char*)calloc(num_pixels * 3, sizeof(char));

    /* TODO: Moving the camera
       Should be done after sphere is implemented so that we can see the result
       Two Vec3's. One for translation and one for euler angle. Done?
       First calculate ray direction, then rotate it with euler angle.
       Rotate and translate the position of the pixel to use as ray starting point.
       Implement euler-angle-to-matrix conversion.
     */

    /*
      Camera related things: fov, resolution, position, orientation.
      Be cool and throw them into a class?
      What design pattern is this?
     */
    /*
      Array of objects
      Array of lights?
     */
    Sphere spheres[MAX_SPHERES];
    int sphere_count = 0;
    sphere_count = initSpheres(spheres);
    
    // TODO: hammersley sampling doesn't seem to work horizontally
    int num_samples = 25;
    int num_sets = 3;    
    Samples samples = Samples_default;
    //genRegularSamples(&samples, num_samples, num_sets);
    //genMultijitteredSamples(&samples, num_samples, num_sets);
    genHammersleySamples(&samples, num_samples, num_sets);
    mapSamplesToDisk(&samples);
    //mapSamplesToHemisphere(&samples, 5);

    Camera camera;
    initThinLensCameraDefault(&camera);

    vec3 look_point = {0.75f, 0.0f, -4.0f};
    vec3 up_vec = {0.0f, 1.0f, 0.0f};
    cameraLookAt(&camera, cam_position, look_point, up_vec);

    float fov = 90.0f;
    float frame_length = 2.0f * (sin(fov/2.0f) * camera.focal_pt_dist);
    float frame_height = frame_length * (float)(frame_res_height)/(float)(frame_res_width);    
    float pixel_length = frame_length/(float)(frame_res_width);

    float focal_length = 4.0f;
    float lens_radius = 1.0f;
    for(int i = 0; i < num_pixels; i++)
    {
        vec3 color = {0.0f, 0.0f, 0.0f};
        for(int p = 0; p < num_samples; p++)
        {
            vec2 sample, disk_sample, viewplane_coord;
            getNextSample(sample, &samples);
            viewplane_coord[0] = -frame_length/2 + pixel_length * ((float)(i % frame_res_width) + sample[0]);            
            viewplane_coord[1] = frame_height/2 - pixel_length * ((float)(i / frame_res_width) + sample[1]);            
            Ray ray;
            switch(camera.camera_type)
            {
            case Pinhole:
                calcRayCam(&ray, viewplane_coord, &camera);
                break;
            case ThinLens:
                getNextDiskSample(disk_sample, &samples);            
                calcRayThinLens(&ray, viewplane_coord, disk_sample, focal_length, lens_radius, &camera);
                break;
            }
            
            float min_t = TMAX,  tmp_t = TMAX;
            ShadeRec min_sr, tmp_sr;
            for(int i = 0; i < sphere_count; i++)
            {
                tmp_t = rayIntersectSphere(&tmp_sr, spheres[i], ray);
                if(tmp_t < min_t)
                {
                    min_t = tmp_t;
                    min_sr = tmp_sr;
                }
            }
            if(min_t < TMAX)
            {
                vec3_add(color, color, min_sr.color);
            }
        }
        vec3_scale(color, color, 1.0f/num_samples);
        image[i*3] = (char)(color[0]);
        image[i*3 + 1] = (char)(color[1]);
        image[i*3 + 2] = (char)(color[2]);
    }
    displayImage(window, viewport, image, frame_res_width, frame_res_height);
    free(image);
    freeSamples(&samples);

    while(true)
    {

    }
    return 0;
}
