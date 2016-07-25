/*
  Things that can be swapped:
  camera --  perspective vs orthographic
  objects -- implicit surfaces vs explicit surfaces

  Features to consider:
  being able to scale the camera and viewplane 

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
#include "util/vec.h"
#include "util/ray.h"
#include "util/mat.h"
#include "util/constants.h"
#include "gl/glcode.h"
#include "sampling.h"
#include "camera.h"
#include "lights.h"
#include "sceneobj.h"
#include "shading.h"
#include "buildscene.h"
#include "intersect.h"
#include "trace.h"
#include "config.h"

#define SHOW_PROGRESS 1

bool EXIT = false;

vec3 cam_position = {0.0f, 0.0f, 0.0f};

// TODO: 
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
        case GLFW_KEY_Q:
            EXIT = true;
            break;
        }
    }
}


int main()
{
    ConfigParams params;
    parseConfigFile(&params);
    
    GLFWwindow* window = initWindow(params.window_width, params.window_height);
    glfwSetKeyCallback(window, keyCallback);
    
    GlViewport viewport;
    initViewport(&viewport);
    
    unsigned char *image;
    int frame_res_width = params.image_width, frame_res_height = params.image_height;
    int num_pixels = frame_res_width * frame_res_height;
    image = (unsigned char*)calloc(num_pixels * 3, sizeof(char));
    
    // Samples
    setNumSamplesAndSets(params.num_samples, params.num_sample_sets);    // This sets the number of samples and sets for every 
                                                    // sample struct that follows
    srand((unsigned int)time(NULL));    
    Samples2D unit_square_samples = Samples2D_default;
    Samples2D disk_samples = Samples2D_default;
    Samples3D h_samples = Samples3D_default;
    genMultijitteredSamples(&unit_square_samples);
    mapSamplesToDisk(&disk_samples, &unit_square_samples);
    mapSamplesToHemisphere(&h_samples, &disk_samples, 1);

    // Scene data structures
    SceneLights scene_lights;
    SceneObjects scene_objects = SceneObjects_create();
    SceneMaterials scene_materials = SceneMaterials_create();
    SceneMeshes scene_meshes = SceneMeshes_create();
    initScene(&scene_objects, &scene_lights, &scene_materials, &scene_meshes, params.file_name, params.accel_type);

    // Camera
    Camera camera;
    initPinholeCameraDefault(&camera);
    //initThinLensCameraDefault(&camera, DEFAULT_FOCAL_LENGTH, DEFAULT_LENS_RADIUS);
    
    vec3 position = {0.0f, 0.0f, 1.0f};
    vec3 look_point = {0.0f, 0.0f, -3.0f};
    vec3 up_vec = {0.0f, 1.0f, 0.0f};
    cameraLookAt(&camera, position, look_point, up_vec);

    // TODO: throw these in a struct for camera scaling
    float fov = 90.0f;
    float frame_length = 2.0f * (sin(fov/2.0f) * camera.focal_pt_dist);
    float frame_height = frame_length * (float)(frame_res_height)/(float)(frame_res_width);    
    float pixel_length = frame_length/(float)(frame_res_width);

    // Set trace function
    float (*trace)(vec3, int, const vec3, const Ray, const SceneObjects*, const SceneLights*);
    trace = getTraceFunc(params.trace_type);

    time_t startTime, endTime;
    time(&startTime);        

    int prev_percent = 0;
    //drawSamples(image, &disk_samples, frame_res_width, frame_res_height, num_pixels);
    //drawHemisphereSamples2D(image, &h_samples, frame_res_width, frame_res_height, num_pixels);
//#if 0
    for(int i = 0; i < num_pixels; i++)
    {
        vec3 color = {0.0f, 0.0f, 0.0f};
        for(int p = 0; p < NUM_SAMPLES; p++)
        {
            // NOTE: put the code below into a function
            vec2 sample, imageplane_coord;
            getNextSample2D(sample, &unit_square_samples);
            imageplane_coord[0] = -frame_length/2 + pixel_length * ((float)(i % frame_res_width) + sample[0]);
            imageplane_coord[1] = frame_height/2 - pixel_length * ((float)(i / frame_res_width) + sample[1]);
            
            Ray ray;
            calcCameraRay(&ray, imageplane_coord, &camera);

            // Hemisphere sample for ambient occlusion
            vec3 h_sample;
            getNextSample3D(h_sample, &h_samples);

            vec3 radiance;
            trace(radiance, params.max_depth, h_sample, ray, &scene_objects, &scene_lights);
            vec3_add(color, color, radiance);
        }        
        vec3_scale(color, color, 1.0f/NUM_SAMPLES);
        maxToOne(color, color);
        image[i*3] = (char)(color[0] * 255.0f);
        image[i*3 + 1] = (char)(color[1] * 255.0f);
        image[i*3 + 2] = (char)(color[2] * 255.0f);
        
        if(SHOW_PROGRESS)
        {
            int cur_percent = (float)i / (float)(num_pixels) * 100.0f;
            if(cur_percent > prev_percent)
            {
                prev_percent = cur_percent;
                printf("%d%%\n", cur_percent);                
            }
        }
        //displayImage(window, viewport, image, frame_res_width, frame_res_height);        
    }
//#endif 
    time(&endTime);
    double sec = difftime(endTime, startTime);
    printf("%f seconds.\n", sec);

    displayImage(window, viewport, image, frame_res_width, frame_res_height);

    // Clean up
    freeSamples2D(&unit_square_samples);
    freeSamples2D(&disk_samples);
    freeSamples3D(&h_samples);
    freeSceneObjects(&scene_objects);
    freeSceneLights(&scene_lights);
    Camera_destroy(&camera);
    SceneMaterials_destroy(&scene_materials);
    SceneMeshes_destroy(&scene_meshes);

    while(!EXIT)
    {
        glfwPollEvents();
    }
    free(image);    
    return 0;
}
