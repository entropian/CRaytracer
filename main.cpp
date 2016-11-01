
#include <GL/glew.h>
#include <GLFW/glfw3.h>

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
#include "scene/scene.h"
#include "shading.h"
#include "buildscene.h"
#include "intersect.h"
#include "trace.h"
#include "config.h"
#include "texture.h"
#include "noise.h"
#include "imagefile.h"
#include "photonmap.h"
#include "projmap.h"

#define SHOW_PROGRESS 1
#define PROG
#define CORNELL_BOX

extern double g_traversal_time;

bool EXIT = false;
int MAX_DEPTH = 0;

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
    MAX_DEPTH = params.max_depth;
    
    GLFWwindow* window = initWindow(params.window_width, params.window_height);
    glfwSetKeyCallback(window, keyCallback);
    
    GlViewport viewport;
    initViewport(&viewport);
    
    unsigned char *image;
    int frame_res_width = params.image_width, frame_res_height = params.image_height;
    int num_pixels = frame_res_width * frame_res_height;
    image = (unsigned char*)calloc(num_pixels * 3, sizeof(char));

    float* color_buffer = (float*)calloc(num_pixels * 3, sizeof(float));

    //LatticeNoise_init(CUBIC, 5, 1.0f, 2.0f);
    
    // Samples
    setNumSamplesAndSets(params.num_samples, params.num_sample_sets);    // This sets the number of samples and sets for every 
                                                                         // sample struct that follows

#ifdef PROG
    setInterleaved(true);
    unsigned char* set_buffer = (unsigned char*)malloc(sizeof(unsigned char) * num_pixels);
    for(unsigned int i = 0; i < num_pixels; i++)
    {
        set_buffer[i] = (unsigned char)(rand() % params.num_sample_sets);
    }    
#endif
    srand((unsigned int)time(NULL));    
    Samples2D unit_square_samples = getDefaultSamples2D();
    Samples2D disk_samples = getDefaultSamples2D();
    Samples3D h_samples = getDefaultSamples3D();
    genMultijitteredSamples(&unit_square_samples);
    mapSamplesToDisk(&disk_samples, &unit_square_samples);
    mapSamplesToHemisphere(&h_samples, &disk_samples, 1);

    // Scene data structures
    Scene scene = Scene_create();
    initScene(&scene, params.file_name, params.accel_type);
    vec3 bsphere_centers[MAX_BOUNDING_SPHERES];
    float bsphere_radii[MAX_BOUNDING_SPHERES];
    // NOTE: calcCausticBoundingSpheres before building accel
    int num_bsphere = calcCausticBoundingSpheres(bsphere_centers, bsphere_radii, &(scene.objects));
    buildSceneAccel(&scene);

    buildProjMap(&(scene.lights), bsphere_centers, bsphere_radii, num_bsphere);

    SceneLights *sl = &(scene.lights);
    PointLight *point_light = NULL;    
    for(int i = 0; i < sl->num_lights; i++)
    {
        if(sl->light_types[i] == POINTLIGHT)
        {
            point_light = (PointLight*)(sl->light_ptrs[i]);
        }
    }

    if(point_light)
    {
        for(int i = 0; i < THETA_ROW; i++)
        {
            for(int j = 0; j < PHI_COLUMN; j++)
            { 
                printf("%d ", point_light->proj_map[i*PHI_COLUMN + j]);
            }
            printf("\n");
        }
    }

    // Photon map

    const int num_photons = params.pm_config.num_photons;
    const int max_bounce = params.pm_config.photon_depth;
    const int num_caustic_photons = params.pm_config.num_caustic_photons;
    Photonmap photon_map, caustic_map;    
    bool photon_map_status = false;
    if(params.photon_map && params.trace_type == WHITTED)
    {
        photon_map_status = true;
        Photonmap_init(&photon_map, num_photons, max_bounce);
        Photonmap_init(&caustic_map, num_caustic_photons, max_bounce);
        emitPhotons(&photon_map, &(scene.objects), &(scene.lights));
        emitCaustics(&caustic_map, &(scene.objects), &(scene.lights));
        Photonmap_balance(&photon_map);
        Photonmap_balance(&caustic_map);            
    }
    PhotonQueryVars query_vars;    
    query_vars.nphotons = params.pm_config.num_estimate;
    query_vars.photon_radius = params.pm_config.photon_radius;
    query_vars.caustic_radius = params.pm_config.caustic_radius;    

    // Camera
    Camera camera;
    initPinholeCameraDefault(&camera);
    //initThinLensCameraDefault(&camera, DEFAULT_FOCAL_LENGTH, DEFAULT_LENS_RADIUS);

#ifdef CORNELL_BOX    
    vec3 position = {278.0f, 273.0f, 800.0f};
    vec3 look_point = {278.0f, 273.0f, 0.0f};
#else
    vec3 position = {0.0f, 2.0f, 5.0f};
    vec3 look_point = {0.0f, 0.0f, 0.0f};
#endif
    vec3 up_vec = {0.0f, 1.0f, 0.0f};
    cameraLookAt(&camera, position, look_point, up_vec);

    // TODO: throw these in a struct for camera scaling
    float fov = 70.0f / 180.0f * PI;
    float frame_length = 2.0f * (sin(fov/2.0f) * camera.focal_pt_dist);
    float frame_height = frame_length * (float)(frame_res_height)/(float)(frame_res_width);    
    float pixel_length = frame_length/(float)(frame_res_width);

    // Set trace function
    float (*trace)(vec3, int, const vec3, const Ray, const SceneObjects*, const SceneLights*, const int);
    trace = getTraceFunc(params.trace_type);

    double start_time, end_time;
    start_time = glfwGetTime(); 
    int prev_percent = 0;
//#if 0
#ifdef PROG
    for(unsigned int p = 0; p < params.num_samples; p++)
    {
        for(int i = 0; i < num_pixels; i++)
        {
            int sample_index = calcInterleavedSampleIndex(p, set_buffer[i]);
            vec3 color = {0.0f, 0.0f, 0.0f};
            vec2 sample, imageplane_coord;
            getSample2D(sample, &unit_square_samples, sample_index);
            imageplane_coord[0] = -frame_length/2 + pixel_length * ((float)(i % frame_res_width) + sample[0]);
            imageplane_coord[1] = frame_height/2 - pixel_length * ((float)(i / frame_res_width) + sample[1]);
            
            Ray ray;
            calcCameraRay(&ray, imageplane_coord, &camera, sample_index);

            // Hemisphere sample for ambient occlusion
            vec3 h_sample;
            getSample3D(h_sample, &h_samples, sample_index);            

            vec3 radiance;
            trace(radiance, params.max_depth, h_sample, ray, &(scene.objects), &(scene.lights), sample_index);
            vec3_add(color, color, radiance);

            // Photon map
            if(photon_map_status)
            {
				vec3 pm_color = {0.0f, 0.0f, 0.0f};
                calcPhotonmapComponent(pm_color, h_sample, query_vars, &photon_map, &caustic_map, &(scene.objects), ray);
                vec3_add(color, color, pm_color);
            }

            // NEW
            color_buffer[i*3] += color[0];
            color_buffer[i*3 + 1] += color[1];
            color_buffer[i*3 + 2] += color[2];
            vec3_assign(color, color_buffer[i*3], color_buffer[i*3 +1], color_buffer[i*3 + 2]);
            vec3_scale(color, color, 1/(float)(p+1));
            maxToOne(color, color);            

            image[i*3] = (char)(color[0] * 255.0f);
            image[i*3 + 1] = (char)(color[1] * 255.0f);
            image[i*3 + 2] = (char)(color[2] * 255.0f);            
        }        
        if(SHOW_PROGRESS)
        {
            displayImage(window, viewport, image, frame_res_width, frame_res_height);                
            //int cur_percent = (int)((float)i / (float)(num_pixels) * 100.0f);
            int cur_percent = (int)((float)p / (float)(params.num_samples) * 100.0f);
            if(cur_percent > prev_percent)
            {
                prev_percent = cur_percent;
                printf("%d%%\n", cur_percent);

            }
        }
    }

#else
    for(int i = 0; i < num_pixels; i++)
    {
        vec3 color = {0.0f, 0.0f, 0.0f};
        for(unsigned int p = 0; p < params.num_samples; p++)
        {
            int sample_index = calcNextSampleIndex();
            vec2 sample, imageplane_coord;
            //getNextSample2D(sample, &unit_square_samples);
            getSample2D(sample, &unit_square_samples, sample_index);
            imageplane_coord[0] = -frame_length/2 + pixel_length * ((float)(i % frame_res_width) + sample[0]);
            imageplane_coord[1] = frame_height/2 - pixel_length * ((float)(i / frame_res_width) + sample[1]);
            
            Ray ray;
            calcCameraRay(&ray, imageplane_coord, &camera, sample_index);

            // Hemisphere sample for ambient occlusion
            vec3 h_sample;
            getSample3D(h_sample, &h_samples, sample_index);

            vec3 radiance;
            trace(radiance, params.max_depth, h_sample, ray, &(scene.objects), &(scene.lights), sample_index);
            vec3_add(color, color, radiance);
        }        
        vec3_scale(color, color, 1.0f/params.num_samples);
        maxToOne(color, color);
        image[i*3] = (char)(color[0] * 255.0f);
        image[i*3 + 1] = (char)(color[1] * 255.0f);
        image[i*3 + 2] = (char)(color[2] * 255.0f);

        if(SHOW_PROGRESS)
        {
            int cur_percent = (int)((float)i / (float)(num_pixels) * 100.0f);
            if(cur_percent > prev_percent)
            {
                prev_percent = cur_percent;
                printf("%d%%\n", cur_percent);
                displayImage(window, viewport, image, frame_res_width, frame_res_height);                
            }
        }
    }
#endif
//#endif 
    end_time = glfwGetTime();
    double sec = end_time - start_time;
    printf("%f seconds.\n", sec);

    displayImage(window, viewport, image, frame_res_width, frame_res_height);
    printf("Traversal time = %f\n", g_traversal_time);

    if(params.image_save)
    {
        PPM_write("output.ppm", image, num_pixels * 3, params.image_width, params.image_height);
        EXIT = true;
    }

    // Clean up
    freeSamples2D(&unit_square_samples);
    freeSamples2D(&disk_samples);
    freeSamples3D(&h_samples);
    Scene_destroy(&scene);
    Camera_destroy(&camera);
    Photonmap_destroy(&photon_map);
    Photonmap_destroy(&caustic_map);
    
    double frames_per_sec = 10.0;
    double time_between_frames = 1.0 / frames_per_sec;
    double current_time, last_draw_time = 0.0;
    while(!EXIT)
    {
        current_time = glfwGetTime();
        if((current_time - last_draw_time) >= time_between_frames)
        {
            displayImage(window, viewport, image, frame_res_width, frame_res_height);
            last_draw_time = current_time;
        }
        glfwPollEvents();
    }
    free(image);    
    return 0;
}
